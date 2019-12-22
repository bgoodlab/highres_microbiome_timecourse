import sys
import pylab
import numpy
import parse_midas_data
import parse_timecourse_data
import stats_utils
import barcode_utils
import gzip

import os.path
import config

import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint
import matplotlib.colors as mcolors

mpl.rcParams['font.size'] = 7
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

readcloud_bins = numpy.logspace(1,4,30)
#readcloud_bins = numpy.logspace(1,4,60)
readcloud_centers = readcloud_bins[1:]


####################################################
#
# Set up Figure (3 panels, arranged in 1x3 grid)
#
####################################################



pangenome_species = parse_midas_data.parse_pangenome_species()
#pangenome_species = ['Bacteroides_vulgatus_57955']

# long gene map
# long gene = (species, gene) tuple
# id = int
longgene_id_map = {}
id_longgene_map = []

all_samples = parse_timecourse_data.morteza_samples    
desired_samples = ['6037','6037.2','6037.3']
colors = ['#253494', '#41b6c4', '#a1dab4']
sample_time_map = parse_timecourse_data.parse_sample_time_map()

total_barcode_map = barcode_utils.parse_total_barcode_map()

pylab.figure(1,figsize=(3.42,3))

fig = pylab.gcf()

outer_grid  = gridspec.GridSpec(3,1, height_ratios=[0.05,0.5,1], hspace=0.1)

legend_axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(legend_axis)

legend_axis.set_ylim([0,1])
legend_axis.set_xlim([0,1])

legend_axis.spines['top'].set_visible(False)
legend_axis.spines['right'].set_visible(False)
legend_axis.spines['left'].set_visible(False)
legend_axis.spines['bottom'].set_visible(False)

legend_axis.set_xticks([])
legend_axis.set_yticks([])
   

density_axis = plt.Subplot(fig, outer_grid[1])
fig.add_subplot(density_axis)

density_axis.set_ylabel('Total reads\n contributed')
density_axis.semilogx([1e-01],[1e-01])
density_axis.set_xlim([1e01,1e04])
    
density_axis.spines['top'].set_visible(False)
density_axis.spines['right'].set_visible(False)
density_axis.get_xaxis().tick_bottom()
density_axis.get_yaxis().tick_left()

species_axis = plt.Subplot(fig, outer_grid[2])
fig.add_subplot(species_axis)

species_axis.set_xlabel('Read pairs per barcode')

species_axis.set_ylabel('Effective # species')
species_axis.set_ylim([0,20])

species_axis.spines['top'].set_visible(False)
species_axis.spines['right'].set_visible(False)
species_axis.get_xaxis().tick_bottom()
species_axis.get_yaxis().tick_left()

pylab.figure(2,figsize=(7,7))

fig2 = pylab.gcf()

supplemental_grid  = gridspec.GridSpec(5,4, hspace=0.8, wspace=0.1)

density_axes = []
for sample_idx in xrange(0,len(all_samples)):   

	sample_name = all_samples[sample_idx]

	row_idx = (sample_idx / 4)
	col_idx = (sample_idx % 4)
	
	print row_idx, col_idx

	current_density_axis = plt.Subplot(fig2, supplemental_grid[row_idx,col_idx])
	fig2.add_subplot(current_density_axis)

	if col_idx==0:
		current_density_axis.set_ylabel('Total reads\n contributed')
	
	current_density_axis.semilogx([1e-01],[1e-01])
	current_density_axis.set_xlim([1e01,9e03])
    
	if (row_idx == 4):
		current_density_axis.set_xlabel('Read pairs per barcode')
	else:
		current_density_axis.set_xticklabels([])
		
	current_density_axis.set_yticks([])
	
	total_reads = total_barcode_map[sample_name][1]*300*1e-09
	day = sample_time_map[sample_name]+1    
    
	current_density_axis.set_title('Day %d (%d Gbp)' % (day, total_reads),fontsize=7)
	
	density_axes.append(current_density_axis)
	
fig2.savefig(parse_midas_data.analysis_directory+'supplemental_barcode_statistics.pdf',bbox_inches='tight')

for sample_idx in xrange(0,len(all_samples)):

    sample_name = all_samples[sample_idx]
    
    
    sys.stderr.write("Processing sample %s...\n" % sample_name)
    sys.stderr.write("Loading depth map...\n")
    barcode_depth_map = barcode_utils.parse_barcode_depth_map(sample_name)
    sys.stderr.write("Done!\n")
        
    
    # create barcode->species map
    sys.stderr.write("Collating species barcodes...\n")
    # first create intermediate data structure:
    # barcode_id->longgene->count
    barcode_longgene_weight_map = {} 
    for species_name in pangenome_species:
        
        # Don't use new species yet!
        if species_name=='new_species':
            continue
        
        # Make sure barcodes exist for this timepoint.
        # BG: aside from bugs, shouldn't this always be true? 
        if not barcode_utils.barcodes_exist(species_name, sample_name):
            continue
         
        # Load barcodes      
        allele_barcode_map, allele_error_map = barcode_utils.parse_allele_barcode_tuples(species_name, sample_name)

        for allele in allele_barcode_map.keys():
        
            if allele.endswith('|A') or allele.endswith('|R'):
                # a SNP allele, don't include
                continue
        
            if len(allele_barcode_map)==0:
                continue
        
            longgene = species_name
            if longgene not in longgene_id_map:
                longgene_id_map[longgene] = len(id_longgene_map)
                id_longgene_map.append(longgene)
            
            longgene_id = longgene_id_map[longgene]
        
            for barcode_id, barcode_weight in allele_barcode_map[allele]:
                if barcode_id not in barcode_longgene_weight_map:
                    barcode_longgene_weight_map[barcode_id] = {}
                if longgene_id not in barcode_longgene_weight_map[barcode_id]: 
                    barcode_longgene_weight_map[barcode_id][longgene_id] = 0.0
                barcode_longgene_weight_map[barcode_id][longgene_id] += barcode_weight

    sys.stderr.write("Done!\n")
    
    if len(barcode_longgene_weight_map)==0:
        continue
    
    # Plot distributions of things
    total_readss = [[] for d in readcloud_centers]
    total_mapped_readss = [[] for d in readcloud_centers]
    
    all_total_reads = 0
    all_total_mapped_reads = 0
    
    num_speciess = [[] for d in readcloud_centers]
    effective_num_speciess = [[] for d in readcloud_centers]
    
    for barcode_id in barcode_longgene_weight_map:
        
        weights = numpy.array(barcode_longgene_weight_map[barcode_id].values())
        
        total_reads = barcode_depth_map[barcode_id]
        total_mapped_reads = weights.sum()/2
        
        fs = weights*1.0/weights.sum()
        
        num_species = len(fs)*1.0
        effective_num_species = 1.0/(fs*fs).mean()**0.5
        
        if (total_reads>=10) and (total_reads<1e04):
            
            bin_idx = numpy.digitize([total_reads],bins=readcloud_bins)[0]-1
            
            # get bin idx
            #bin_idx = numpy.nonzero((total_reads>=readcloud_bins)*(total_reads<readcloud_bins))[0]
            
            total_readss[bin_idx].append(total_reads)
            total_mapped_readss[bin_idx].append(total_mapped_reads)
            num_speciess[bin_idx].append(num_species)            
            effective_num_speciess[bin_idx].append(effective_num_species)
    
    avg_num_species = numpy.zeros_like(readcloud_centers)
    avg_effective_num_species = numpy.zeros_like(avg_num_species)
    total_reads_in_bin = numpy.zeros_like(avg_num_species)
    
    for bin_idx in xrange(0,len(total_readss)):
            
        total_readss[bin_idx] = numpy.array(total_readss[bin_idx])
        total_mapped_readss[bin_idx] = numpy.array(total_mapped_readss[bin_idx])
        num_speciess[bin_idx] = numpy.array(num_speciess[bin_idx])
        effective_num_speciess[bin_idx] = numpy.array(effective_num_speciess[bin_idx])
    
        total_reads_in_bin[bin_idx] = total_readss[bin_idx].sum()
        avg_num_species[bin_idx] = numpy.median(num_speciess[bin_idx])
        avg_effective_num_species[bin_idx] = numpy.median(effective_num_speciess[bin_idx])
        
    total_reads = total_barcode_map[sample_name][1]
    
    line, = density_axes[sample_idx].semilogx(readcloud_centers, total_reads_in_bin*1.0/total_reads_in_bin.sum(),'r-')
    
    # Time to plot!
    if sample_name in desired_samples:
    	desired_sample_idx = desired_samples.index(sample_name)
    	color = colors[desired_sample_idx]
        
    	print sample_name, total_reads
    	line, = density_axis.semilogx(readcloud_centers, total_reads_in_bin*1.0/total_reads_in_bin.sum(),'-',color=color)
    #species_axis.semilogx(readcloud_centers, avg_num_species,'r.-')
    	species_axis.semilogx(readcloud_centers, avg_effective_num_species,'-',color=color)
    
    	legend_axis.plot([-1],[0],'-',color=color,label=('Day %d (%dGbp)' % (sample_time_map[sample_name]+1, long(total_reads*300*1e-09))))
    
density_axis.set_yticks([])
density_axis.set_xticklabels([])

species_axis.set_xlim([1e01,1e04])
species_axis.semilogx([1e0-01,1e-01])

legend_axis.legend(loc='lower center',frameon=False,fontsize=7,numpoints=1,ncol=4,handlelength=1)   

for idx in xrange(0,len(density_axes)):
	row_idx = (sample_idx / 4)
	col_idx = (sample_idx % 4)
	
	if (row_idx == 4):
		pass
	else:
		current_density_axis.set_xticklabels([])
		
    
fig.savefig(parse_midas_data.analysis_directory+'barcode_statistics.pdf',bbox_inches='tight')
fig2.savefig(parse_midas_data.analysis_directory+'supplemental_barcode_statistics.pdf',bbox_inches='tight')

