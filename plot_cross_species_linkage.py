import sys
import numpy
import pylab
import parse_midas_data
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint
import figure_utils
import parse_timecourse_data

good_species_list = parse_midas_data.parse_good_species_list()
good_species_pretty_list = figure_utils.get_pretty_species_names(good_species_list)
pretty_species_name_map = {species_name: pretty_species_name for species_name, pretty_species_name in zip(good_species_list, good_species_pretty_list)}


mpl.rcParams['font.size'] = 7
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

filename = 'fixation_barcode_output_core.txt'
file = open(filename,"r")

fixation_linkage_data = {}

line=file.readline().strip()
while line!="":    
    items = line.split("|")
    species_name = items[0]
    contig = items[1]
    gene_name = items[2]
    location = long(items[3])
    
    if species_name not in fixation_linkage_data:
        fixation_linkage_data[species_name] = []
    
    snp = (contig, gene_name, location)
    
    snp_linkage_data = {}
    for i in xrange(0,2):    
        line = file.readline()
        items = line.split(",")
    
        allele = items[0].strip()
        total_barcodes = long(items[1])
    
        allele_linkages = []
        if len(items)>3:
        
            for item in items[3:]:
                subitems = item.split("|")
                other_species = subitems[0].strip()
                other_gene = subitems[1].strip()
                observed = long(subitems[2])
                expected = float(subitems[3])
            
                allele_linkages.append((other_species, other_gene, observed, expected))
        
        snp_linkage_data[allele] = allele_linkages
    
    fixation_linkage_data[species_name].append((snp, snp_linkage_data))
    line = file.readline()

species_list = fixation_linkage_data.keys()
species_list.sort()

print len(species_list)

pylab.figure(1,figsize=(7,1.5))
fig = pylab.gcf()
outer_grid  = gridspec.GridSpec(2,1,height_ratios=[0.05, 0.7],hspace=0.2)

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
   

legend_axis.plot([-2],[-1],'s',color='0.7',markersize=3,label='SNVs tested ($\leq 1000$)',markeredgewidth=0)
legend_axis.plot([-2],[-1],'s',color='b',markersize=3,label='Confirmed linked to correct core genome',markeredgewidth=0)
legend_axis.plot([-2],[-1],'s',color='r',markersize=3,label='Different core genome',markeredgewidth=0)


legend_axis.legend(loc='lower right',frameon=False,fontsize=7,numpoints=1,ncol=5,handlelength=1)   

log_axes = []

log_axis = plt.Subplot(fig, outer_grid[1])
fig.add_subplot(log_axis)
log_axis.set_ylim([0.5,1e03])
log_axis.set_ylabel('SNV differences\nfrom baseline',fontsize=7)
log_axes.append(log_axis)

#log_axis = plt.Subplot(fig, outer_grid[2])
#fig.add_subplot(log_axis)
#log_axis.set_ylim([0.5,1e03])
#log_axis.set_ylabel('SNV differences\nfrom previous',fontsize=7)
#log_axes.append(log_axis)

all_total = 0
all_pos_confirmed = 0
all_neg_confirmed = 0

non_replacement_total = 0
non_replacement_pos_confirmed = 0
non_replacement_neg_confirmed = 0

import calculate_barcode_within_species_fixations as calculate_within_species_fixations

epochs = [parse_timecourse_data.initial_epoch_intervals] #,parse_timecourse_data.previous_epoch_intervals]

for epoch_idx in xrange(0,len(epochs)):

    
    for species_idx in xrange(0,len(species_list)):
        species_name = species_list[species_idx]
        print species_name
        num_total = 0
        num_pos_confirmed = 0
        num_neg_confirmed = 0
    
    
        snp_changes, gene_changes = calculate_within_species_fixations.load_within_species_fixations(species_name,allowed_epochs=epochs[epoch_idx])
    
        allowed_snp_changes = set()
        for contig, location, gene_name, variant_type, initial_freq, final_freq in snp_changes:
            allowed_snp_changes.add((contig, gene_name, location))
    
        
        for snp, snp_linkage_data in fixation_linkage_data[species_name]:
            if snp not in allowed_snp_changes:
                continue
        
            allele_pos_confirmations = []
            allele_neg_confirmations = []
        
            for allele in sorted(snp_linkage_data):
            
                pos_confirmed = False
                neg_confirmed = False
            
                if len(snp_linkage_data[allele])>0:     
                
                    num_pos_barcodes = 0
                    num_neg_barcodes = 0
            
                    for other_species, other_gene, observed, expected in snp_linkage_data[allele]:
                        if (other_species==species_name):
                            num_pos_barcodes += observed
                        else:
                            num_neg_barcodes += observed
            
                    if num_pos_barcodes>2*num_neg_barcodes:
                        pos_confirmed=True
                    else:
                        neg_confirmed=True
            
            
                allele_pos_confirmations.append( pos_confirmed)
                allele_neg_confirmations.append( neg_confirmed)
            
            num_total += 1
            
            
            if any(allele_neg_confirmations):
                num_neg_confirmed+=1
                if neg_confirmed and len(snp_linkage_data[allele])<100:
                    print snp, allele
                    print snp_linkage_data[allele]
            
            elif all(allele_pos_confirmations):
                num_pos_confirmed+=1
            else:
                pass
            
        all_total += num_total
        all_pos_confirmed += num_pos_confirmed
        all_neg_confirmed += num_neg_confirmed

        if num_total<100:
            non_replacement_total += num_total
            non_replacement_pos_confirmed += num_pos_confirmed
            non_replacement_neg_confirmed += num_neg_confirmed
        
    
        print species_name, num_total, num_pos_confirmed, num_neg_confirmed

        log_axis = log_axes[epoch_idx]

        if num_total>0:
            log_axis.bar([species_idx-0.25], [num_total],width=0.5,color='0.7',log=True,edgecolor='none')

        if num_pos_confirmed>0:
            log_axis.bar([species_idx-0.25], [num_pos_confirmed],width=0.25,color='b',log=True,edgecolor='none')
    
        if num_neg_confirmed>0:
            log_axis.bar([species_idx], [num_neg_confirmed],width=0.25,color='r',log=True,edgecolor='none')
    
    log_axis.set_xlim([-1,len(species_list)])  

    log_axis.set_xticks( numpy.arange(0,len(species_list)))

    xticklabels = []
    for species_name in species_list:
        xticklabels.append(pretty_species_name_map[species_name])
    

    if epoch_idx==(len(epochs)-1):
        log_axis.set_xticklabels( xticklabels,rotation=90)
    else:
        log_axis.set_xticklabels([])
    log_axis.tick_params(axis='x', labelsize=6,direction='out',length=3,pad=1)
    log_axis.get_xaxis().tick_bottom()
    log_axis.get_yaxis().tick_left()
    #log_axis.spines['top'].set_visible(False)
    #log_axis.spines['right'].set_visible(False)
    
    
print 'all', all_total, all_pos_confirmed, (all_pos_confirmed*1.0/all_total), all_neg_confirmed, all_neg_confirmed*1.0/(all_neg_confirmed+all_pos_confirmed)
print 'modification_only', non_replacement_total, non_replacement_pos_confirmed, (non_replacement_pos_confirmed*1.0/non_replacement_total), non_replacement_neg_confirmed, (non_replacement_neg_confirmed+1.0)*1.0/(non_replacement_neg_confirmed+non_replacement_pos_confirmed+1.0)


fig = pylab.gcf()
fig.savefig('%s/cross_species_sweep_linkage.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
          
        