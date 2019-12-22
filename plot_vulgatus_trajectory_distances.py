###############################
#
# Rest of script begins here
#
################################
import matplotlib  
matplotlib.use('Agg') 
import pylab
import numpy
import sys
from math import log10
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from numpy.random import binomial, random_sample
import bz2
import parse_midas_data
import parse_timecourse_data
import matplotlib
import matplotlib.pyplot as plt
import timecourse_utils
import parse_patric
import calculate_preexisting_snps
import cluster_utils
import core_gene_utils
################################################################################
#
# Standard header to read in argument information
#
################################################################################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
parser.add_argument("-fstar", "--freq-threshold", type=float, help="Frequency has to exceed to be included",default=0.2)
parser.add_argument("--fraction-covered", type=float, help="Fraction of timepoints with sufficient coverage",default=0.5)

args = parser.parse_args()

debug = args.debug
chunk_size = args.chunk_size
fstar = args.freq_threshold
fraction_covered = args.fraction_covered
################################################################################

      


mpl.rcParams['font.size'] = 7.0
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'
    
pylab.figure(figsize=(7,2))
fig = pylab.gcf()
outer_grid = gridspec.GridSpec(1,2,width_ratios=[1,1],wspace=0.1)

# Inconsistency axis
barcode_axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(barcode_axis)

barcode_axis.set_xlabel('$d_{ij}$ (unique read clouds)')
barcode_axis.set_ylabel('Fraction SNV pairs')

barcode_axis.spines['top'].set_visible(False)
barcode_axis.spines['right'].set_visible(False)
barcode_axis.get_xaxis().tick_bottom()
barcode_axis.get_yaxis().tick_left()
barcode_axis.set_xlim([1e-01,1e03])
 

read_axis = plt.Subplot(fig, outer_grid[1])
fig.add_subplot(read_axis)
read_axis.set_xlabel('$d_{ij}$ (raw reads)')

read_axis.spines['top'].set_visible(False)
read_axis.spines['right'].set_visible(False)
read_axis.get_xaxis().tick_bottom()
read_axis.get_yaxis().tick_left()
read_axis.set_xlim([1e-01,1e03])

distance_axes = [barcode_axis, read_axis]

# Real version:        
snps_directories = [parse_midas_data.barcode_directory+"/barcode_snps/", parse_midas_data.default_snps_directory]
# Debug version:
#snps_directories = [parse_midas_data.barcode_directory+"/barcode_snps/", parse_midas_data.barcode_directory+"/barcode_snps/"]
                         
    
# Mininum coverage for frequency estimation vs interpolation 
min_coverage = 10
desired_samples = parse_timecourse_data.morteza_samples
snp_samples = desired_samples
species_name = "Bacteroides_vulgatus_57955"

# Load snv pair linkage classifications
file = open("%s_snv_snv_linkage_categories.txt" % species_name,"r")
snv_pair_linkage_categories = []
target_snvs = set()
for line in file:
    items = line.split()
    contig = items[0]
    focal_location = long(items[1])
    target_location = long(items[2])
    category = long(items[3])
    snv_pair_linkage_categories.append((contig,focal_location,target_location,category))
    
    target_snvs.add((contig, focal_location))
    target_snvs.add((contig, target_location))
    
     
# Output filename
filename = parse_midas_data.analysis_directory+('%s_trajectory_distances.pdf' % species_name)

sys.stderr.write("Loading core genes...\n")
core_genes = core_gene_utils.parse_core_genes(species_name)
personal_core_genes = core_gene_utils.parse_personal_core_genes(species_name)
non_shared_genes = personal_core_genes
shared_pangenome_genes = core_gene_utils.parse_shared_genes(species_name)
sys.stderr.write("Done! %d core genes and %d shared genes and %d non-shared genes\n" % (len(core_genes), len(shared_pangenome_genes), len(non_shared_genes)))


species_coverage_matrix, samples, species = parse_midas_data.parse_global_marker_gene_coverages()
sample_time_map = parse_timecourse_data.parse_sample_time_map()
ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, samples)
samples = numpy.array(samples)[sample_idxs]
species_coverage_matrix = species_coverage_matrix[:,sample_idxs]
total_coverage = species_coverage_matrix.sum(axis=0)
species_freq_matrix = numpy.clip(species_coverage_matrix*1.0/total_coverage,0, 2)    

for panel_idx in xrange(0,len(snps_directories)):
    snps_directory = snps_directories[panel_idx]
    axis = distance_axes[panel_idx]
    
    sample_size = len(snp_samples)
    sys.stderr.write("Proceeding with %d temporal samples!\n" % sample_size)

    cluster_As = [] # (construct a # sites x # samples x # bases (4) array)   
    cluster_Ds = []
    location_idx_map = {}
    final_line_number = 0
    while final_line_number >= 0:
    
        sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
        dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_variant_types=set(['1D','2D','3D','4D']),chunk_size=chunk_size,allowed_samples=snp_samples, initial_line_number=final_line_number, allowed_genes=non_shared_genes, snps_directory=snps_directory)
        sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
        snp_samples = dummy_samples
        
        for gene_name in allele_counts_map.keys():
    
            #if gene_name!='435590.9.peg.242':
            #    continue
    
            for var_type in allele_counts_map[gene_name].keys():
                locations = allele_counts_map[gene_name][var_type]['locations']
                allele_counts = allele_counts_map[gene_name][var_type]['alleles']
                if len(allele_counts)==0:
                    continue
                
                alts = allele_counts[:,:,0]
                depths = allele_counts.sum(axis=2)
                freqs = alts*1.0/(depths+(depths==0))
                
                for snp_idx in xrange(0,len(locations)):
                
                    contig,location = locations[snp_idx]
                    
                    if (contig,location) not in target_snvs:
                        continue
                
                    insufficient_coverage = ((depths[snp_idx,:]>0).sum() < fraction_covered*depths.shape[1])
                    low_frequency = ((freqs[snp_idx]<fstar).all() or (freqs[snp_idx]>(1-fstar)).all())
                
                    if insufficient_coverage or low_frequency:
                        continue
                    
                    location_idx_map[(contig,location)] = len(cluster_As)
                    cluster_As.append(alts[snp_idx,:])
                    cluster_Ds.append(depths[snp_idx,:])
                    
    cluster_As = numpy.array(cluster_As)
    cluster_Ds = numpy.array(cluster_Ds)
    
    sys.stderr.write("Calculating distance matrix between %d SNVs...\n" % cluster_As.shape[0])              
    distance_matrix, distance_matrix_1, distance_matrix_2 = cluster_utils.calculate_distance_matrix(cluster_As, cluster_Ds)   
    sys.stderr.write("Done!\n")              
    
    print numpy.median(distance_matrix)
    print distance_matrix.max()
    print distance_matrix.min()
    
    linkage_distance_distributions = {1:[], 2:[],3:[],4:[]}
    linkage_labels = {1: 'Perfect LD (strong)',2:'Perfect LD (weak)',3:'Complete LD',4:'Four haplotypes'}
    for contig, focal_location, target_location,category in snv_pair_linkage_categories:
    
        if (contig,focal_location) not in location_idx_map:
            continue
        
        if (contig,target_location) not in location_idx_map:
            continue
            
        i = location_idx_map[(contig,focal_location)]
        j = location_idx_map[(contig,target_location)]
        dij = distance_matrix[i,j]
        
        linkage_distance_distributions[category].append(dij)
        
        
    bins = numpy.logspace(-1,3,100)
    ds = numpy.array(bins,copy=True)[1:]
    bins[0] = 0
    bins[-1] = 1e08
    
    for category in sorted(linkage_distance_distributions,reverse=True):
        
        hist,bin_edges = numpy.histogram(linkage_distance_distributions[category],bins=bins)
        
        axis.semilogx(ds,hist*1.0/hist.sum(),label=linkage_labels[category])
        
        #axis.hist(linkage_distance_distributions[category], bins=bins,histtype='step',label=linkage_labels[category])
    axis.semilogx([1e-02],[1e-02],'k.')     
    if panel_idx==0:
        axis.legend(frameon=False,loc='upper left')


read_axis.set_yticks([])
barcode_axis.set_yticks([])

        
sys.stderr.write("Saving final image...\t")
fig.savefig(filename, bbox_inches='tight', transparent=True)
pylab.close(fig)
sys.stderr.write("Done!\n")

