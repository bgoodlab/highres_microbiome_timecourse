###############################
#
# Import tools
#
################################
import matplotlib  
matplotlib.use('Agg') 
import pylab
import numpy
import sys
from math import log10, floor, ceil,log
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from numpy.random import binomial, random_sample, shuffle,random
import bz2
import parse_midas_data
import parse_timecourse_data
import matplotlib
import matplotlib.pyplot as plt
import timecourse_utils
import parse_patric
#import calculate_preexisting_snps
import calculate_snp_prevalences
import cluster_utils
import core_gene_utils
#import calculate_within_species_fixations
import calculate_barcode_within_species_fixations as calculate_within_species_fixations
import figure_utils
import config
import diversity_utils
import sfs_utils
import stats_utils


desired_samples = parse_timecourse_data.morteza_samples
snps_directory=parse_midas_data.barcode_directory+"/barcode_snps/"


pylab.figure(1,figsize=(5, 3))
fig = pylab.gcf()
coverage_axis = pylab.gca()

coverage_axis.set_ylabel('Genome-wide coverage, $\overline{\mathcal{D}}_t$')
coverage_axis.set_xlabel('Species relative abundance')

# First 


desired_samples = parse_timecourse_data.morteza_samples    

species_coverage_matrix, samples, species = parse_midas_data.parse_global_marker_gene_coverages()
sample_time_map = parse_timecourse_data.parse_sample_time_map()
species_ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, samples,min_time=-1)
samples = numpy.array(samples)[sample_idxs]
species_coverage_matrix = species_coverage_matrix[:,sample_idxs]
total_coverage = species_coverage_matrix.sum(axis=0)
species_freq_matrix = numpy.clip(species_coverage_matrix*1.0/total_coverage,0, 2)    


for species_idx in xrange(0,len(species)):        

    species_name = species[species_idx]
    
    sample_coverage_map = parse_midas_data.parse_median_coverage_map(species_name,snps_directory=snps_directory)
    
    
    if len(sample_coverage_map) == 0:
        continue
        
    for sample_idx in xrange(0,len(desired_samples)):
        
        sample_name = desired_samples[sample_idx]
        abundance = species_freq_matrix[species_idx,sample_idx]
        
        if sample_name not in sample_coverage_map:
            continue
            
        Dbar = sample_coverage_map[sample_name]
        
        if sample_name in parse_timecourse_data.highcoverage_samples:
            color='b'
            zorder = 2
        else:
            color='0.7'
            zorder = 1
        
        
        coverage_axis.loglog([abundance],[Dbar],'.',color=color,zorder=zorder)


coverage_axis.loglog([3e-04,1],[10,10],'k:')
        
coverage_axis.set_ylim([1,1e03])
coverage_axis.set_xlim([3e-04,1])



coverage_axis.plot([1],[0.1],'.',color='b',label='Full lane samples')
coverage_axis.plot([1],[0.1],'.',color='0.7',label='1/3 lane samples')
coverage_axis.legend(loc='upper left',frameon=False,numpoints=1)
fig.savefig(parse_midas_data.analysis_directory+'supplemental_species_coverage.pdf',bbox_inches='tight')