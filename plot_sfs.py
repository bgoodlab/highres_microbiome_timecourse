import matplotlib  
matplotlib.use('Agg') 
import parse_midas_data
import sample_utils
import pylab
import sys
import numpy
from numpy.random import normal
#from calculate_pi_matrix import calculate_self_pis
import diversity_utils
import stats_utils
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint,binomial


import config
import sfs_utils
import stats_utils
import figure_utils

fontsize = 6
mpl.rcParams['font.size'] = fontsize
mpl.rcParams['lines.linewidth'] = 1.0
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'


min_coverage = config.min_median_coverage

species_name = "Bacteroides_vulgatus_57955"
#species_name = "Bacteroides_uniformis_57318"
species_name = "Eubacterium_eligens_61678"
species_name = "Blautia_wexlerae_56130"
sample_1 = '1021' # complicated polyploid
sample_2 = '1022.1'
sample_3 = '1014.2'
sample_4 = '4021A'


haploid_color = '#08519c'
light_haploid_color = '#6699CC'
diploid_color = '#de2d26' #'#fb6a4a' #
transition_color = '#756bb1'


# Load SNP information for species_name
sys.stderr.write("Loading SFSs for %s...\t" % species_name)
samples, sfs_map = parse_midas_data.parse_within_sample_sfs(species_name, allowed_variant_types=set(['4D'])) 
sys.stderr.write("Done!\n")

# Load genomic coverage distributions
sample_coverage_histograms, samples = parse_midas_data.parse_coverage_distribution(species_name)
median_coverages = numpy.array([stats_utils.calculate_nonzero_median_from_histogram(sample_coverage_histogram) for sample_coverage_histogram in sample_coverage_histograms])
sample_coverage_map = {samples[i]: median_coverages[i] for i in xrange(0,len(samples))}
samples = numpy.array(samples)

median_coverages = numpy.array([sample_coverage_map[samples[i]] for i in xrange(0,len(samples))])

# Only plot samples above a certain depth threshold
desired_samples = samples[(median_coverages>=min_coverage)]
desired_median_coverages = numpy.array([sample_coverage_map[sample] for sample in desired_samples])

pylab.figure(1,figsize=(1.2,3))
fig = pylab.gcf()
# make three panels
outer_grid  = gridspec.GridSpec(1,1)

sfs_grid = gridspec.GridSpecFromSubplotSpec(4, 1, height_ratios=[1,1,1,1],
                subplot_spec=outer_grid[0], hspace=0.25)
                
sfs_axis_1 = plt.Subplot(fig, sfs_grid[0])
fig.add_subplot(sfs_axis_1)

sfs_axis_1.set_title('Sample 1 ($\\overline{D}=%d$)' % sample_coverage_map[sample_1],fontsize=5,y=0.9)
sfs_axis_1.set_xticks([10*i for i in xrange(0,11)])
sfs_axis_1.set_xticklabels([])
sfs_axis_1.set_xlim([50,100])
sfs_axis_1.set_yticks([])
sfs_axis_1.xaxis.tick_bottom()


sfs_axis_2 = plt.Subplot(fig, sfs_grid[1])
fig.add_subplot(sfs_axis_2)

sfs_axis_2.set_title('Sample 2 ($\\overline{D}=%d$)' % sample_coverage_map[sample_2],fontsize=5,y=0.9)
sfs_axis_2.set_xticks([10*i for i in xrange(0,11)])
sfs_axis_2.set_xticklabels([])
sfs_axis_2.set_xlim([50,100])
sfs_axis_2.set_yticks([])
sfs_axis_2.xaxis.tick_bottom()

#sfs_axis_2.set_ylabel('Fraction of 4D sites')

sfs_axis_3 = plt.Subplot(fig, sfs_grid[2])
fig.add_subplot(sfs_axis_3)

sfs_axis_3.set_title('Sample 3 ($\\overline{D}=%d$)' % sample_coverage_map[sample_3],fontsize=5,y=0.9)
sfs_axis_3.set_xticks([10*i for i in xrange(0,11)])
sfs_axis_3.set_xticklabels([])
sfs_axis_3.set_xlim([50,100])
sfs_axis_3.set_yticks([])
sfs_axis_3.set_ylabel('                  Fraction of synonymous sites')
sfs_axis_3.xaxis.tick_bottom()


sfs_axis_4 = plt.Subplot(fig, sfs_grid[3])
fig.add_subplot(sfs_axis_4)

sfs_axis_4.set_title('Sample 4 ($\\overline{D}=%d$)' % sample_coverage_map[sample_4],fontsize=5,y=0.9)
sfs_axis_4.set_xlabel('Major allele freq (%)')

sfs_axis_4.set_xticks([10*i for i in xrange(0,11)])
sfs_axis_4.set_xlim([50,100])
sfs_axis_4.set_yticks([])
sfs_axis_4.xaxis.tick_bottom()


###################################
#
# Plot example SFSs
#
###################################

# Sample 1
fs,pfs = sfs_utils.calculate_binned_sfs_from_sfs_map(sfs_map[sample_1],folding='major')
df = fs[1]-fs[0]

within_sites, between_sites, total_sites = sfs_utils.calculate_polymorphism_rates_from_sfs_map(sfs_map[sample_1])

    

between_line = between_sites*1.0/total_sites/((fs>0.2)*(fs<0.5)).sum()
pmax = pfs[(fs>0.5)*(fs<0.90)].max()

within_rate = within_sites*1.0/total_sites
print "Sample 1: within =", within_rate, "avg-distance =", between_sites*1.0/total_sites

sfs_axis_1.fill_between([80,100],[0,0],[1,1],color='0.8')

#sfs_axis_1.fill_between([20,100],[0,0],[1,1],color='0.8')
sfs_axis_1.bar((fs-df/2)*100,pfs,width=df, edgecolor=light_haploid_color, color=light_haploid_color)
line, = sfs_axis_1.plot([20,80], [between_line,between_line], 'k-',linewidth=0.35)
line.set_dashes((1.5,1))
sfs_axis_1.set_ylim([0,pmax*2])

# Sample 2
fs,pfs = sfs_utils.calculate_binned_sfs_from_sfs_map(sfs_map[sample_2],folding='major')
df = fs[1]-fs[0]
within_sites, between_sites, total_sites = sfs_utils.calculate_polymorphism_rates_from_sfs_map(sfs_map[sample_2])
between_line = between_sites*1.0/total_sites/((fs>0.2)*(fs<0.5)).sum()

within_rate = within_sites*1.0/total_sites
print "Sample 2: within =", within_rate, "avg-distance =", between_sites*1.0/total_sites

pmax = pfs[(fs>0.5)*(fs<0.90)].max() 
#pmax = between_line
sfs_axis_2.fill_between([80,100],[0,0],[1,1],color='0.8')

sfs_axis_2.bar((fs-df/2)*100,pfs,width=df,edgecolor=light_haploid_color, color=light_haploid_color)
line, = sfs_axis_2.plot([20,80], [between_line,between_line], 'k-',linewidth=0.35)
line.set_dashes((1.5,1))

sfs_axis_2.set_ylim([0,pmax*2])

# Sample 3
fs,pfs = sfs_utils.calculate_binned_sfs_from_sfs_map(sfs_map[sample_3],folding='major')
df = fs[1]-fs[0]
within_sites, between_sites, total_sites = sfs_utils.calculate_polymorphism_rates_from_sfs_map(sfs_map[sample_3])
between_line = between_sites*1.0/total_sites/((fs>0.2)*(fs<0.5)).sum()
pmax = pfs[(fs>0.5)*(fs<0.90)].max() 
#pmax = between_line

within_rate = within_sites*1.0/total_sites
print "Sample 3: within =", within_rate, "avg-distance =", between_sites*1.0/total_sites
sfs_axis_3.fill_between([80,100],[0,0],[1,1],color='0.8')

sfs_axis_3.bar((fs-df/2)*100,pfs,width=df,edgecolor=haploid_color,color=haploid_color)
line, = sfs_axis_3.plot([20,80], [between_line,between_line], 'k-',linewidth=0.35)
line.set_dashes((1.5,1))

sfs_axis_3.set_ylim([0,pmax*2])

# Sample 3
fs,pfs = sfs_utils.calculate_binned_sfs_from_sfs_map(sfs_map[sample_4],folding='major')
df = fs[1]-fs[0]
within_sites, between_sites, total_sites = sfs_utils.calculate_polymorphism_rates_from_sfs_map(sfs_map[sample_4])
between_line = between_sites*1.0/total_sites/((fs>0.2)*(fs<0.5)).sum()
pmax = pfs[(fs>0.5)*(fs<0.90)].max() 
#pmax = between_line

within_rate = within_sites*1.0/total_sites
print "Sample 4: within =", within_rate, "avg-distance =", between_sites*1.0/total_sites
sfs_axis_4.fill_between([80,100],[0,0],[1,1],color='0.8')

sfs_axis_4.bar((fs-df/2)*100,pfs,width=df,edgecolor=haploid_color,color=haploid_color)
line, = sfs_axis_4.plot([20,80], [between_line,between_line], 'k-',linewidth=0.35)
line.set_dashes((1.5,1))

sfs_axis_4.set_ylim([0,pmax*2])
filename = parse_midas_data.analysis_directory+('%s_sfs.pdf' % species_name) 

fig.savefig(filename, bbox_inches='tight')