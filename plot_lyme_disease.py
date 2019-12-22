import matplotlib  
matplotlib.use('Agg') 
import parse_midas_data
import parse_timecourse_data
import bacterial_phylogeny_utils
import pylab
import sys
import numpy
from math import log10, fabs, log
import figure_utils

import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint

mpl.rcParams['font.size'] = 7
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

species_coverage_matrix, samples, species_list = parse_midas_data.parse_global_marker_gene_coverages()
sample_time_map = parse_timecourse_data.parse_sample_time_map()
ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, samples, min_time=-1000)

max_dt = 2.0


dts = numpy.diff(ts)*1.0

left_dts = numpy.hstack([[500.0],dts/4])
right_dts = numpy.hstack([dts/4,[500.0]])

dts = numpy.fmin(left_dts, right_dts)
dts = numpy.fmin(dts,max_dt/2)

print "ts:", ts

print "dts: ", dts

#ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, parse_timecourse_data.highcoverage_samples)

species_coverage_matrix = species_coverage_matrix[:,sample_idxs]

total_coverage = species_coverage_matrix.sum(axis=0)
species_freq_matrix = numpy.clip(species_coverage_matrix*1.0/total_coverage,0, 2)    

pylab.figure(figsize=(5,2))
pylab.xlabel('Time (days)')
pylab.ylabel('Relative abundance')
for species_idx in xrange(0,len(species_list)):
	
	species_name = species_list[species_idx]
	print species_name
	if species_name.startswith('Borrelia'):
		print species_name
		pylab.plot(ts,species_freq_matrix[species_idx,:],'k.-')
		
pylab.savefig('lyme_abundance.pdf',bbox_inches='tight')