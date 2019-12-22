import config
import parse_timecourse_data
import parse_midas_data
import gzip
import pylab
import numpy
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import sys

mpl.rcParams['font.size'] = 7
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

samples = parse_timecourse_data.morteza_samples

species_coverage_matrix, samples, species = parse_midas_data.parse_global_marker_gene_coverages()
sample_time_map = parse_timecourse_data.parse_sample_time_map()
ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, samples)
samples = numpy.array(samples)[sample_idxs]

species_freq_matrix = species_coverage_matrix*1.0/species_coverage_matrix.sum(axis=0)

desired_species = sys.argv[1]

for species_idx in xrange(0,len(species)):
    species_name = species[species_idx]
    if species_name.startswith(desired_species):
        pylab.plot(ts, species_freq_matrix[species_idx,sample_idxs],'k.-')
        
pylab.show()