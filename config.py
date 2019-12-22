###############################################################################
#
# Set up default source and output directories
#
###############################################################################
import os.path 
from math import log10

data_directory = os.path.expanduser("~/highres_microbiome_timecourse_data_022119/")
# Pre-022019 version
#data_directory = os.path.expanduser("~/highres_microbiome_timecourse_data/")
# V. old version
#data_directory = os.path.expanduser("~/highres_microbiome_timecourse_data_old/")

analysis_directory = os.path.expanduser("~/highres_microbiome_timecourse_analysis/")
scripts_directory = os.path.expanduser("~/Dropbox/highres_microbiome_timecourse/")
patric_directory = os.path.expanduser("~/patric_db/")
uniref_directory = os.path.expanduser("~/uniref_db/")
midas_directory = os.path.expanduser("~/midas_db/")

barcode_directory = os.path.expanduser("~/highres_microbiome_timecourse_barcode_data_022119/")
# Pre-022019 version
#barcode_directory = os.path.expanduser("~/highres_microbiome_timecourse_barcode_data/")
#barcode_directory = os.path.expanduser("~/highres_microbiome_timecourse_barcode_data_old/")

humann2_directory = os.path.expanduser("~/highres_microbiome_timecourse_humann2_data/")
# We use this one to debug because it was the first one we looked at
debug_species_name = 'Bacteroides_uniformis_57318'

good_species_min_coverage = 10
good_species_min_prevalence = 3

min_median_coverage = 20

# Genome-wide coverage must be at least this high 
# otherwise timepoint is masked for all sites
barcode_min_median_coverage = 5
# At least one timepoint must have "good" coverage at the site
# otherwise entire trajectory is masked for that site
barcode_min_good_coverage = 10
# A single timepoint must have at least this coverage
# otherwise it is masked for this site
barcode_min_fixation_coverage = 5
# Expected error rate must be at least this small,
# otherwise will not be checked for SNV difference
barcode_max_null_fixation_pvalue = 1e-03

consensus_lower_threshold = 0.2
consensus_upper_threshold = 0.8
fixation_min_change = consensus_upper_threshold-consensus_lower_threshold
fixation_log10_depth_ratio_threshold = log10(3)

threshold_within_between_fraction = 0.1
threshold_pi = 1e-03

modification_difference_threshold = 20
replacement_difference_threshold = 500

twin_modification_difference_threshold = 1000
twin_replacement_difference_threshold = 1000

# Older, permissive version
#gainloss_max_absent_copynum = 0.05
#gainloss_min_normal_copynum = 0.5
#gainloss_max_normal_copynum = 2

# Settings we used for final version of PLoS Bio paper
gainloss_max_absent_copynum = 0.05
gainloss_min_normal_copynum = 0.6
gainloss_max_normal_copynum = 1.2

core_genome_min_marker_coverage = 10 # coverage to use a point to assess core genome
core_genome_min_copynum = 0.5
core_genome_max_copynum = 2 # BG: should we use a maximum for "core genome"? I'm going to go w/ yes for now
core_genome_min_prevalence = 0.7
core_genome_min_samples = 3

shared_genome_min_copynum = 3

snv_sweep_min_present_fraction = 0.7
gene_sweep_max_any_copynum = gainloss_max_normal_copynum


# Default parameters for pipe snps
# (Initial filtering for snps, done during postprocessing)
pipe_snps_min_samples=4
pipe_snps_min_nonzero_median_coverage=5
pipe_snps_lower_depth_factor=0.3
pipe_snps_upper_depth_factor=3

parse_snps_min_freq = 0.1

between_host_min_sample_size = 33
within_host_min_sample_size = 3
within_host_min_haploid_sample_size = 10

cluster_distance_threshold_barcodes = 3.5
#cluster_distance_threshold_barcodes = 3.5 #3 # 4
cluster_distance_threshold_reads = 7 #8

from parse_timecourse_data import *
