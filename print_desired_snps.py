###############################
#
# Rest of script begins here
#
################################
import numpy
import sys
from math import log10
import bz2
import parse_midas_data
import parse_timecourse_data
import timecourse_utils

min_coverage = 5

################################################################################
#
# Standard header to read in argument information
#
################################################################################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("output_filename", help="output file")
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
args = parser.parse_args()

output_filename = args.output_filename
debug = args.debug
chunk_size = args.chunk_size
################################################################################

sample_time_map = parse_timecourse_data.parse_sample_time_map()
species_names = parse_midas_data.parse_good_species_list()
desired_samples = parse_timecourse_data.morteza_samples

if debug:
    species_names = [species_names[0]]

output_items = {}
for species_idx in xrange(0,len(species_names)):
    species_name = species_names[species_idx]
    
    sys.stderr.write("Processing %s...\n" % species_name)
       
    times = []
    alt_matrix = []
    depth_matrix = []
    snp_infos = []

    final_line_number = 0
    while final_line_number >= 0:
    
        sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
        samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_variant_types=set(['1D','2D','3D','4D']),chunk_size=chunk_size,allowed_samples=desired_samples, initial_line_number=final_line_number)
        sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))

        if len(samples)<5:
            continue
            
        sample_ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, samples)

        if len(times)==0:
            times = sample_ts

        # Calculate fixation matrix
        sys.stderr.write("Calculating allele freqs...\n")
        chunk_alts, chunk_depths, chunk_snp_infos = timecourse_utils.calculate_read_count_matrix(allele_counts_map, passed_sites_map, allowed_variant_types=set(['1D','2D','3D','4D']))    
        sys.stderr.write("Done!\n")

        chunk_alts = chunk_alts[:,sample_idxs]
        chunk_depths = chunk_depths[:,sample_idxs]
        
        # should be polarized already
        
        # now see which sites to track        
        desired_sites = ((chunk_alts>(0.1*chunk_depths)).sum(axis=1)>2.5)*((chunk_depths>0).sum(axis=1)>10)
    
        chunk_alts = chunk_alts[desired_sites,:]
        chunk_depths = chunk_depths[desired_sites,:]
        chunk_allele_freqs = chunk_alts*1.0/(chunk_depths+(chunk_depths==0))
                      
        if desired_sites.sum()>0:
            alt_matrix.append(chunk_alts)
            depth_matrix.append(chunk_depths)
            desired_site_idxs = numpy.nonzero(desired_sites)[0]
            for idx in desired_site_idxs:
                snp_infos.append(chunk_snp_infos[idx])
                    
    sys.stderr.write("Done loading SNVs!\n")
           
    if len(alt_matrix)>0:     
        alt_matrix = numpy.vstack(alt_matrix)
        depth_matrix = numpy.vstack(depth_matrix) 
    else:
        alt_matrix = numpy.array([])
        depth_matrix = numpy.array([])
    
    for mutation_idx in xrange(0,len(snp_infos)):
        
        chromosome, location, gene_name, variant_type, dummy_polarization = snp_infos[mutation_idx]
        
        alts = alt_matrix[mutation_idx,:]
        depths = depth_matrix[mutation_idx,:]
        
        if (depths>=min_coverage).sum() < 2:
            continue
        
        freqs = alts*1.0/(depths+(depths==0))
        masked_times = times[depths>=min_coverage]
        masked_freqs = freqs[depths>=min_coverage]
        masked_depths = depths[depths>=min_coverage]
        
        if (masked_freqs>0.1).sum() < 2.5:
            continue
        
        if True: 
            output_str = "\t".join([species_name, snp_infos[mutation_idx][0], str(snp_infos[mutation_idx][1])])
        
            output_key = (species_idx, chromosome, location)
            output_items[output_key] = output_str

output_file = open(output_filename,"w")        
output_file.write( "\t".join(["species_id","contig","pos"]) )
output_file.write("\n")
for output_key in sorted(output_items.keys()):
    output_file.write(output_items[output_key])
    output_file.write("\n")                     
output_file.close()       
sys.stderr.write("Done!\n")

