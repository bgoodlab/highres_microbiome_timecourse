###########################################################################
#
# The purpose of this script is to figure out which genes (and therefore which species)
# are linked to a given allele
#
###########################################################################


import sys
import numpy
import parse_midas_data
import parse_timecourse_data
import stats_utils
import barcode_utils
import gzip
import collections
import midas_db_utils

import os.path
import config
import cPickle as pickle
from math import ceil, fabs
import sys
import argparse
from scipy.special import gammaln as loggamma
from scipy.special import gammainc
from scipy.special import expm1
from math import log,exp
from scipy.stats import poisson

from numpy.random import shuffle

pangenome_species = parse_midas_data.parse_pangenome_species()

genome_length_map = midas_db_utils.parse_genome_length_map()

    
###########################################################################
#
# Standard header to read in argument information
#
###########################################################################


gamete_idx_map = {('R','R'): 0, ('R','A'):1, ('A','R'):2, ('A','A'):3  }
allele_idx_map = {'R':0, 'A':1}

# long snp map
# long snp = (species, "contig|location") tuple (allele is stored separately)
# id = int

shortsnp_id_map = {}
id_shortsnp_map = []

species_id_map = {}
id_species_map = []

contig_id_map = {}
id_contig_map = []

min_S = 10

desired_species = sys.argv[1]
desired_speciess = [desired_species]
shortsnp_counter_map = {}


distance_bins = numpy.hstack([numpy.arange(0,1e05), numpy.logspace(5,7,20)])
distance_centers = numpy.array(distance_bins[0:-1],copy=True)
distance_bins-=0.5


desired_samples = parse_timecourse_data.morteza_samples
#desired_samples = ['6037'] 

G_function_map = barcode_utils.parse_G_function_map()
    
for sample_idx in xrange(0,len(desired_samples)):
    # Process each sample separately
    # (combine results at end)
    
    sample_name = desired_samples[sample_idx]
        
    G_function = G_function_map[sample_name]
       
    sys.stderr.write("Processing sample %s...\n" % sample_name)
    #sys.stderr.write("Loading depth map...\n")
    barcode_depth_map = barcode_utils.parse_barcode_depth_map(sample_name)
    #sys.stderr.write("Done!\n")
        
    # create barcode->species map
    sys.stderr.write("Collating species barcodes...\n")
        
    shortsnp_barcode_map = {} 
    barcode_shortsnp_map = {}

    for species_name in desired_speciess:
        
        # Make sure barcodes exist for this species at this timepoint.
        if not barcode_utils.barcodes_exist(species_name, sample_name):
            continue
         
        # Load barcodes      
        allele_barcode_map = barcode_utils.parse_allele_barcode_tuples(species_name, sample_name)

        for allele in allele_barcode_map:
        
            if not (allele[-1]=='A' or allele[-1]=='R'):
                # Not a snp, ignore
                continue
                    
            snp_str = allele[:-2]
            polarization = allele[-1]
            allele_idx = allele_idx_map[polarization]
            
            items = snp_str.split("|")
            contig = (species_name, items[0])
            location = long(items[1])
            
            if species_name not in species_id_map:
                species_id = len(id_species_map)
                id_species_map.append(species_name)
                species_id_map[species_name] = species_id    
            species_id = species_id_map[species_name]
            
            if contig not in contig_id_map:
                contig_id = len(id_contig_map)
                id_contig_map.append(contig)
                contig_id_map[contig] = contig_id
            contig_id = contig_id_map[contig]
            
            shortsnp = (species_id, contig_id, location)
                
            if shortsnp not in shortsnp_id_map:
                shortsnp_id = len(id_shortsnp_map)
                shortsnp_id_map[shortsnp] = shortsnp_id
                id_shortsnp_map.append(shortsnp)
                shortsnp_counter_map[shortsnp_id] =  collections.Counter() 
                 
            shortsnp_id = shortsnp_id_map[shortsnp]
                                   
            # Add in the barcodes
            for barcode_id, barcode_weight in allele_barcode_map[allele]:
                
                # don't include barcode if too little coverage
                # Poor man's error correction
                if barcode_id not in barcode_depth_map: 
                    continue
                
                D = barcode_depth_map[barcode_id]
                G = G_function(D)
                
                if G>10:
                    continue
                
                        
                if shortsnp_id not in shortsnp_barcode_map:
                    shortsnp_barcode_map[shortsnp_id] = []
                        
                # Add barcode to list
                shortsnp_barcode_map[shortsnp_id].append( barcode_id)
                    
                if barcode_id not in barcode_shortsnp_map:
                    barcode_shortsnp_map[barcode_id] = []
                
                barcode_shortsnp_map[barcode_id].append( shortsnp_id)
                
    sys.stderr.write("Done!\n")            
    
    sys.stderr.write("Counting shared barcodes for sites!\t")                
    for focal_shortsnp_id in shortsnp_barcode_map:
        for barcode_id in shortsnp_barcode_map[focal_shortsnp_id]: 
            shortsnp_counter_map[    focal_shortsnp_id].update(barcode_shortsnp_map[barcode_id])
        if focal_shortsnp_id not in shortsnp_counter_map[focal_shortsnp_id]:
            sys.stderr.write("Problem already!\n")
    sys.stderr.write("Done!\n")
    
    sys.stderr.write("Tallying pairs...\t")
    total_pairs = 0
    for focal_shortsnp_id in xrange(0,len(id_shortsnp_map)):
        total_pairs += len(shortsnp_counter_map[focal_shortsnp_id])
    sys.stderr.write("%0.1fe06\n" % (total_pairs*1e-06))
                
sys.stderr.write("Done!\n")

sys.stderr.write("Postprocessing...\n")
B_counts = numpy.zeros_like(distance_centers)
S_counts = numpy.zeros_like(distance_centers)

B_distribution = numpy.array([shortsnp_counter_map[focal_shortsnp_id][focal_shortsnp_id] for focal_shortsnp_id in shortsnp_counter_map])
B_distribution.sort()

lower_B = B_distribution[long(0.1*len(B_distribution))]-1
upper_B = B_distribution[long(0.9*len(B_distribution))]+1

all_contig_location_map = {}
for focal_shortsnp_id in shortsnp_counter_map:
        
    focal_species_id, focal_contig_id, focal_location = id_shortsnp_map[focal_shortsnp_id]

    if focal_contig_id not in all_contig_location_map:
        all_contig_location_map[focal_contig_id] = []
    all_contig_location_map[focal_contig_id].append(focal_location)
    
for contig_id in all_contig_location_map:
    all_contig_location_map[contig_id] = numpy.array(all_contig_location_map[contig_id])
sys.stderr.write("Done!\n")

sys.stderr.write("Binning counts by distance...\n")

for focal_shortsnp_id in shortsnp_counter_map:
    
    B = shortsnp_counter_map[focal_shortsnp_id][focal_shortsnp_id]
    
    if B<min_S:
        continue
    
    if B<=lower_B or B>=upper_B:
        continue
    
    focal_species_id, focal_contig_id, focal_location = id_shortsnp_map[focal_shortsnp_id]
    
    
    approx_genome_length = genome_length_map[ id_species_map[focal_species_id]]
    
    target_locations = []
    target_Ss = []
    
    for target_shortsnp_id in shortsnp_counter_map[focal_shortsnp_id]:
    
        target_species_id, target_contig_id, target_location = id_shortsnp_map[target_shortsnp_id]
        
        if target_contig_id!=focal_contig_id:
            continue
            
        S = shortsnp_counter_map[focal_shortsnp_id][target_shortsnp_id]
        
        target_locations.append(target_location)
        target_Ss.append(S)
        
    target_locations = numpy.array(target_locations)
    target_Ss = numpy.array(target_Ss)
    target_Bs = numpy.ones_like(target_Ss)*B
    distances = numpy.fabs(focal_location-target_locations)
    
    distances = numpy.fmin(distances, approx_genome_length-distances)
        
    distance_idxs = numpy.digitize(distances, bins=distance_bins)-1
    
    numpy.add.at(S_counts, distance_idxs, target_Ss)
    
    
    all_distances = numpy.fabs(focal_location-all_contig_location_map[focal_contig_id])
    all_distances = numpy.fmin(all_distances, approx_genome_length-all_distances)
        
    all_distance_idxs = numpy.digitize(all_distances, bins=distance_bins)-1
    
    numpy.add.at(B_counts, all_distance_idxs, B)
                
sys.stderr.write("Done!\n")

output_filename = "%s/new_barcode_trajectories/%s_distance_counts.txt" % (config.barcode_directory, desired_species)

file = open(output_filename,"w")
file.write(", ".join([str(d) for d in distance_centers]))
file.write("\n")
file.write(", ".join([str(B) for B in B_counts]))
file.write("\n")
file.write(", ".join([str(S) for S in S_counts]))