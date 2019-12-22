###########################################################################
#
# This script calculates the snp-snp .p files that are ultimately used to make Fig. 4 
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

import os.path
import config
import cPickle as pickle
from math import ceil
import sys
import argparse
from scipy.special import gammaln as loggamma
from scipy.special import gammainc
from scipy.special import expm1
from math import log,exp
from scipy.stats import poisson

from numpy.random import shuffle
from stats_utils import calculate_poisson_logP

pangenome_species = parse_midas_data.parse_pangenome_species()

    
###########################################################################
#
# Standard header to read in argument information
#
###########################################################################

min_S = 11.5 # (this means later we can look for things where 4th gamete is at least 3)
Pstar = 0.1

gamete_idx_map = barcode_utils.gamete_idx_map
allele_idx_map = barcode_utils.allele_idx_map

# shortsnp = (species, contig, location)
shortsnp_id_map = {}
id_shortsnp_map = []

# long snp = (shortsnp_id, allele) 
longsnp_id_map = {}
id_longsnp_map = []


longsnp_shortsnp_map = [] # from longsnp_id to shortsnp_id

desired_speciess = sys.argv[1:]

#desired_samples = parse_timecourse_data.highcoverage_samples
#desired_samples = [parse_timecourse_data.highcoverage_end]
desired_samples = parse_timecourse_data.morteza_samples
#desired_samples = ['6037','6037.3']
#desired_samples = ['6037.3'] 
#for testing
#desired_samples = parse_timecourse_data.morteza_samples[1:3]
#desired_samples = [parse_timecourse_data.highcoverage_antibiotic, parse_timecourse_data.highcoverage_end]
#desired_samples = [parse_timecourse_data.highcoverage_start_2, parse_timecourse_data.highcoverage_antibiotic, parse_timecourse_data.highcoverage_postantibiotic, parse_timecourse_data.highcoverage_end]

G_function_map = barcode_utils.parse_G_function_map()
total_barcode_map = barcode_utils.parse_total_barcode_map()

Btots = numpy.array([total_barcode_map[sample_name][0] for sample_name in desired_samples])*1.0

shortsnp_Bt_matrix_map = {} # shortsnp x 2 x T
shortsnp_Gt_matrix_map = {} # shortsnp_id x 2 x T 

shortsnp_S_matrix_map = {} # map from (shortsnp_id, tuple) to 2x2

shortsnp_error_map = {} # map from shortsnp id to 2 (# of errors)

longsnp_counter_map = []

tracked_shortsnps = barcode_utils.parse_tracked_snps(set(desired_speciess))
for shortsnp in tracked_shortsnps:

    if shortsnp not in shortsnp_id_map:
        shortsnp_id = len(id_shortsnp_map)
        shortsnp_id_map[shortsnp] = shortsnp_id
        id_shortsnp_map.append(shortsnp)
        shortsnp_Bt_matrix_map[shortsnp_id] = numpy.zeros((2,len(desired_samples)))*1.0
        shortsnp_Gt_matrix_map[shortsnp_id] = numpy.zeros((2,len(desired_samples)))*1.0
        shortsnp_error_map[shortsnp_id] = numpy.array([0,0])
        
# Counts # of shared barcodes between shortsnps 
# (for efficiency, we actually calculate an approximate versin, ignoring multiple
#  alleles of same snp in a single barcode. This gets caught in the next step)
initial_S_matrix = numpy.zeros((len(id_shortsnp_map),len(id_shortsnp_map))) 
    
for sample_idx in xrange(0,len(desired_samples)):
    # Process each sample separately
    # (combine results at end)
    
    sample_name = desired_samples[sample_idx]
        
    G_function = G_function_map[sample_name]
    Btot = total_barcode_map[sample_name][0]
       
    sys.stderr.write("Processing sample %s...\n" % sample_name)
    #sys.stderr.write("Loading depth map...\n")
    barcode_depth_map = barcode_utils.parse_barcode_depth_map(sample_name)
    #sys.stderr.write("Done!\n")
        
    # create barcode->species map
    sys.stderr.write("Collating species barcodes...\n")
        
    longsnp_barcode_map = {} 
    barcode_shortsnp_map = {}

    for species_name in desired_speciess:
        
        # Make sure barcodes exist for this species at this timepoint.
        if not barcode_utils.barcodes_exist(species_name, sample_name):
            continue
         
        # Load barcodes      
        allele_barcode_map, allele_error_map = barcode_utils.parse_allele_barcode_tuples(species_name, sample_name)
        
        for allele in allele_barcode_map:
        
            if not (allele[-1]=='A' or allele[-1]=='R'):
                # Not a snp, ignore
                continue
                    
            snp_str = allele[:-2]
            polarization = allele[-1]
            allele_idx = allele_idx_map[polarization]
            
            items = snp_str.split("|")
            contig = items[0]
            location = long(items[1])
            shortsnp = (species_name, contig, location)
            
            if shortsnp not in shortsnp_id_map:
            	continue
            	 
            shortsnp_id = shortsnp_id_map[shortsnp]
            
            longsnp = (shortsnp_id, allele_idx)
            
            if longsnp not in longsnp_id_map:
                longsnp_id = len(id_longsnp_map)
                id_longsnp_map.append(longsnp)
                longsnp_id_map[longsnp] = longsnp_id
                longsnp_shortsnp_map.append(shortsnp_id)
                longsnp_counter_map.append( collections.Counter() )
                
            longsnp_id = longsnp_id_map[longsnp]
                                       
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
                
                shortsnp_Bt_matrix_map[shortsnp_id][ allele_idx, sample_idx] += 1
                shortsnp_Gt_matrix_map[shortsnp_id][ allele_idx, sample_idx] += G
                        
                if longsnp_id not in longsnp_barcode_map:
                    longsnp_barcode_map[longsnp_id] = set()
                        
                # Add barcode to list
                longsnp_barcode_map[longsnp_id].add( barcode_id)
                    
                if barcode_id not in barcode_shortsnp_map:
                    barcode_shortsnp_map[barcode_id] = []
                
                barcode_shortsnp_map[barcode_id].append(shortsnp_id)
                
                if barcode_weight > 1.5:
                	# Only count as a non-error if we had the capacity
                	# to observe an error in the first place
                	shortsnp_error_map[shortsnp_id][1] += 1
            
            # Add in barcode errors
            for barcode_id, barcode_weight in allele_error_map[allele]:
                
                # don't include barcode if too little coverage
                # Poor man's error correction
                if barcode_id not in barcode_depth_map: 
                    continue
                
                D = barcode_depth_map[barcode_id]
                G = G_function(D)
                
                if G>10:
                    continue   
                    
                shortsnp_error_map[shortsnp_id][0] += 1
                
    sys.stderr.write("Done!\n") 
    sys.stderr.write("Pruning barcode->shortsnp entries...\t")
    for barcode_id in barcode_shortsnp_map:
    	
    	unique_elements, counts_elements = numpy.unique(barcode_shortsnp_map[barcode_id], return_counts=True)           
    	
    	barcode_shortsnp_map[barcode_id] = (unique_elements, counts_elements)
    sys.stderr.write("Done!\n")
    sys.stderr.write("Counting shared barcodes for sites!\t")                
    for focal_longsnp_id in longsnp_barcode_map:
        focal_shortsnp_id = longsnp_shortsnp_map[focal_longsnp_id]
        for barcode_id in longsnp_barcode_map[focal_longsnp_id]: 
        	# crude version
        	unique_elements, counts_elements = barcode_shortsnp_map[barcode_id]
        	initial_S_matrix[focal_shortsnp_id, unique_elements] += counts_elements
        	
    sys.stderr.write("Done!\n")
    
    sys.stderr.write("Tallying pairs...\t")
    total_pairs = initial_S_matrix.shape[0]*initial_S_matrix.shape[1]
    sys.stderr.write("%0.1fe06\n" % (total_pairs*1e-06))
                
sys.stderr.write("Done!\n")

sys.stderr.write("Calculating Qs...\t")
shortsnp_Qt_matrix_map = {}
shortsnp_B_matrix_map = {}
for shortsnp_id in xrange(0,len(id_shortsnp_map)):
    shortsnp_Qt_matrix_map[shortsnp_id] =  shortsnp_Bt_matrix_map[shortsnp_id]*1.0/Btots[None,:]
    shortsnp_B_matrix_map[shortsnp_id] =  shortsnp_Bt_matrix_map[shortsnp_id].sum(axis=1)
    

sys.stderr.write("Done!\n")

sys.stderr.write("Calculating significantly linked sites...\n")


# Now go through and check for significance
total_shortsnps = len(id_shortsnp_map)
log_corrected_Pstar = log(Pstar/total_shortsnps/total_shortsnps)

num_threshold_rejected = 0
num_model_rejected = 0
num_passed = 0
num_bad_coverage = 0
for focal_shortsnp_id in xrange(0,len(id_shortsnp_map)):
    
    B = shortsnp_B_matrix_map[focal_shortsnp_id].sum()
    target_shortsnp_ids = numpy.nonzero(initial_S_matrix[focal_shortsnp_id,:]>min_S)[0]
    
    if B<=min_S:
        
        num_bad_coverage+=len(target_shortsnp_ids)
        continue
    
    shortsnp_S_matrix_map[focal_shortsnp_id] = {}
    
    for target_shortsnp_id in target_shortsnp_ids:
        
        S = initial_S_matrix[focal_shortsnp_id, target_shortsnp_id]
    
        if S<=min_S:
        	# Should never get here anymore
            num_threshold_rejected += 1
        else:
        
            Qjbt = shortsnp_Qt_matrix_map[target_shortsnp_id]
            Giat = shortsnp_Gt_matrix_map[focal_shortsnp_id]
        
            S0s = numpy.einsum('at,bt', Giat, Qjbt)
            S0 = S0s.sum()
        
            logP = calculate_poisson_logP(S,S0)
                    
            if (logP < log_corrected_Pstar) or (focal_shortsnp_id==target_shortsnp_id):        
                num_passed += 1
                
                shortsnp_S_matrix_map[ focal_shortsnp_id][target_shortsnp_id] = numpy.zeros((2, 2,2))
                shortsnp_S_matrix_map[focal_shortsnp_id][target_shortsnp_id][1] = S0s
            	
            	#print S,S0,logP
            	    
            else:
                num_model_rejected += 1
                
        if target_shortsnp_id==focal_shortsnp_id:
            if target_shortsnp_id not in shortsnp_S_matrix_map[focal_shortsnp_id]:
                sys.stderr.write("Problem! %g, %g\n" % (S,S0))

del initial_S_matrix
sys.stderr.write("Done!\n")
sys.stderr.write("%d bad coverage rejected, %d threshold rejected, %d model rejected, %d passed!\n" % (num_bad_coverage, num_threshold_rejected, num_model_rejected, num_passed))

# Now we have to take another pass to fill in per gamete Ss for shortsnps that are significant

sys.stderr.write("Second pass...\n")
for sample_idx in xrange(0,len(desired_samples)):
    # Process each sample separately
    # (combine results at end)
    
    sample_name = desired_samples[sample_idx]    
    sys.stderr.write("Processing sample %s...\n" % sample_name)
    #sys.stderr.write("Loading depth map...\n")
    barcode_depth_map = barcode_utils.parse_barcode_depth_map(sample_name)
    #sys.stderr.write("Done!\n")
        
    if len(barcode_depth_map)==0:
        continue

    # create barcode->species map
    sys.stderr.write("Collating species barcodes...\n")
        
    longsnp_barcode_map = {} 
    barcode_longsnp_map = {}

    G_function = G_function_map[sample_name]

    for species_name in desired_speciess:
        
        # Make sure barcodes exist for this species at this timepoint.
        if not barcode_utils.barcodes_exist(species_name, sample_name):
            sys.stderr.write("Barcodes don't exist for %s, %s\n" % (species_name, sample_name))
            continue
         
        # Load barcodes     
        allele_barcode_map, allele_error_map = barcode_utils.parse_allele_barcode_tuples(species_name, sample_name)

        for allele in allele_barcode_map:
        
            if not (allele[-1]=='A' or allele[-1]=='R'):
                # Not a snp, ignore
                continue
                    
            snp_str = allele[:-2]
            polarization = allele[-1]
            allele_idx = allele_idx_map[polarization]
            
            items = snp_str.split("|")
            contig = items[0]
            location = long(items[1])
            shortsnp = (species_name, contig, location)
            
            if shortsnp not in shortsnp_id_map:
            	continue
            	
            shortsnp_id = shortsnp_id_map[shortsnp]
            longsnp = (shortsnp_id, allele_idx)
            longsnp_id = longsnp_id_map[longsnp]
                                       
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
                        
                if longsnp_id not in longsnp_barcode_map:
                    longsnp_barcode_map[longsnp_id] = [] 
                                            
                # Add barcode to list
                longsnp_barcode_map[longsnp_id].append( barcode_id)
                    
                if barcode_id not in barcode_longsnp_map:
                    barcode_longsnp_map[barcode_id] = []
                
                barcode_longsnp_map[barcode_id].append( longsnp_id)
                
    sys.stderr.write("Done!\n")            
    # Make sure we loaded some barcodes...
    if len(barcode_longsnp_map)==0:
        continue
    
    sys.stderr.write("Counting shared barcodes for sites!\t")                
    for focal_longsnp_id in longsnp_barcode_map:
    
        focal_shortsnp_id, focal_allele = id_longsnp_map[focal_longsnp_id]
        
        if focal_shortsnp_id not in shortsnp_S_matrix_map:
            continue
        
        longsnp_counts = collections.Counter()
        
        for barcode_id in longsnp_barcode_map[focal_longsnp_id]: 
            longsnp_counts.update( barcode_longsnp_map[barcode_id])
            
        for target_longsnp_id in longsnp_counts:
            
            S = longsnp_counts[target_longsnp_id]
            
            target_shortsnp_id, target_allele = id_longsnp_map[target_longsnp_id]
            
            if target_shortsnp_id in shortsnp_S_matrix_map[focal_shortsnp_id]:
                
                gamete = (focal_allele, target_allele)
                
                shortsnp_S_matrix_map[focal_shortsnp_id][target_shortsnp_id][0][gamete] += S
            
    sys.stderr.write("Done!\n")
  
sys.stderr.write("Done!\n")

data_to_save = {'id_longsnp_map': id_shortsnp_map, 'allele_idx_map': allele_idx_map, 'B': shortsnp_B_matrix_map, 'S': shortsnp_S_matrix_map, 'E': shortsnp_error_map}

output_filename = "%s/new_barcode_trajectories/%s_snp_barcodes.p" % (config.barcode_directory, "_".join(desired_speciess))

sys.stderr.write("Saving data...\t")
pickle.dump( data_to_save, open( output_filename, "wb" ) )
sys.stderr.write("Done!\n")

