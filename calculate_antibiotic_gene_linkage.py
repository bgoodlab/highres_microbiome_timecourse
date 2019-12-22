import sys
import pylab
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

print "Looking for linkage between antibiotic gene subtypes and core genes"
print "Created with calculate_antibiotic_gene_linkage.py"

pangenome_species = parse_midas_data.parse_pangenome_species()

# long gene map
# long gene = (species, gene) tuple
# id = int
longgene_id_map = {}
id_longgene_map = []


import core_gene_utils

target_total_trajectory_map = {}
target_longgene_observed_trajectory_map = {}
target_longgene_expected_trajectory_map = {}

desired_samples = parse_timecourse_data.morteza_samples
#desired_samples = [parse_timecourse_data.highcoverage_start_1]
#desired_samples = desired_samples[0:2]
#desired_samples = ['6037', '1014.2', '6041']
for sample_idx in xrange(0,len(desired_samples)):

    sample_name = desired_samples[sample_idx]
    
    sys.stderr.write("Processing sample %s...\n" % sample_name)

    longgene_id_counter_map = {}
    total_barcode_count_map = {}

    Btot = 0
    longgene_id_B_map = {}
    
    # create barcode->species map
    sys.stderr.write("Collating species barcodes...\n")
    # first create intermediate data structure:
    # barcode_id->set of longgene ids
    barcode_longgene_ids_map = {} 
    barcode_depth_map = {}
    for species_name in pangenome_species:
        #sys.stderr.write("%s\n" % species_name)
        # Don't use new species yet!
        if species_name=='new_species':
            continue
        
        # Make sure barcodes exist for this timepoint.
        # BG: aside from bugs, shouldn't this always be true? 
        if not barcode_utils.barcodes_exist(species_name, sample_name):
            continue
        
        # Load core genes across those species
        core_longgenes = set()
        core_genes = core_gene_utils.parse_core_genes(species_name)
        #core_genes = core_gene_utils.parse_personal_core_genes(species_name)
        for gene_name in core_genes:
            core_longgenes.add((species_name, gene_name))
         
        # Load barcodes      
        allele_barcode_map = barcode_utils.parse_allele_barcode_tuples(species_name, sample_name)

        for allele in allele_barcode_map.keys():
        
            if allele.endswith('|A') or allele.endswith('|R'):
                # a SNP allele, don't include
                continue
                    
            # form longgene
            longgene = (species_name, allele)
            
            if longgene not in core_longgenes:
                continue
            
            # update longgene id
            if longgene not in longgene_id_map:
                longgene_id_map[longgene] = len(id_longgene_map)
                id_longgene_map.append(longgene)
                
            longgene_id = longgene_id_map[longgene]
    
            if longgene_id not in longgene_id_B_map:
                longgene_id_B_map[longgene_id] = 0
            
    
            for barcode_id, barcode_weight in allele_barcode_map[allele]:
                if barcode_id not in barcode_longgene_ids_map:
                    barcode_longgene_ids_map[barcode_id] = set()
                    barcode_depth_map[barcode_id] = 0
                    
                barcode_longgene_ids_map[barcode_id].add(longgene_id)
                barcode_depth_map[barcode_id] += barcode_weight
                
                Btot += 1
                longgene_id_B_map[longgene_id] += 1
    
    
    longgene_id_q_map = {longgene_id: longgene_id_B_map[longgene_id]*1.0/Btot for longgene_id in longgene_id_B_map}
                
    sys.stderr.write("Done!\n")
    
    if len(barcode_longgene_ids_map)==0:
        continue
    
    
    
    
    sys.stderr.write("Looping through antibiotic genes...\n")
           
    antibiotic_filename = "%s/antibiotic_barcodes/%s_antibiotic_corrected_barcodes.txt.gz" % (config.barcode_directory, sample_name) 
    antibiotic_file = gzip.open(antibiotic_filename,"r")
    antibiotic_file.readline() # header
    
    for line in antibiotic_file:
        line = line.strip()
        items = line.split("\t")
        gene_name = items[0].strip()
        antibiotic_annotation = items[1].strip()
        rpm = float(items[2])
        
        # Only look at tetracycline genes
        if not antibiotic_annotation.startswith('tet'):
            continue
            
        gene_name = antibiotic_annotation
        
        barcode_weight_strs = items[3].strip().split(", ")

        # Only look at ABX genes w/ at least two barcodes
        if len(barcode_weight_strs) < 1.5:
            continue
                                
        num_barcodes = 0
        
        if gene_name not in longgene_id_counter_map:
            longgene_id_counter_map[gene_name] = collections.Counter()
            total_barcode_count_map[gene_name] = 0
            
        for barcode_weight_str in barcode_weight_strs:
            subitems = barcode_weight_str.split(":")
            barcode_id = long(subitems[0])
            
            if barcode_id in barcode_longgene_ids_map:
            
                total_barcode_count_map[gene_name] += 1
                longgene_id_counter_map[gene_name].update( barcode_longgene_ids_map[barcode_id] )


    for gene_name in sorted(longgene_id_counter_map):

        num_barcodes = total_barcode_count_map[gene_name]

        # Only look at ABX genes w/ at least two barcodes that pass
        if num_barcodes < 1.5:
            continue
            
        if gene_name not in target_total_trajectory_map:
            target_total_trajectory_map[gene_name] = numpy.zeros(len(desired_samples))*1.0
            target_longgene_observed_trajectory_map[gene_name] = {}
            target_longgene_expected_trajectory_map[gene_name] = {}
            
        target_total_trajectory_map[gene_name][sample_idx] = num_barcodes
    
        # Tuple of key,counts sorted in descending order of counts    
        longgene_id_count_list = longgene_id_counter_map[gene_name].most_common()
        
        # Print output
        
        longgene_output_strs = []
        for longgene_id, longgene_count in longgene_id_count_list:
        
            longgene_q = longgene_id_q_map[longgene_id]
            # An estimate of a second number.. too small?
            expected_longgene_count = num_barcodes*longgene_q
            
            # Only output things with at least 3 counts
            if longgene_count < 2.5:
                continue
            # And with at least 1% barcode coverage        
            #if longgene_count*1.0/num_barcodes < 0.01:
            #    continue
        
            if longgene_count < 10*expected_longgene_count:
                continue
        
            if longgene_id not in target_longgene_observed_trajectory_map[gene_name]:
                target_longgene_observed_trajectory_map[gene_name][longgene_id] = numpy.zeros(len(desired_samples))*1.0
                target_longgene_expected_trajectory_map[gene_name][longgene_id] = numpy.zeros(len(desired_samples))*1.0
            
            target_longgene_observed_trajectory_map[gene_name][longgene_id][sample_idx] = longgene_count   
            target_longgene_expected_trajectory_map[gene_name][longgene_id][sample_idx] = expected_longgene_count   
            
                
sys.stderr.write("Done looping over samples!\n")   

# Now print output
for gene_name in sorted(target_total_trajectory_map):
    print gene_name 
    print " ".join([str(d) for d in target_total_trajectory_map])
    
    for longgene_id in sorted(target_longgene_observed_trajectory_map[gene_name]):
        longgene_items = id_longgene_map[longgene_id]
        longgene_str = "%s|%s" % longgene_items
        
        observed_trajectory = target_longgene_observed_trajectory_map[gene_name][longgene_id]
        expected_trajectory = target_longgene_expected_trajectory_map[gene_name][longgene_id]
        
        # Only print things with at least two surprising events!
        if (observed_trajectory>10*expected_trajectory).sum() < 0.5*len(expected_trajectory):
            continue
        
        print longgene_str
        print " ".join([str(d) for d in observed_trajectory])
        print " ".join([str(d) for d in expected_trajectory])
        
    print "---"            
sys.stderr.write("Done!\n")
