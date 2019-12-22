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

pangenome_species = parse_midas_data.parse_pangenome_species()

# long gene map
# long gene = (species, gene) tuple
# id = int
longgene_id_map = {}
id_longgene_map = []

longgene_id_counter_map = {}
total_barcode_count_map = {}

#desired_samples = parse_timecourse_data.morteza_samples
#desired_samples = [parse_timecourse_data.highcoverage_start_1]
#desired_samples = desired_samples[0:2]
desired_samples = [parse_timecourse_data.highcoverage_antibiotic]

longgene_antibiotic_map = {}
antibiotic_barcode_count_map = {}
all_barcode_ids = set()

for sample_name in desired_samples:

    sys.stderr.write("Processing sample %s...\n" % sample_name)
    
    # create barcode->species map
    sys.stderr.write("Collating antibiotic barcodes...\n")
           
    antibiotic_filename = "%s/antibiotic_barcodes/%s_antibiotic_corrected_barcodes.txt.gz" % (config.barcode_directory, sample_name) 
    antibiotic_file = gzip.open(antibiotic_filename,"r")
    antibiotic_file.readline() # header
    barcode_antibiotic_map = {}
    
    total_barcodes = 0
    
    for line in antibiotic_file:
        line = line.strip()
        items = line.split("\t")
        gene_name = items[0].strip()
        antibiotic_annotation = items[1].strip()
        #antibiotic = antibiotic_annotation.split("__")[0]
        antibiotic = "%s:%s" % (antibiotic_annotation, gene_name)
        rpm = float(items[2])
        
        barcode_weight_strs = items[3].strip().split(", ")
    
        for barcode_weight_str in barcode_weight_strs:
            subitems = barcode_weight_str.split(":")
            barcode_id = long(subitems[0])
            all_barcode_ids.add(barcode_id)
            if barcode_id not in barcode_antibiotic_map:
                barcode_antibiotic_map[barcode_id] = set()
            barcode_antibiotic_map[barcode_id].add(antibiotic)
        
        if antibiotic not in antibiotic_barcode_count_map:    
            antibiotic_barcode_count_map[antibiotic] = 0
        antibiotic_barcode_count_map[antibiotic] += len(barcode_weight_strs)
        
    sys.stderr.write("Loaded %d barcodes!\n" % len(barcode_antibiotic_map))                
    sys.stderr.write("Looping through species barcodes...\n")
    for species_name in pangenome_species[:1]:
        
        # Don't use new species yet!
        if species_name=='new_species':
            continue
        
        # Make sure barcodes exist for this timepoint.
        # BG: aside from bugs, shouldn't this always be true? 
        if not barcode_utils.barcodes_exist(species_name, sample_name):
            continue
         
        # Load barcodes      
        allele_barcode_map = barcode_utils.parse_allele_barcode_tuples(species_name, sample_name)

        for allele in allele_barcode_map.keys():
        
            if not (allele.endswith('|A') or allele.endswith('|R')):
                # a SNP allele, don't include
                continue
        
            if len(allele_barcode_map)==0:
                continue
            
            # form longgene
            longgene = (species_name, allele)
            
            for barcode_id, barcode_weight in allele_barcode_map[allele]:
                
                all_barcode_ids.add(barcode_id)
                if barcode_id in barcode_antibiotic_map:
                    # Got a match!
                    
                    if longgene not in longgene_antibiotic_map:
                        longgene_antibiotic_map[longgene] = collections.Counter()
                        total_barcode_count_map[longgene] = len(allele_barcode_map[allele])
                    
                    longgene_antibiotic_map[longgene].update( barcode_antibiotic_map[barcode_id])
                    
# Time to print!            
for longgene in sorted(longgene_antibiotic_map):

    species, gene_name = longgene
    
    num_barcodes = total_barcode_count_map[longgene]

    # Only look at ABX genes w/ at least two barcodes that pass
    if num_barcodes < 1.5:
        continue
        
    antibiotic_count_list = longgene_antibiotic_map[longgene].most_common()
        
    # Print output
    antibiotic_output_strs = []
    for antibiotic, antibiotic_count in antibiotic_count_list:
        # Only output things with at least 2 counts
        if antibiotic_count < 2.5:
            continue
                
        if antibiotic_count*1.0/num_barcodes < 0.05:
            continue
                
        antibiotic_output_strs.append("|".join([antibiotic, str(antibiotic_count)]))
     
    if len(antibiotic_output_strs) < 0.5:
        continue
        
    antibiotic_output_str = ", ".join(antibiotic_output_strs)
    output_str = "\t".join([species, gene_name, "%d" % num_barcodes, antibiotic_output_str])
    print output_str    

print len(all_barcode_ids), antibiotic_barcode_count_map           
sys.stderr.write("Done!\n")
