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

#desired_samples = parse_timecourse_data.morteza_samples
#desired_samples = [parse_timecourse_data.highcoverage_start_1]
#desired_samples = desired_samples[0:2]
desired_samples = [parse_timecourse_data.highcoverage_antibiotic]
for sample_name in desired_samples:

    sys.stderr.write("Processing sample %s...\n" % sample_name)
    
    # create barcode->species map
    sys.stderr.write("Collating antibiotic barcodes...\n")
    
    antibiotic_filename = "%s/antibiotic_barcodes/%s_antibiotic_corrected_barcodes.txt.gz" % (config.barcode_directory, sample_name) 
    antibiotic_file = gzip.open(antibiotic_filename,"r")
    antibiotic_file.readline() # header
    
    barcode_antibiotic_map = {}
    
    
    for line in antibiotic_file:
        line = line.strip()
        items = line.split("\t")
        gene_name = items[0].strip()
        antibiotic_annotation = items[1].strip()
        rpm = float(items[2])
        #antibiotic = antibiotic_annotation.split("__")[0]
        antibiotic = "%s:%s" % (antibiotic_annotation, gene_name)
        barcode_weight_strs = items[3].strip().split(", ")

        # Only look at ABX genes w/ at least two barcodes
        if len(barcode_weight_strs) < 1.5:
            continue
                                
        for barcode_weight_str in barcode_weight_strs:
            subitems = barcode_weight_str.split(":")
            barcode_id = long(subitems[0])
            
            if barcode_id not in barcode_antibiotic_map:
                barcode_antibiotic_map[barcode_id] = set()
                
            barcode_antibiotic_map[barcode_id].add(antibiotic)
    antibiotic_file.close()
    sys.stderr.write("Done!\n")
    
    sys.stderr.write("Looping through antibiotic genes...\n")
    
    output_str = "\t".join(["Gene name", "Antibiotic annotation", "RPM", "Total num Barcodes", "linked antibiotic class|num barcodes"])
    print output_str
           
    antibiotic_filename = "%s/antibiotic_barcodes/%s_antibiotic_corrected_barcodes.txt.gz" % (config.barcode_directory, sample_name) 
    antibiotic_file = gzip.open(antibiotic_filename,"r")
    antibiotic_file.readline() # header
    
    for line in antibiotic_file:
        line = line.strip()
        items = line.split("\t")
        gene_name = items[0].strip()
        antibiotic_annotation = items[1].strip()
        rpm = float(items[2])
        
        barcode_weight_strs = items[3].strip().split(", ")

        # Only look at ABX genes w/ at least two barcodes
        if len(barcode_weight_strs) < 1.5:
            continue
                                
        num_barcodes = 0
        longgene_id_counter = collections.Counter()
        for barcode_weight_str in barcode_weight_strs:
            subitems = barcode_weight_str.split(":")
            barcode_id = long(subitems[0])
            
            if barcode_id in barcode_antibiotic_map:
                num_barcodes += 1
                longgene_id_counter.update( barcode_antibiotic_map[barcode_id] )
        
        # Only look at ABX genes w/ at least two barcodes that pass
        if num_barcodes < 1.5:
            continue
        
        longgene_id_count_list = longgene_id_counter.most_common()
        
        # Print output

        
        longgene_output_strs = []
        for longgene_id, longgene_count in longgene_id_count_list:
            # Only output things with at least 2 counts
            if longgene_count < 2.5:
                continue
                
            if longgene_count*1.0/num_barcodes < 0.05:
                continue
                
            antibiotic_name = longgene_id
            longgene_output_strs.append("|".join([antibiotic_name, str(longgene_count)]))
        
        # Only print things linked to more than one class
        if len(longgene_output_strs)<1.5 and (not antibiotic_annotation.startswith('beta')):
            continue
        
        
        longgene_output_str = ", ".join(longgene_output_strs[1:])
        output_str = "\t".join([gene_name, antibiotic_annotation, "%g" % rpm, "%d" % num_barcodes, longgene_output_str])
        print output_str    
            
    sys.stderr.write("Done!\n")
