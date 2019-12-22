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
import os
import config
import cPickle as pickle

desired_samples = parse_timecourse_data.morteza_samples
#desired_samples = [parse_timecourse_data.highcoverage_start_1, parse_timecourse_data.highcoverage_antibiotic]
for sample_name in desired_samples:

    sys.stderr.write("Processing sample %s...\n" % sample_name)
    
    sys.stderr.write("Loading error corrected barcodes...\n")
    barcode_error_map = barcode_utils.parse_barcode_error_correction_map(sample_name)
    sys.stderr.write("Done!\n")
  
    if len(barcode_error_map)==0:
        continue
    
    # First correct the all_barcodes file
    sys.stderr.write("Correcting barcode strs...\n")
    barcode_filename = "%s%s/output/all_barcodes.gz" % (config.barcode_directory, sample_name)    
    # Open uncorrected barcode file
    barcode_file = gzip.GzipFile(barcode_filename,"r")
    # Read header
    line = barcode_file.readline() # skip header
      
    
    # dictionary of barcode_string -> corrected_barcode_id 
    corrected_barcode_id_map = {} 
    # For each uncorrected barcode 
    for line in barcode_file:
        items = line.split()
        barcode_id = long(items[0])
        barcode_str = items[1]
        barcode_count = long(items[2])
        
        if barcode_id not in barcode_error_map:
            # This barcode error corrects to itself
            corrected_barcode_id = barcode_id
            
        else:
            # This barcode error corrects to something else
            corrected_barcode_id = barcode_error_map[barcode_id]
            
        corrected_barcode_id_map[barcode_str] = corrected_barcode_id
             
    barcode_file.close()
    sys.stderr.write("Done!\n")
    
    sys.stderr.write("Correcting antibiotic reads...\n")
    
    antibiotic_filename = "%s/antibiotic_barcodes/%s_antibiotic_barcodes.txt.gz" % (config.barcode_directory, sample_name)  
    new_antibiotic_filename = "%s/antibiotic_barcodes/%s_antibiotic_corrected_barcodes.txt.gz" % (config.barcode_directory, sample_name) 
    
    antibiotic_file = gzip.open(antibiotic_filename,"r")
    new_antibiotic_file = gzip.open(new_antibiotic_filename,"w")
    
    header_line = "\t".join(['Gene name', 'ABX Annotation', 'RPM', 'Corrected Barcode IDs'])
    new_antibiotic_file.write("%s\n" % header_line)
    
    for line in antibiotic_file:
        items = line.split("\t")
        gene_name = items[0].strip()
        subtype = items[1].strip()
        rpm = float(items[2])
        barcode_strs = items[3].strip().split(",")
        
        barcode_weight_map = {}
        for barcode_str in barcode_strs:
            
            if barcode_str in corrected_barcode_id_map:
                barcode_id = corrected_barcode_id_map[barcode_str]
                
                if barcode_id not in barcode_weight_map:
                    barcode_weight_map[barcode_id] = 0
                
                barcode_weight_map[barcode_id] += 1
        
        if len(barcode_weight_map) < 1.5:
            continue
                    
        barcode_weight_str = ", ".join(["%d:%d" % (barcode_id,barcode_weight_map[barcode_id]) for barcode_id in sorted(barcode_weight_map)])
        
        output_str = "\t".join([gene_name, subtype, "%g" % rpm, barcode_weight_str])
        new_antibiotic_file.write("%s\n" % output_str)
        
        
    antibiotic_file.close()
    new_antibiotic_file.close()
    sys.stderr.write("Done!\n")
    
    
    