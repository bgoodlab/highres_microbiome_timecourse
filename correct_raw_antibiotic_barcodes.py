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
#desired_samples = ['6037.3']
print sorted(desired_samples)
for sample_name in sorted(desired_samples):

    sys.stderr.write("Processing sample %s...\n" % sample_name)
    
    sys.stderr.write("Loading barcode error map...\n")
    raw_barcode_error_map, raw_barcode_depth_map = barcode_utils.parse_raw_barcode_error_correction_map(sample_name)
    sys.stderr.write("Done!\n")
    
    del raw_barcode_depth_map
    
    # Now load corrected barcodes that got mapped to. 
    sys.stderr.write("Loading all mapped barcodes...\n")
    barcode_depth_map, barcode_id_map = barcode_utils.parse_all_mapped_barcodes(sample_name)
    sys.stderr.write("Done!\n")
    
    # We can now find the corrected id by going from raw_barcode_error_map
    # to barcode_id_map
    
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
            
            corrected_barcode_str = raw_barcode_error_map[barcode_str]
            
            if corrected_barcode_str in barcode_id_map:
                
                barcode_id = barcode_id_map[corrected_barcode_str]
                
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
    
    
    