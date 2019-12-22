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
#desired_samples = ['6037.3']
#desired_samples = [parse_timecourse_data.highcoverage_start_1, parse_timecourse_data.highcoverage_antibiotic]
#desired_samples = ['4026.2']
min_barcode_depth = 10

print "sample_name unfiltered_total_barcodes unfiltered_total_barcode_depth total_barcodes total_barcode_depth"
for sample_name in desired_samples:

    sys.stderr.write("Processing sample %s...\n" % sample_name)
    
    sys.stderr.write("Loading error corrected barcodes...\n")
    raw_barcode_error_map, raw_barcode_depth_map = barcode_utils.parse_raw_barcode_error_correction_map(sample_name)

    
    if len(raw_barcode_error_map)==0:
        sys.stderr.write("Couldn't load raw barcode error map for %s!\n" % sample_name) 
        continue
    sys.stderr.write("Done!\n")
  
    sys.stderr.write("Correcting barcodes...\n")
    corrected_barcode_depth_map = {}
    for barcode_str in raw_barcode_error_map:
        corrected_barcode_str = raw_barcode_error_map[barcode_str]
        
        if corrected_barcode_str not in corrected_barcode_depth_map:
            corrected_barcode_depth_map[corrected_barcode_str] = 0
            
            corrected_barcode_depth_map[corrected_barcode_str] += raw_barcode_depth_map[corrected_barcode_str]
        
    sys.stderr.write("Calculating barcode blacklist...\n")
    corrected_barcode_str_blacklist = barcode_utils.calculate_barcode_blacklist(corrected_barcode_depth_map.iterkeys())
    sys.stderr.write("Done!\n")
    
    total_barcodes = 0
    total_barcode_depth = 0
    
    unfiltered_total_barcodes = 0
    unfiltered_total_barcode_depth = 0
    
    for corrected_barcode_str in corrected_barcode_depth_map:
        
        
        depth = corrected_barcode_depth_map[corrected_barcode_str]
        
        unfiltered_total_barcodes+=1
        unfiltered_total_barcode_depth+=depth
        
        if corrected_barcode_str in corrected_barcode_str_blacklist:
            continue
            
        if depth<10:
            continue
            
        if depth>9999:
            continue
        
        total_barcodes += 1
        total_barcode_depth += depth
        
    print sample_name, unfiltered_total_barcodes, unfiltered_total_barcode_depth, total_barcodes, total_barcode_depth
    
    