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
#desired_samples = ['4026.2']
min_barcode_depth = 10
for sample_name in desired_samples:

    sys.stderr.write("Processing sample %s...\n" % sample_name)
    
    sys.stderr.write("Loading error corrected barcodes...\n")
    raw_barcode_error_map, raw_barcode_depth_map = barcode_utils.parse_raw_barcode_error_correction_map(sample_name)

    
    if len(raw_barcode_error_map)==0:
        sys.stderr.write("Couldn't load raw barcode error map for %s!\n" % sample_name) 
        continue
    sys.stderr.write("Done!\n")
  
    sys.stderr.write("Calculating barcode blacklist...\n")
    corrected_barcode_str_blacklist = barcode_utils.calculate_barcode_blacklist(raw_barcode_depth_map.iterkeys())
    sys.stderr.write("Done!\n")
    
    # First correct the all_barcodes file
    sys.stderr.write("Correcting all barcodes file...\n")
    barcode_filename = "%s%s/output/all_barcodes.gz" % (config.barcode_directory, sample_name)
    new_barcode_filename = "%s%s/output/all_corrected_barcodes.gz" % (config.barcode_directory, sample_name)
        
    # Open uncorrected barcode file
    barcode_file = gzip.GzipFile(barcode_filename,"r")
    # Read header
    header_line = barcode_file.readline() # skip header
      
    
    # map from barcode str to barcode id
    barcode_id_map = {} 
    # map from barcode_id to barcode_str
    barcode_str_map = {}
    sys.stderr.write("Loading all barcodes file...\n")
    for line in barcode_file:
        items = line.split()
        barcode_id = long(items[0])
        barcode_str = items[1].strip()
        barcode_count = long(items[2])
        
        if barcode_id not in barcode_str_map:
        
            if barcode_str not in raw_barcode_error_map:
                print "Weird thing:", barcode_id, barcode_str, len(barcode_str)
                sys.exit(1)
        
            barcode_str_map[barcode_id] = barcode_str
            barcode_id_map[barcode_str] = barcode_id
            
        
        
    # Done
    barcode_file.close()
    sys.stderr.write("Done!\n")
    
    # Now calculate barcode_error_map (# ID to corrected ID)
    corrected_barcode_str_map = []
    corrected_barcode_id_map = {}
    barcode_error_map = {}
    corrected_barcode_id_blacklist = set()
    for barcode_id in barcode_str_map:
        
        barcode_str = barcode_str_map[barcode_id]
        corrected_barcode_str = raw_barcode_error_map[barcode_str]
        
        if corrected_barcode_str not in corrected_barcode_id_map:
        
            # Create a new entry for that corrected barcode str
        
            corrected_barcode_id = len(corrected_barcode_str_map)
            corrected_barcode_str_map.append(corrected_barcode_str)
            corrected_barcode_id_map[corrected_barcode_str] = corrected_barcode_id
            
            if corrected_barcode_str in corrected_barcode_str_blacklist:
                corrected_barcode_id_blacklist.add(corrected_barcode_id)
            
        corrected_barcode_id = corrected_barcode_id_map[corrected_barcode_str]
        barcode_error_map[barcode_id] = corrected_barcode_id
        
        #if corrected_barcode_str not in barcode_id_map:
        #    # corrected barcode wasn't assigned a barcode
        #    sys.stderr.write("Weird thing happened! %s %s\n" % (barcode_str, corrected_barcode_str))
            
        #    barcode_id_map[corrected_barcode_str] = barcode_id
        #    barcode_str_map[barcode_id] = corrected_barcode_str
            
    
    # Now calculate:
    # dictionary of corrected_barcode_id -> [count, corrected_barcode_str (or "" if not set yet)] 
    #corrected_barcode_stats = {}
    # map from corrected_barcode_id -> count
    corrected_barcode_depth_map = {} 
    for barcode_id in barcode_error_map:
        
        corrected_barcode_id = barcode_error_map[barcode_id]
        if corrected_barcode_id not in corrected_barcode_depth_map:
            corrected_barcode_str = corrected_barcode_str_map[corrected_barcode_id]
            corrected_barcode_depth = raw_barcode_depth_map[corrected_barcode_str]
            corrected_barcode_depth_map[corrected_barcode_id] = corrected_barcode_depth
            
            if corrected_barcode_depth < min_barcode_depth:
                corrected_barcode_id_blacklist.add(corrected_barcode_id)
        
        
    del raw_barcode_depth_map
    del raw_barcode_error_map
    del corrected_barcode_str_blacklist
    
    # From here on, same as previous file!
    # now print stuff out!
    # Create corrected barcode file  
    sys.stderr.write("Writing new barcode counts file...\n")
    new_barcode_file = gzip.GzipFile(new_barcode_filename,"w")
    # Write header
    new_barcode_file.write(header_line)
    for barcode_id in sorted(corrected_barcode_depth_map):
        
        barcode_counts = corrected_barcode_depth_map[barcode_id]
        barcode_str = corrected_barcode_str_map[barcode_id]
        
        # Don't output the bad ones!
        if barcode_id in corrected_barcode_id_blacklist:
            continue
            
        output_str = "%s\t%s\t%s\n" % (barcode_id, barcode_str, barcode_counts)
        new_barcode_file.write(output_str)
    
    new_barcode_file.close()        
    sys.stderr.write("Done!\n")
    
         
    # Now correct barcodes for each species separately
    for filename in os.listdir('%s%s/output' % (config.barcode_directory, sample_name)):
    
        barcode_filename = "%s%s/output/%s" % (config.barcode_directory, sample_name, filename)
        
        if not barcode_filename.endswith('.barcodes.gz'):
            continue
        
        sys.stderr.write("Correcting %s...\n" % filename)  
            
        filename_items = barcode_filename.split(".")
        new_barcode_filename = ".".join(filename_items[:-2]+["corrected_barcodes", "gz"])
        
        # Open uncorrected barcode file
        barcode_file = gzip.GzipFile(barcode_filename,"r")
        # Read header
        line = barcode_file.readline() 
        
        # Create corrected barcode file
        new_barcode_file = gzip.GzipFile(new_barcode_filename,"w")
        # Write header
        new_barcode_file.write(line)
        
        # Loop through each tracked allele in the uncorrected barcode file
        for line in barcode_file:
            line = line.strip()
            items = line.split("\t")
            allele = items[0].strip()
            
            if len(items)>1:
        
                # If there are barcodes that mapped to this allele
                # run error correction
        
                barcode_weights = {}
                
                # Parse the barcodes
                barcode_items = items[1].split(",")
                for barcode_item in barcode_items:
                
                    # Break barcode into barcode_id and # reads ("weight")
                    barcode_subitems = barcode_item.split(":")
                    original_barcode_id = long(barcode_subitems[0])
                    barcode_weight = long(barcode_subitems[1])
                
                    # Error correct if necessary
                    corrected_barcode_id = barcode_error_map[original_barcode_id]
                    
                    # Don't output the bad ones!
                    if corrected_barcode_id in corrected_barcode_id_blacklist:
                        continue
                    
                    # Add entry to barcode_weights dictionary
                    if corrected_barcode_id not in barcode_weights:
                        barcode_weights[corrected_barcode_id] = 0
                     
                    barcode_weights[corrected_barcode_id] += barcode_weight
                
                barcode_weight_str = ", ".join(["%d:%d" % (corrected_barcode_id, barcode_weights[corrected_barcode_id]) for corrected_barcode_id in barcode_weights])
            
            else:
            
                # No barcodes mapped to this allele at this timepoint
                # No correction necessary.
            
                barcode_weight_str = ""
            
            new_barcode_file.write("%s\t%s\n" % (allele, barcode_weight_str))
    