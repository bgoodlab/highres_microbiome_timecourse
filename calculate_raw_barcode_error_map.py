#!/usr/bin/python

#########################################################
#   Barcode error correction
#   to run: python calculate_barcode_error_map.py
#
#   Stephen Martis
#   BHG: edited to work with native MIDAS output paths / files
#        and made some speed / memory usage optimizations
#########################################################

import gzip
import sys
import csv
from datetime import datetime
import Levenshtein

input_file = sys.stdin
output_file = sys.stdout



#Barcode edit distribution/error correction
startTime = datetime.now()

#input_fastq_filename = "%s%s/output/all_barcodes.gz" % (config.barcode_directory, sample_name)
#output_barcode_map_filename = "%s%s/output/barcode_map.gz" % (config.barcode_directory, sample_name)

    
# Populate these two dictionaries with barcode data    
barcode_id_map = {} # map from barcode str to barcode id
barcode_str_map = {} # map from barcode id to barcode str
barcode_depth_map = {} # map from barcode id to depth (total number of read counts)
current_id = 0
    
sys.stderr.write("Loading barcodes from fastq...\n")
line = input_file.readline().strip()
while line!="":
        
    # Find barcode str
    read_name = line.split()[0]
    barcode_str = read_name.split(":")[-1]
    barcode_str = barcode_str.split(",")[0]
        
    if barcode_str not in barcode_id_map:
        barcode_id = current_id
        current_id += 1
        barcode_str_map[barcode_id] = barcode_str
        barcode_depth_map[barcode_id] = 0
    else:
        barcode_id = barcode_id_map[barcode_str]
        
    barcode_depth_map[barcode_id] += 1
        
    # Advance record
    input_file.readline() # DNA sequence
    input_file.readline() # +
    input_file.readline() # quality score
    line = input_file.readline() # next read header
        
        
sys.stderr.write("Done! %s\n" % str(datetime.now()-startTime))
    
# Now everything can proceed like last time
sys.stderr.write("Sorting barcodes in descending order of counts...\n")
sorted_barcode_ids = sorted(barcode_depth_map.iterkeys(), key=lambda id: barcode_depth_map[id], reverse=True)
sys.stderr.write("Done! %s\n" % str(datetime.now()-startTime))
    
    
output_file.write("\t".join(['barcode_str','corrected_barcode_str', 'num_reads']))
output_file.write("\n")
    
sys.stderr.write("Building error map...\n")
inverted_deletion_map = {}
num_processed = 0
num_errors = 0
for barcode_id in sorted_barcode_ids:
    num_processed += 1
        
    barcode_depth = barcode_depth_map[barcode_id]
    barcode_str = barcode_str_map[barcode_id]
            
    if num_processed%100000==0:
        sys.stderr.write("%d errors in %d barcodes\n" % (num_errors, num_processed))
    barcode_str = barcode_str_map[barcode_id]
        
    deletion_ids = []
        
    match = False
    matched_barcode_id = -1
    for i in range(len(barcode_str)):
        
        new_barcode_str = barcode_str[:i] + barcode_str[(i+1):]
            
        edit_barcode_id = hash(new_barcode_str)
            
        if edit_barcode_id in inverted_deletion_map:
            # potential match
            # check to see if it has right levenstein distance
            for other_barcode_id in inverted_deletion_map[edit_barcode_id]:
                if Levenshtein.distance(barcode_str_map[barcode_id], barcode_str_map[other_barcode_id])==1:
                    match = True
                    matched_barcode_id = other_barcode_id
                    break
            
        if match:
            break
        else:
            deletion_ids.append(edit_barcode_id)
        
    if match:
        num_errors += 1
        # This barcode corrects to something else!
            
        corrected_barcode_str = barcode_str_map[matched_barcode_id]
        
            
    else:
        # This barcode is a new template.
        # Add it to database
        for edit_barcode_id in deletion_ids:
            if edit_barcode_id not in inverted_deletion_map:
                inverted_deletion_map[edit_barcode_id] = []
            inverted_deletion_map[edit_barcode_id].append(barcode_id)
                
        corrected_barcode_str = barcode_str
            
    output_file.write("%s\t%s\t%d\n" % (barcode_str, corrected_barcode_str, barcode_depth))
        
output_file.close()
            
sys.stderr.write("Done! %s\n" % str(datetime.now()-startTime))
sys.stderr.write("%d errors in %d barcodes\n" % (num_errors, num_processed))
    
