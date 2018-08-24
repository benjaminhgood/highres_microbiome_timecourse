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

# MIDAS specific imports
import config
import parse_midas_data
import parse_timecourse_data

for sample_name in parse_timecourse_data.morteza_samples:
    sys.stderr.write("Processing %s...\n" % sample_name)

    #Barcode edit distribution/error correction
    startTime = datetime.now()

    input_barcode_filename = "%s%s/output/all_barcodes.gz" % (config.barcode_directory, sample_name)
    output_barcode_map_filename = "%s%s/output/barcode_map.gz" % (config.barcode_directory, sample_name)

    
    # Populate these two dictionaries with barcode data    
    barcode_id_map = {} # map from barcode str to barcode id
    barcode_str_map = {} # map from barcode id to barcode str
    barcode_depth_map = {} # map from barcode id to depth (total number of read counts)
     
    
    sys.stderr.write("Loading barcodes...\n")
    file = gzip.GzipFile(input_barcode_filename,"r")
    file.readline() # header
    for line in file:
        
        items = line.split("\t")
    
        barcode_id = long(items[0])
        barcode_str = items[1]
        depth = long(items[2])
        
        barcode_id_map[barcode_str] = barcode_id
        barcode_str_map[barcode_id] = barcode_str
        barcode_depth_map[barcode_id] = depth
        
    file.close()
    sys.stderr.write("Done! %s\n" % str(datetime.now()-startTime))
    
    sys.stderr.write("Sorting barcodes in descending order of counts...\n")
    sorted_barcode_ids = sorted(barcode_depth_map.iterkeys(), key=lambda id: barcode_depth_map[id], reverse=True)
    sys.stderr.write("Done! %s\n" % str(datetime.now()-startTime))
    
    
    output_file = gzip.GzipFile(output_barcode_map_filename,"w")    
    output_file.write("\t".join(['barcode_id','error_corrected_barcode_id']))
    output_file.write("\n")
    
    inverted_deletion_map = {}
    sys.stderr.write("Building error map...\n")
    num_processed = 0
    num_errors = 0
    for barcode_id in sorted_barcode_ids:
        num_processed += 1
        
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
            output_file.write("%d\t%d\n" % (barcode_id, matched_barcode_id))
        else:
            # This barcode is a new template.
            # Add it to database
            for edit_barcode_id in deletion_ids:
                if edit_barcode_id not in inverted_deletion_map:
                    inverted_deletion_map[edit_barcode_id] = []
                inverted_deletion_map[edit_barcode_id].append(barcode_id)
                    
    
    sys.stderr.write("Done! %s\n" % str(datetime.now()-startTime))
    sys.stderr.write("%d errors in %d barcodes\n" % (num_errors, num_processed))
    output_file.close()
    
