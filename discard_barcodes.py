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

blacklisted_barcode_strs = set(['AAAAAAAAAAAAAAAA', 'CCCCCCCCCCCCCCCC', 'GGGGGGGGGGGGGGGG', 'TTTTTTTTTTTTTTTT'])

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
    sys.stderr.write("Correcting all barcodes file...\n")
    barcode_filename = "%s%s/output/all_barcodes.gz" % (config.barcode_directory, sample_name)
    new_barcode_filename = "%s%s/output/all_corrected_barcodes.gz" % (config.barcode_directory, sample_name)
        
    # Open uncorrected barcode file
    barcode_file = gzip.GzipFile(barcode_filename,"r")
    # Read header
    line = barcode_file.readline() # skip header
      
    # Create corrected barcode file  
    new_barcode_file = gzip.GzipFile(new_barcode_filename,"w")
    # Write header
    new_barcode_file.write(line)
    
    # Now merge barcodes
    
    # dictionary of corrected_barcode_id -> [count, barcode_str (or "" if not set yet)] 
    corrected_barcode_stats = {}
    # For each uncorrected barcode 
    for line in barcode_file:
        items = line.split()
        barcode_id = long(items[0])
        barcode_str = items[1]
        barcode_count = long(items[2])
        
        if barcode_id in barcode_error_map:
            # This barcode error corrects to something else.
            # Discard!
            continue
            
        if barcode_str in blacklisted_barcode_strs:
            # This barcode has a very weird form (e.g. all A's)
            # Discard!
            continue
        
        if 'N' in barcode_str:
            # A base in this barcode was hard to resolve
            # Discard!
            continue
            
        # If you made it here, you are pretty good!    
        
        # This barcode error corrects to itself
        corrected_barcode_id = barcode_id
            
        if corrected_barcode_id not in corrected_barcode_stats:
            # not in barcode stats so far, so create entry
            corrected_barcode_stats[corrected_barcode_id] = [0, barcode_str]
        else:
            # already has an entry, so just set string
            corrected_barcode_stats[corrected_barcode_id][1] = barcode_str
            
        # now add in reads
        corrected_barcode_stats[corrected_barcode_id][0] += barcode_count
        
    barcode_file.close()
    
    # now print stuff out!
    for barcode_id in sorted(corrected_barcode_stats):
        
        barcode_counts = corrected_barcode_stats[barcode_id][0]
        barcode_str = corrected_barcode_stats[barcode_id][1]
        
        if barcode_str=="":
            sys.stderr.write("No barcode string for %d with %d counts. Shouldn't happen!\n" % (barcode_id, barcode_counts))
            
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
                    if original_barcode_id in barcode_error_map:
                        # Discard!
                        continue
                        #corrected_barcode_id = barcode_error_map[original_barcode_id]
                    else:
                        corrected_barcode_id = original_barcode_id
                    
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
    