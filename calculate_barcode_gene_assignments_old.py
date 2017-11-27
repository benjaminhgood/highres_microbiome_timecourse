import sys
import pylab
import numpy
import parse_midas_data
import parse_timecourse_data
import stats_utils
import barcode_utils
import gzip

import os.path
import config

pangenome_species = parse_midas_data.parse_pangenome_species()

# long gene map
# long gene = (species, gene) tuple
# id = int
longgene_id_map = {}
id_longgene_map = []

# Create ids for species
species_id_map = {pangenome_species[i]: i for i in xrange(0,len(pangenome_species))}
id_species_map = list(pangenome_species)

# Run this algorithm separately for each sample. 
# BG: should we combine across samples? Pros: more data. Cons: could have switching across species between timepoints. 
desired_samples = parse_timecourse_data.morteza_samples
#desired_samples = [parse_timecourse_data.highcoverage_start_1]
for sample_name in desired_samples:

    sys.stderr.write("Processing sample %s...\n" % sample_name)
    
    # create barcode->species map
    sys.stderr.write("Collating species barcodes...\n")
    # first create intermediate data structure:
    # barcode_id->species_id-> # times observed
    barcode_species_weight_map = {} 
    for species_name in pangenome_species:
        
        # Don't use new species yet!
        if species_name=='new_species':
            continue
        
        #sys.stderr.write(species_name+"\n")
            
        species_id = species_id_map[species_name]
        
        # Make sure barcodes exist for this timepoint.
        # BG: aside from bugs, shouldn't this always be true? 
        if not barcode_utils.barcodes_exist(species_name, sample_name):
            continue
         
        # Load barcodes      
        allele_barcode_map = barcode_utils.parse_allele_barcodes(species_name, sample_name)

        for allele in allele_barcode_map.keys():
            for barcode_id in allele_barcode_map[allele]:
                if barcode_id not in barcode_species_weight_map:
                    barcode_species_weight_map[barcode_id] = {}
                if species_id not in barcode_species_weight_map[barcode_id]: 
                    barcode_species_weight_map[barcode_id][species_id] = 0.0
                barcode_species_weight_map[barcode_id][species_id] += 1

    sys.stderr.write("Done!\n")
    sys.stderr.write("Pruning barcode->species map...\n")
    # now calculate barcode_species_map using "unique barcodes"
    barcode_species_map = {}
    non_unique_barcodes = 0
    unique_barcodes = 0
    for barcode_id in barcode_species_weight_map.keys():
        if len(barcode_species_weight_map[barcode_id]) > 1.5:
            non_unique_barcodes+=1
        else:
            unique_barcodes+=1
            barcode_species_map[barcode_id] = barcode_species_weight_map[barcode_id].values()[0]
    sys.stderr.write("Done!\n")
    sys.stderr.write("%d unique barcodes\n" % unique_barcodes)
    sys.stderr.write("%d non-unique barcodes\n" % non_unique_barcodes)
    
    if len(barcode_species_map)==0:
        continue
    
    sys.stderr.write("Looping through unassigned genes...\n")
           
    # Now go through genes in new_species (from de novo assembly)
    # Re-written to be in online form. (ugly but more memory efficient)
    barcode_filename = "%s%s/output/new_species.barcodes.gz" % (config.barcode_directory, sample_name)
    
    barcode_file = gzip.GzipFile(barcode_filename,"r")
    barcode_file.readline() # skip header
    
    all_genes = 0
    unique_genes = 0
    non_unique_genes = 0
    for line in barcode_file:
        line = line.strip()
        items = line.split("\t")
        gene_name = items[0].strip()
        
        all_genes += 1
        
        # Make sure there are barcodes to look at
        if len(items)==1:
            continue
            
        # Get barcode ids    
        barcode_items = items[1].split(",")
        barcode_ids = []
        for barcode_item in barcode_items:
            barcode_subitems = barcode_item.split(":")
            barcode_id = long(barcode_subitems[0])
            barcode_ids.append(barcode_id)
                
        species_weight_map = {}
        for barcode_id in barcode_ids:
            if barcode_id in barcode_species_map:
                species_id = barcode_species_map[barcode_id]
                if species_id not in species_weight_map:
                    species_weight_map[species_id] = 0
                species_weight_map[species_id] += 1
                
        # no barcodes belonging to species
        if len(species_weight_map)==0:
            continue
            
        if len(species_weight_map)==1:
            unique_genes += 1
        else:
            non_unique_genes += 1
            #print species_weight_map.values()
            
    sys.stderr.write("Done!\n")
    sys.stderr.write("%d total genes\n" % all_genes)
    sys.stderr.write("%d unique genes\n" % unique_genes)
    sys.stderr.write("%d nonunique genes\n" % non_unique_genes)

