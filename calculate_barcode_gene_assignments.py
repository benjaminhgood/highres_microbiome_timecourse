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

pangenome_species = parse_midas_data.parse_pangenome_species()

# long gene map
# long gene = (species, gene) tuple
# id = int
longgene_id_map = {}
id_longgene_map = []

total_candidate_genes = set()

# Run this algorithm separately for each sample. 
# BG: should we combine across samples? Pros: more data. Cons: could have switching across species between timepoints. 
desired_samples = parse_timecourse_data.morteza_samples
#desired_samples = [parse_timecourse_data.highcoverage_start_1]
for sample_name in desired_samples:

    sys.stderr.write("Processing sample %s...\n" % sample_name)
    
    # create barcode->species map
    sys.stderr.write("Collating species barcodes...\n")
    # first create intermediate data structure:
    # barcode_id->set of longgene ids
    barcode_longgene_ids_map = {} 
    for species_name in pangenome_species:
        
        # Don't use new species yet!
        if species_name=='new_species':
            continue
        
        # Make sure barcodes exist for this timepoint.
        # BG: aside from bugs, shouldn't this always be true? 
        if not barcode_utils.barcodes_exist(species_name, sample_name):
            continue
         
        # Load barcodes      
        allele_barcode_map = barcode_utils.parse_allele_barcodes(species_name, sample_name)

        for allele in allele_barcode_map.keys():
        
            if allele.endswith('|A') or allele.endswith('|R'):
                # a SNP allele, don't include
                continue
        
            if len(allele_barcode_map)==0:
                continue
            
            # form longgene
            longgene = (species_name, allele)
            # update longgene id
            if longgene not in longgene_id_map:
                longgene_id_map[longgene] = len(id_longgene_map)
                id_longgene_map.append(longgene)
            
            longgene_id = longgene_id_map[longgene]
    
            for barcode_id in allele_barcode_map[allele]:
                if barcode_id not in barcode_longgene_ids_map:
                    barcode_longgene_ids_map[barcode_id] = set()
                barcode_longgene_ids_map[barcode_id].add(longgene_id)
                
    sys.stderr.write("Done!\n")
    
    if len(barcode_longgene_ids_map)==0:
        continue
    
    sys.stderr.write("Looping through unassigned genes...\n")
           
    # Now go through genes in new_species (from de novo assembly)
    # Re-written to be in online form. (ugly but more memory efficient)
    barcode_filename = "%s%s/output/new_species.barcodes.gz" % (config.barcode_directory, sample_name)
    
    barcode_file = gzip.GzipFile(barcode_filename,"r")
    barcode_file.readline() # skip header
    
    all_genes = 0
    barcode_genes = 0
    unique_genes = 0
    non_unique_genes = 0
    candidate_genes = set()
    for line in barcode_file:
        line = line.strip()
        items = line.split("\t")
        gene_name = items[0].strip()
        
        all_genes += 1
        
        if (all_genes % 100000) == 0:
            print all_genes, barcode_genes, len(candidate_genes)
        
        # Make sure there are barcodes to look at
        if len(items)==1:
            continue
            
        # Get barcode ids    
        barcode_items = items[1].split(",")
        
        # Can't have 5 good barcodes if less than 5 barcodes
        if len(barcode_items) < 5:
            continue
        
        barcode_ids = []
        for barcode_item in barcode_items:
            barcode_subitems = barcode_item.split(":")
            barcode_id = long(barcode_subitems[0])
            barcode_ids.append(barcode_id)
        
                
        num_barcodes = 0
        
        for barcode_id in barcode_ids:
            if barcode_id in barcode_longgene_ids_map:
                num_barcodes += 1
                    
        # only look at genes with a decent number of barcodes            
        if num_barcodes<5:
            continue
        
        barcode_genes += 1 
        
        longgene_id_set = set()
        num_duplicates = 0
        for barcode_id in barcode_ids:
            if barcode_id in barcode_longgene_ids_map:
                if not longgene_id_set.isdisjoint( barcode_longgene_ids_map[barcode_id] ):
                    num_duplicates += 1
                
                longgene_id_set.update( barcode_longgene_ids_map[barcode_id] )
        
        if num_duplicates >= 0.5*num_barcodes:    
            candidate_genes.add(gene_name)
    
    barcode_file.close()
    
    print len(candidate_genes), "pre-candidate genes"
    
    #total_candidate_genes.update(candidate_genes)
    #continue 
     
    new_candidate_genes = set()
    
    barcode_file = gzip.GzipFile(barcode_filename,"r")
    barcode_file.readline() # skip header
    for line in barcode_file:
        line = line.strip()
        items = line.split("\t")
        gene_name = items[0].strip()
        
        if gene_name not in candidate_genes:
            continue
        
        # Get barcode ids    
        barcode_items = items[1].split(",")
        
        barcode_ids = []
        for barcode_item in barcode_items:
            barcode_subitems = barcode_item.split(":")
            barcode_id = long(barcode_subitems[0])
            barcode_ids.append(barcode_id)
                
        num_barcodes = 0
        longgene_id_counter = collections.Counter()
        for barcode_id in barcode_ids:
            if barcode_id in barcode_longgene_ids_map:
                num_barcodes += 1
                longgene_id_counter.update( barcode_longgene_ids_map[barcode_id] )
        
        # calculate max longgene weight:
        max_longgene_id, max_longgene_weight = longgene_id_counter.most_common(1)[0]
        
        if max_longgene_weight < 2:
            continue
            
        if max_longgene_weight < 0.5*num_barcodes:
            continue
            
        new_candidate_genes.add(gene_name)
        
        
        # make a deeper species->gene->weight_map
        #species_gene_weight_map = {}
        #for longgene,weight in longgene_weight_map.iteritems():
        #    species_name, gene_name = longgene
        #    if species_name not in species_gene_weight_map:
        #        species_gene_weight_map[species_name] = {}
        #    if gene_name not in species_gene_weight_map[species_name]:
        #        species_gene_weight_map[species_name][gene_name] = 0
            
        #    species_gene_weight_map[species_name][gene_name]+=weight 
            
        #if len(species_weight_map)==1:
        #    unique_genes += 1
        #else:
        #    non_unique_genes += 1
        #    print num_barcodes, max_species_weight, max_longgene_weight
        #    print species_gene_weight_map
            #print species_weight_map.values()
    
    total_candidate_genes.update(new_candidate_genes)
            
    sys.stderr.write("Done!\n")
    sys.stderr.write("%d total genes\n" % all_genes)
    sys.stderr.write("%d barcode covered genes\n" % barcode_genes)
    sys.stderr.write("%d candidate genes\n" % len(total_candidate_genes))
    
    
sys.stderr.write("Done!\n")
sys.stderr.write("%d candidate genes\n" % len(total_candidate_genes))
