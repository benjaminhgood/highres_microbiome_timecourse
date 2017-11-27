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

# Run this algorithm separately for each sample. 
# BG: should we combine across samples? Pros: more data. Cons: could have switching across species between timepoints. 
desired_samples = parse_timecourse_data.morteza_samples
#desired_samples = [parse_timecourse_data.highcoverage_start_1]
for sample_name in desired_samples:

    sys.stderr.write("Processing sample %s...\n" % sample_name)
    
    # create barcode->species map
    sys.stderr.write("Collating species barcodes...\n")
    # first create intermediate data structure:
    # barcode_id->longgene->count
    barcode_longgene_weight_map = {} 
    for species_name in pangenome_species:
        
        # Don't use new species yet!
        if species_name=='new_species':
            continue
        
        # Make sure barcodes exist for this timepoint.
        # BG: aside from bugs, shouldn't this always be true? 
        if not barcode_utils.barcodes_exist(species_name, sample_name):
            continue
         
        # Load barcodes      
        allele_barcode_map = barcode_utils.parse_allele_barcode_tuples(species_name, sample_name)

        for allele in allele_barcode_map.keys():
        
            if allele.endswith('|A') or allele.endswith('|R'):
                # a SNP allele, don't include
                continue
        
            if len(allele_barcode_map)==0:
                continue
        
            longgene = species_name
            if longgene not in longgene_id_map:
                longgene_id_map[longgene] = len(id_longgene_map)
                id_longgene_map.append(longgene)
            
            longgene_id = longgene_id_map[longgene]
        
            for barcode_id, barcode_weight in allele_barcode_map[allele]:
                if barcode_id not in barcode_longgene_weight_map:
                    barcode_longgene_weight_map[barcode_id] = {}
                if longgene_id not in barcode_longgene_weight_map[barcode_id]: 
                    barcode_longgene_weight_map[barcode_id][longgene_id] = 0.0
                barcode_longgene_weight_map[barcode_id][longgene_id] += barcode_weight

    sys.stderr.write("Done!\n")
    
    if len(barcode_longgene_weight_map)==0:
        continue
    
    # Plot distributions of things
    total_weights = []
    max_weights = []
    for barcode_id in barcode_longgene_weight_map:
        weights = numpy.array(barcode_longgene_weight_map[barcode_id].values())
        
        total_weights.append( weights.sum() )
        max_weights.append( weights.max() )
        
    total_weights = numpy.array( total_weights )
    max_weights = numpy.array( max_weights )
    
    max_fractions = max_weights*1.0/total_weights
    
    print len(max_fractions), total_weights.sum(), numpy.median(total_weights), numpy.median(max_fractions)
    
    weight_bins = numpy.logspace(0,4,100)
    weight_bins[0] = -1
    weights = weight_bins[1:] 
    
    #ns, dummy_bins = numpy.histogram(total_weights, bins=weight_bins)
    #cdf = numpy.cumsum(ns)
    #survivals = cdf[-1]-cdf
    #pylab.step(weights, survivals,'-')
    
       