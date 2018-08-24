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
import cPickle as pickle
from math import ceil
import sys
import argparse

from numpy.random import shuffle

pangenome_species = parse_midas_data.parse_pangenome_species()

###########################################################################
#
# Standard header to read in argument information
#
###########################################################################import argparse
parser = argparse.ArgumentParser()
parser.add_argument("gene_filename", help="list of de novo genes to check")
parser.add_argument("--bootstrapped", help="Perform search on bootstrapped barcodes", action="store_true")
args = parser.parse_args()


gene_filename = args.gene_filename
bootstrapped = args.bootstrapped
corrected = True

# create output filename
gene_filename_items = gene_filename.split(".")
output_filename = gene_filename_items[0]
if bootstrapped:
    output_filename += ".bootstrapped"
output_filename += ".p"

# load candidate genes from supplied file
gene_names = []
gene_file = open(gene_filename,"r")
for line in gene_file:
    gene_name = line.strip()
    gene_names.append(gene_name)
gene_file.close()

# species map
# id = int
species_id_map = {}
id_species_map = []

# long gene map
# long gene = (species, gene) tuple
# id = int
longgene_id_map = {}
id_longgene_map = []

desired_samples = parse_timecourse_data.morteza_samples
#desired_samples = [parse_timecourse_data.highcoverage_start_1]
#desired_samples = desired_samples[0:2]

# make the samples go in temporal order
sample_time_map = parse_timecourse_data.parse_sample_time_map()
species_times, species_time_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, desired_samples)
desired_samples = numpy.array(desired_samples)[species_time_idxs]

###############
#
processed_gene_barcode_timecourse = {}
#
# The final post-processed data structure.
# keys are unassigned genes. values are dictionaries with two elements:
# 
# 'all':  
# values are {'all': all_barcodes, 'species': {}}

# Memory requirements are big, so we have to break things up into smaller batches
# of unassigned genes. The output is then grouped into a single file.
batch_size = 5000

num_batches = long(ceil(len(gene_names)*1.0/batch_size))

if (len(gene_names)%batch_size)*1.0/batch_size < 0.1:
    # just expand the batch size a little bit
    batch_size = long(ceil(len(gene_names)*1.0/(num_batches-1)))
    num_batches = long(ceil(len(gene_names)*1.0/batch_size))

sys.stderr.write("Divided %d genes into %d batches of size %d.\n" % (len(gene_names), num_batches, batch_size))

for batch in xrange(0,num_batches):
    
    sys.stderr.write("Processing batch %d...\n" % batch)
    
    gene_barcode_timecourse = {}
    for i in xrange(0,len(gene_names)):
        if (i/batch_size) == batch:
            gene_name = gene_names[i]
            gene_barcode_timecourse[gene_name] = {'all':{}, 'longgenes':{}, 'species':{}}

    sys.stderr.write("(%d genes)\n" % len(gene_barcode_timecourse))
    
    for sample_idx in xrange(0,len(desired_samples)):
        sample_name = desired_samples[sample_idx]

        sys.stderr.write("Processing sample %s...\n" % sample_name)
         
        # create barcode->species map
        sys.stderr.write("Collating species barcodes...\n")
        # first create intermediate data structure:
        # barcode_id->set of longgene ids
        barcode_longgene_ids_map = {} 
        barcode_depth_map = {}
        for species_name in pangenome_species:
        
            # Don't use new species yet!
            if species_name=='new_species':
                continue
        
            # Make sure barcodes exist for this timepoint.
            if not barcode_utils.barcodes_exist(species_name, sample_name):
                continue
         
            # Load barcodes      
            allele_barcode_map = barcode_utils.parse_allele_barcode_tuples(species_name, sample_name,corrected=corrected, bootstrapped=bootstrapped)

            for allele in allele_barcode_map.keys():
        
                # Only look at barcodes that map to genes
                if allele.endswith('|A') or allele.endswith('|R'):
                    continue
        
                # Only look at genes that have barcodes
                if len(allele_barcode_map)==0:
                    continue
            
                # Get gene ID
                longgene = (species_name, allele)
                # update longgene id
                if longgene not in longgene_id_map:
                    longgene_id_map[longgene] = len(id_longgene_map)
                    id_longgene_map.append(longgene)
                longgene_id = longgene_id_map[longgene]
    
                # Get species ID
                if species_name not in species_id_map:
                    species_id_map[species_name] = len(id_species_map)
                    id_species_map.append(species_name)
                species_id = species_id_map[species_name]
            
                for barcode_id, barcode_weight in allele_barcode_map[allele]:
                
                    if barcode_id not in barcode_longgene_ids_map:
                        barcode_longgene_ids_map[barcode_id] = set()
                        barcode_depth_map[barcode_id] = 0
                    
                    barcode_longgene_ids_map[barcode_id].add(longgene_id)
                    barcode_depth_map[barcode_id] += barcode_weight
                
        sys.stderr.write("Done! Loaded %d barcodes\n" % len(barcode_longgene_ids_map))
        sys.stderr.write("(other count: %d)\n" % len(barcode_depth_map))
        
        if len(barcode_longgene_ids_map)==0:
            continue
    
        sys.stderr.write("Postprocessing barcode map!\n")
        # Remake barcode maps...
        new_barcode_longgene_ids_map = {}
        new_barcode_species_ids_map = {}
        
        for barcode_id in barcode_longgene_ids_map:
            
            if barcode_depth_map[barcode_id] > 2.5: # Poor man's error correction
                
                for longgene_id in barcode_longgene_ids_map[barcode_id]:
                    
                    # get species ID
                    species_id = species_id_map[id_longgene_map[longgene_id][0]]
        
                    
                    if barcode_id not in new_barcode_longgene_ids_map:
                        new_barcode_longgene_ids_map[barcode_id] = set()
                        new_barcode_species_ids_map[barcode_id] = set() 
                        
                    new_barcode_longgene_ids_map[barcode_id].add(longgene_id)
                    new_barcode_species_ids_map[barcode_id].add(species_id)
        
        barcode_longgene_ids_map = new_barcode_longgene_ids_map
        barcode_species_ids_map = new_barcode_species_ids_map
                    
        sys.stderr.write("Done!\n")
        sys.stderr.write("Looping through unassigned genes...\n")
           
        # Now go through genes in new_species (from de novo assembly)
        # Re-written to be in online form. (ugly but more memory efficient)
        barcode_filename = "%s%s/output/new_species.corrected_barcodes.gz" % (config.barcode_directory, sample_name)
    
        barcode_file = gzip.GzipFile(barcode_filename,"r")
        barcode_file.readline() # skip header
    
        all_genes = 0
        for line in barcode_file:
            line = line.strip()
            items = line.split("\t")
            gene_name = items[0].strip()

            if gene_name not in gene_barcode_timecourse:
                continue
        
            all_genes += 1
        
            if (all_genes % 1000) == 0:
                sys.stderr.write("%d...\n" % all_genes)
            
            # Make sure there are barcodes to look at
            if len(items)==1:
                continue
            
            # Get barcode ids that map to this gene    
            barcode_ids = []
            barcode_items = items[1].split(",")
            for barcode_item in barcode_items:
                barcode_subitems = barcode_item.split(":")
                barcode_id = long(barcode_subitems[0])
                
                # If it doesn't map to anywhere we know of, 
                # skip it
                if barcode_id not in barcode_longgene_ids_map:
                    continue
                
                # Poor man's error correction
                # Only look at barcodes that map to multiple places
                if barcode_depth_map[barcode_id] < 2.5:
                    continue
                    
                
                     
                barcode_ids.append(barcode_id)
            
            num_barcodes = len(barcode_ids)
            
            # only create entries if there are at least two barcodes
            if num_barcodes<1.5:
                continue
                
            # Count number of barcodes that map to different species
            species_id_counter = collections.Counter()
            for barcode_id in barcode_ids:
                species_id_counter.update( barcode_species_ids_map[barcode_id] )
                
            # only create an entry if at least one species 
            # has two barcodes that map to it
            max_per_species = max(species_id_counter.values())
            if max_per_species<1.5:
                continue
            
            # Counts number of barcodes that map to different genes
            # within species
            longgene_id_counter = collections.Counter()
            for barcode_id in barcode_ids:
                longgene_id_counter.update( barcode_longgene_ids_map[barcode_id] )
                    
            # create entry for gene        
            gene_barcode_timecourse[gene_name]['all'][sample_idx] = num_barcodes
        
            # Now add in stuff from gene counter
            for longgene_id in longgene_id_counter:
            
                if longgene_id not in gene_barcode_timecourse[gene_name]['longgenes']:
                    gene_barcode_timecourse[gene_name]['longgenes'][longgene_id] = {}
            
                gene_barcode_timecourse[gene_name]['longgenes'][longgene_id][sample_idx] = longgene_id_counter[longgene_id]
            
            # Now add in stuff from species counter
            for species_id in species_id_counter:
            
                if species_id not in gene_barcode_timecourse[gene_name]['species']:
                    gene_barcode_timecourse[gene_name]['species'][species_id] = {}
            
                gene_barcode_timecourse[gene_name]['species'][species_id][sample_idx] = species_id_counter[species_id]    
            
        sys.stderr.write("Done!\n")
    
    # Now fill in missing zeros:
    # Look at all timepoints in order
    # If a gene doesn't have a timepoint in 'all', add it in. 
    # If a sub-gene doesn't have a timepoint, then add it in (0)

    # unassigned gene name -> 'all' | 'species_name' -> 'gene_name' -> vector of counts

    for gene_name in gene_barcode_timecourse:

        # first do "all" barcodes
        all_barcodes = []
        for sample in desired_samples:
            if sample in gene_barcode_timecourse[gene_name]['all']:
                all_barcodes.append( gene_barcode_timecourse[gene_name]['all'][sample] )
            else:
                all_barcodes.append( 0 )
    
        all_barcodes = numpy.array(all_barcodes)
    
        # Want to have at least 5 barcodes at at least one timepoint. 
        if (all_barcodes<5).all():
            continue
         
        processed_gene_barcode_timecourse[gene_name] = {'all': all_barcodes, 'species': {}}
    
        # Then do all barcodes per species
        for species_id in gene_barcode_timecourse[gene_name]['species']:
        
            barcodes = []
            for sample in desired_samples:
        
                #if longgene_id not in gene_barcode_timecourse[gene_name]['longgenes']:
                #    print gene_barcode_timecourse[gene_name]['longgenes'].keys(), longgene_id
        
                if sample in gene_barcode_timecourse[gene_name]['species'][species_id]:
                    barcodes.append( gene_barcode_timecourse[gene_name]['species'][species_id][sample] )
                else:
                    barcodes.append( 0 )
                
            barcodes = numpy.array(barcodes)
        
            # Make sure there is at least one "good" timepoint. 
            if not ((all_barcodes>=5)*(barcodes>=4)*(barcodes>0.75*all_barcodes)).any():
                continue
        
            other_species_name = id_species_map[species_id]
            processed_gene_barcode_timecourse[gene_name]['species'][other_species_name] = {'any': barcodes}
        
        # Then do barcodes per gene 
        for longgene_id in gene_barcode_timecourse[gene_name]['longgenes']:
        
            other_species_name, other_gene_name = id_longgene_map[longgene_id]
        
            if other_species_name not in processed_gene_barcode_timecourse[gene_name]['species']:
                # Don't look at genes in species that don't have enough barcodes overall...
                continue
        
            barcodes = []
            for sample in desired_samples:
        
                #if longgene_id not in gene_barcode_timecourse[gene_name]['longgenes']:
                #    print gene_barcode_timecourse[gene_name]['longgenes'].keys(), longgene_id
        
                if sample in gene_barcode_timecourse[gene_name]['longgenes'][longgene_id]:
                    barcodes.append( gene_barcode_timecourse[gene_name]['longgenes'][longgene_id][sample] )
                else:
                    barcodes.append( 0 )
                    
            barcodes = numpy.array(barcodes)
        
            # Make sure there is at least one "good" timepoint. 
            if not ((all_barcodes>=5)*(barcodes>=3)*(barcodes>0.2*all_barcodes)).any():
                continue
        
            other_species_name, other_gene_name = id_longgene_map[longgene_id]
            
            processed_gene_barcode_timecourse[gene_name]['species'][other_species_name][other_gene_name] = barcodes

        deleted_species = []
        for species_name in processed_gene_barcode_timecourse[gene_name]['species']:
            # No genes to keep track of!
            if len( processed_gene_barcode_timecourse[gene_name]['species'][species_name] ) == 1:
                deleted_species.append(species_name)
            
        for species_name in deleted_species:
            del processed_gene_barcode_timecourse[gene_name]['species'][species_name]

        # No species to keep track of!
        if len( processed_gene_barcode_timecourse[gene_name]['species'] ) == 0:
            del processed_gene_barcode_timecourse[gene_name]

pickle.dump( processed_gene_barcode_timecourse, open( output_filename, "wb" ) )
#print processed_gene_barcode_timecourse    
sys.stderr.write("Done!\n")

