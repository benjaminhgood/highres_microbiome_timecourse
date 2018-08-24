###########################################################################
#
# The purpose of this script is to figure out which genes (and therefore which species)
# are linked to a given allele
#
###########################################################################


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
from scipy.special import gammaln as loggamma
from scipy.special import gammainc
from math import log,exp
from scipy.stats import poisson

from numpy.random import shuffle

pangenome_species = parse_midas_data.parse_pangenome_species()

    
###########################################################################
#
# Standard header to read in argument information
#
###########################################################################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("features_filename", help="list of alleles to check")
args = parser.parse_args()


features_filename = args.features_filename
corrected = True
core_genes_only = True
external_core_genes=False
disable_filter=True
remove_focal_from_targets=True
Pstar = 0.1

# create output filename
features_filename_items = features_filename.split("/")[-1].split(".")
output_filename = config.barcode_directory+"barcode_trajectories/"+features_filename_items[0]
output_filename += ".p"

# load candidate genes from supplied file
feature_allele_map = {} # mapping from feature id to set of allele name
feature_species_map = {} # mapping from feature id to focal species
allele_features_map = {} # mapping from allele name to feature id
feature_idx = -1
allele_names = []
features = []
features_file = open(features_filename,"r")
# each line contains a feature. starts with species name, then TSL of alleles
for line in features_file:
    feature_idx+=1
    items = line.split()
    species_name = items[0]
    if feature_idx not in feature_allele_map:
        features.append(feature_idx)
        feature_allele_map[feature_idx] = set()
        feature_species_map[feature_idx] = species_name
    
    for item in items[1:]:
        allele_name_prefix = item.strip()
        if '|' in allele_name_prefix and not (allele_name_prefix[-1]=='A' or allele_name_prefix[-1]=='R'):
            # it's a SNP.. add both A and R
            sub_allele_names = [allele_name_prefix+'|A', allele_name_prefix+'|R']
        else:
            sub_allele_names = [allele_name_prefix]
            
        for allele_name in sub_allele_names:
            allele_features_map[allele_name] = feature_idx
            feature_allele_map[feature_idx].add(allele_name)
            allele_names.append(allele_name)
features_file.close()


# species map
# id = int
species_id_map = {}
id_species_map = []

# long gene map
# long gene = (species, gene) tuple
# id = int
longgene_id_map = {}
id_longgene_map = []

# core longgenes
# Load core genes across those species
import core_gene_utils
core_longgenes = set()
for species_name in pangenome_species:
    core_genes = core_gene_utils.parse_core_genes(species_name, external_filtering=external_core_genes)
    for gene_name in core_genes:
        core_longgenes.add((species_name, gene_name))

#desired_samples = ['1021']
desired_samples = parse_timecourse_data.highcoverage_samples
#desired_samples = parse_timecourse_data.morteza_samples
#desired_samples = [parse_timecourse_data.highcoverage_start_2, parse_timecourse_data.highcoverage_antibiotic, parse_timecourse_data.highcoverage_postantibiotic, parse_timecourse_data.highcoverage_end]
#desired_samples = desired_samples[0:2]
#desired_samples = parse_timecourse_data.alistipes_onderdonkii_gene_gain_samples
#desired_samples = [parse_timecourse_data.highcoverage_antibiotic]


def calculate_M_from_D(D):
    
    if D<10:
        return D
    else:
        return 10


# Memory requirements are big, so we have to break things up into smaller batches
# of unassigned genes. The output is then grouped into a single file.
batch_size = 5000

num_batches = long(ceil(len(allele_names)*1.0/batch_size))

if num_batches>1.5 and (len(allele_names)%batch_size)*1.0/batch_size < 0.1:
    # just expand the batch size a little bit
    batch_size = long(ceil(len(allele_names)*1.0/(num_batches-1)))
    num_batches = long(ceil(len(allele_names)*1.0/batch_size))

sys.stderr.write("Divided %d features into %d batches of size %d.\n" % (len(features), num_batches, batch_size))

for batch in xrange(0,num_batches):
    
    sys.stderr.write("Processing batch %d...\n" % batch)
    
    feature_barcode_timecourse = {}
    for i in xrange(0,len(features)):
        if (i/batch_size) == batch:
            feature = features[i]
            feature_barcode_timecourse[feature] = {'all':{}, 'longgenes':{}}

    sys.stderr.write("(%d features)\n" % len(feature_barcode_timecourse))
    
    for sample_idx in xrange(0,len(desired_samples)):
        # Process each sample separately
        # (combine results at end)
    
        sample_name = desired_samples[sample_idx]

        num_significant_hits = 0 
        
        sys.stderr.write("Processing sample %s...\n" % sample_name)
        sys.stderr.write("Loading depth map...\n")
        barcode_depth_map = barcode_utils.parse_barcode_depth_map(sample_name,corrected=corrected,min_depth=2.5)
        sys.stderr.write("Done!\n")
        
        if len(barcode_depth_map)==0:
            continue
        
        sys.stderr.write("Postprocessing depth map...\n")
        
        total_barcodes = len(barcode_depth_map)
        avg_M = 0
        for D in barcode_depth_map.itervalues():
            avg_M += calculate_M_from_D(D)
        avg_M /= total_barcodes
        sys.stderr.write("Done!\n")
         
        # create barcode->species map
        sys.stderr.write("Collating species barcodes...\n")
        # first create intermediate data structure:
        # barcode_id->set of longgene ids
        barcode_longgene_ids_map = {} 
        
        focal_feature_barcode_map = {} # Map from feature to barcode ids that that feature has
        
        target_longgene_barcode_fraction_map = {} # map from target longgene to fraction of barcodes in set (i.e., Mavg*pj from notes)
        
        for species_name in pangenome_species:
        
            # 'new_species' is set of de-novo assembled genes
            # do not use these for now
            if species_name=='new_species':
                continue
        
            # Make sure barcodes exist for this species at this timepoint.
            if not barcode_utils.barcodes_exist(species_name, sample_name):
                continue
         
            # Load barcodes      
            allele_barcode_map = barcode_utils.parse_allele_barcode_tuples(species_name, sample_name,corrected=corrected)

            for allele in allele_barcode_map:
        
                if allele in allele_features_map:
                    # One of the alleles we are trying to track...
                    
                    # get feature ID
                    feature = allele_features_map[allele] 
                    
                    # Add in the barcodes
                    for barcode_id, barcode_weight in allele_barcode_map[allele]:
                
                        # don't include barcode if too little coverage
                        # Poor man's error correction
                        if barcode_id not in barcode_depth_map: 
                            continue
                        
                        # Start list of barcodes for that feature
                        if feature not in focal_feature_barcode_map:
                            focal_feature_barcode_map[feature] = set()
                        # Add barcode to list
                        focal_feature_barcode_map[feature].add(barcode_id)
                    
                    if remove_focal_from_targets:
                    # then take allele out of circulation
                        continue           
        
                # Only look at barcodes that map to genes
                if allele.endswith('|A') or allele.endswith('|R'):
                    continue
        
                # Only look at genes that have barcodes
                if len(allele_barcode_map[allele])<0.5:
                    continue
            
                # Get gene ID
                longgene = (species_name, allele)
                
                # If specified, only look at genes in the core genome
                if (core_genes_only) and (longgene not in core_longgenes):
                    continue
                    
                # update longgene id database
                if longgene not in longgene_id_map:
                    longgene_id_map[longgene] = len(id_longgene_map)
                    id_longgene_map.append(longgene)
                longgene_id = longgene_id_map[longgene]
    
                # Get species ID
                if species_name not in species_id_map:
                    species_id_map[species_name] = len(id_species_map)
                    id_species_map.append(species_name)
                species_id = species_id_map[species_name]
            
                # estimated # of fragments per barcode
                total_B = 0

                for barcode_id, barcode_weight in allele_barcode_map[allele]:
                
                    # don't include barcode if too little coverage
                    # Poor man's error correction
                    if barcode_id not in barcode_depth_map: 
                        continue
                    
                    total_B += 1
                        
                    if barcode_id not in barcode_longgene_ids_map:
                        barcode_longgene_ids_map[barcode_id] = set()
                        #barcode_depth_map[barcode_id] = 0
                    
                    barcode_longgene_ids_map[barcode_id].add(longgene_id)
                    #barcode_depth_map[barcode_id] += barcode_weight
                
                target_longgene_barcode_fraction_map[longgene_id] = total_B*1.0/total_barcodes
                
        total_features = len(features)
        total_longgenes = len(target_longgene_barcode_fraction_map)
        log_corrected_Pstar = log(Pstar/total_longgenes/total_features) # correct for multiple hypothesis testing
                
        sys.stderr.write("Done! Loaded %d barcodes\n" % len(barcode_longgene_ids_map))
        sys.stderr.write("(other count: %d)\n" % len(barcode_depth_map))
        sys.stderr.write("Checking linkage with %d target genes\n" % total_longgenes)
        
        
        # Make sure we loaded some barcodes...
        if len(barcode_longgene_ids_map)==0:
            continue
      
        sys.stderr.write("Looping through target features...\n")
        
        for feature in focal_feature_barcode_map:
        
            # B = total # barcodes for this feature
            B = len(focal_feature_barcode_map[feature])
        
            # Get total M
            Mtot = 0
            for barcode_id in focal_feature_barcode_map[feature]:
                Mtot += calculate_M_from_D(barcode_depth_map[barcode_id])
            
            
            # Get barcode ids that map to this gene    
            barcode_ids = []
            for barcode_id in focal_feature_barcode_map[feature]:
                
                # If it doesn't map to anywhere we know of, 
                # skip it (saves some time later)
                if barcode_id not in barcode_longgene_ids_map:
                    continue
                     
                barcode_ids.append(barcode_id)
            
            num_barcodes = len(barcode_ids)
            
            # only create entries if there are at least two barcodes
            if num_barcodes<0.5:
                continue
            
            #sys.stderr.write("Processing %d barcodes for %s...\n" % (num_barcodes,allele))
               
            # Counts number of barcodes that map to different genes
            # within species
            longgene_id_counter = collections.Counter()
            for barcode_id in barcode_ids:
                longgene_id_counter.update( barcode_longgene_ids_map[barcode_id] )
                    
            # create entry for feature        
            feature_barcode_timecourse[feature]['all'][sample_idx] = B
        
            # Now add in stuff from gene counter
            for longgene_id in longgene_id_counter:
            
                # Observed of shared barcodes
                S = longgene_id_counter[longgene_id]
            
                if S>=1: 
                
                    # Calculate how surprised we are...
                
                    # Calculate expected # of shared barcodes
                    S0 = Mtot/avg_M*target_longgene_barcode_fraction_map[longgene_id]    
            
                    if S==1:
                        logP = log(1-exp(-S0))
                    elif S==2:
                        logP = log(1-(1+S0)*exp(-S0))
                    elif (S-1)>100*S0 or (S>100 and S>10*S0):
                        # Use asymptotic approximation
                        # for large argument
                        logP = S*log(S0)-S0+loggamma(S+1)
                    else:
                        # Full calculation
                        logP = log(gammainc(S, S0))
                    
                    if logP < log_corrected_Pstar:
                        
                        num_significant_hits += 1
                        
                        # Surprising! Save it! 
                
                        if longgene_id not in feature_barcode_timecourse[feature]['longgenes']:
                            feature_barcode_timecourse[feature]['longgenes'][longgene_id] = {}
            
                        feature_barcode_timecourse[feature]['longgenes'][longgene_id][sample_idx] = (S, logP)
            
            # Done!
            
        sys.stderr.write("Done! %d significant hits\n" % num_significant_hits)
    
    # Now fill in missing zeros:
    # Look at all timepoints in order
    # If a gene doesn't have a timepoint in 'all', add it in. 
    # If a sub-gene doesn't have a timepoint, then add it in (0)

    # unassigned gene name -> 'all' | 'species_name' -> 'gene_name' -> vector of counts
    processed_feature_barcode_timecourse = {}
    for feature in feature_barcode_timecourse:

        # first do "all" barcodes
        all_barcodes = []
        for sample_idx in xrange(0,len(desired_samples)):
            if sample_idx in feature_barcode_timecourse[feature]['all']:
                all_barcodes.append( feature_barcode_timecourse[feature]['all'][sample_idx] )
            else:
                all_barcodes.append( 0 )
    
        all_barcodes = numpy.array(all_barcodes)
 
        processed_feature_barcode_timecourse[feature] = {'all': all_barcodes, 'species': {}}
        processed_feature_barcode_timecourse[feature]['focal_species'] = feature_species_map[feature]
        processed_feature_barcode_timecourse[feature]['focal_alleles'] = feature_allele_map[feature]
    
        # Then do barcodes per gene 
        for longgene_id in feature_barcode_timecourse[feature]['longgenes']:
        
            barcodes = []
            pvalues = []
            for sample_idx in xrange(0,len(desired_samples)):
        
                if sample_idx in feature_barcode_timecourse[feature]['longgenes'][longgene_id]:
                    barcodes.append( feature_barcode_timecourse[feature]['longgenes'][longgene_id][sample_idx][0] )
                    pvalues.append( feature_barcode_timecourse[feature]['longgenes'][longgene_id][sample_idx][1] )
                    
                else:
                    barcodes.append( 0 )
                    pvalues.append( 0 )
                    
            barcodes = numpy.array(barcodes)
            pvalues = numpy.array(pvalues)
        
            other_species_name, other_gene_name = id_longgene_map[longgene_id]
        
            if other_species_name not in processed_feature_barcode_timecourse[feature]['species']:
                processed_feature_barcode_timecourse[feature]['species'][other_species_name] = {}
            processed_feature_barcode_timecourse[feature]['species'][other_species_name][other_gene_name] = (barcodes, pvalues)

pickle.dump( processed_feature_barcode_timecourse, open( output_filename, "wb" ) )
sys.stderr.write("Done!\n")

