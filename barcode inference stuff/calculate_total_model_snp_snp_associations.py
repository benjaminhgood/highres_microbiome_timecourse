###########################################################################
#
# The purpose of this script is to figure out which genes (and therefore which species)
# are linked to a given allele
#
###########################################################################


import sys
import matplotlib as mpl
mpl.use('Agg')
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
from scipy.special import expm1
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
#features_filename = '/srv/gsfs0/projects/snyder/smartis/highres_microbiome_timecourse_code/barcode_tracking_files/Bacteroides_vulgatus_snps_to_track.txt'
#features_filename = '/srv/gsfs0/projects/snyder/smartis/testfile.txt'
#features_filename = '/srv/gsfs0/projects/snyder/smartis/highres_microbiome_timecourse_code/barcode_tracking_files/Bacteroides_faecis_snps_to_track.txt'
#features_filename = '/srv/gsfs0/projects/snyder/smartis/highres_microbiome_timecourse_code/barcode_tracking_files/Alistipes_onderdonkii_snps_to_track.txt'
corrected = True
core_genes_only = True
external_core_genes=False
disable_filter=True
min_S = 1 # or 1
Pstar = 0.1

# create output filename
features_filename_items = features_filename.split("/")[-1].split(".")
#output_filename = "/srv/gsfs0/projects/snyder/smartis/highres_microbiome_timecourse_code/snp_correlations/"+features_filename_items[0]
#output_filename += ".total.p"
output_filename = 'test_out.txt'


# load candidate SNPs from supplied file

# long snp map
# long snp = (species, "contig|location") tuple (allele is stored separately)
# id = int
longsnp_id_map = {}
id_longsnp_map = []

desired_speciess = set()
snp_file = open(features_filename,"r")
for line in snp_file:
    if line.startswith('#'):
        continue
        
    if line.startswith('>'):
        # a new species
        current_species = line[1:].strip()
        desired_speciess.add(current_species)
        continue
    
    snp_str = line.strip()
    items = snp_str.split("|")
    contig = items[0]
    location = long(items[1])
    longsnp = (contig, location)
    if longsnp not in longsnp_id_map:
        id = len(id_longsnp_map)
        longsnp_id_map[longsnp]=id
        id_longsnp_map.append(longsnp)
        
snp_file.close()

#desired_samples = parse_timecourse_data.highcoverage_samples
#desired_samples = [parse_timecourse_data.highcoverage_end]
desired_samples = parse_timecourse_data.morteza_samples
#for testing
#desired_samples = parse_timecourse_data.morteza_samples[1:3]
#desired_samples = [parse_timecourse_data.highcoverage_antibiotic, parse_timecourse_data.highcoverage_end]
#desired_samples = [parse_timecourse_data.highcoverage_start_2, parse_timecourse_data.highcoverage_antibiotic, parse_timecourse_data.highcoverage_postantibiotic, parse_timecourse_data.highcoverage_end]


#load inferred G(D) from file
f = open('/srv/gsfs0/projects/snyder/smartis/highres_microbiome_timecourse_code/G.pkl', 'r')
sys.stderr.write("Loading G(D)...")
G_all_samples = pickle.load(f)
f.close()

def getG(D, Gdist):
    idx = min(numpy.digitize([D], Gdist[0], right=True)[0], len(Gdist[0])-1)
    return Gdist[1][idx]

gamete_idx_map = {('R','R'): 0, ('R','A'):1, ('A','R'):2, ('A','A'):3  }
allele_idx_map = {'R':0, 'A':1}
# Memory requirements are big, so we have to break things up into smaller batches
# of unassigned genes. The output is then grouped into a single file.
# TODO: Is there a better way to do this?
batch_size = 100000
num_batches = long(ceil(len(id_longsnp_map)*1.0/batch_size))

if num_batches>1.5 and (len(id_longsnp_map)%batch_size)*1.0/batch_size < 0.1:
    # just expand the batch size a little bit
    batch_size = long(ceil(len(id_longsnp_map)*1.0/(num_batches-1)))
    num_batches = long(ceil(len(id_longsnp_map)*1.0/batch_size))

sys.stderr.write("Divided %d snps into %d batches of size %d.\n" % (len(id_longsnp_map), num_batches, batch_size))

snp_barcode_timecourse = {}

for batch in xrange(0,num_batches):
    
    sys.stderr.write("Processing batch %d...\n" % batch)
    
    for i in xrange(0,len(id_longsnp_map)):
        if (i/batch_size) == batch:
            snp_barcode_timecourse[i] = {'all': numpy.array([0,0]), 'longsnps':{}}

    sys.stderr.write("(%d features)\n" % len(snp_barcode_timecourse))
    
    for sample_idx in xrange(0,len(desired_samples)):
        # Process each sample separately
        # (combine results at end)
    
        sample_name = desired_samples[sample_idx]
        
        G_func = G_all_samples[sample_name]

        num_significant_hits = 0 
        
        sys.stderr.write("Processing sample %s...\n" % sample_name)
        #sys.stderr.write("Loading depth map...\n")
        barcode_depth_map = barcode_utils.parse_barcode_depth_map(sample_name,corrected=corrected,min_depth=2.5)
        #sys.stderr.write("Done!\n")
        
        if len(barcode_depth_map)==0:
            continue

        total_barcodes = len(barcode_depth_map) 
        # create barcode->species map
        #sys.stderr.write("Collating species barcodes...\n")
        
        barcode_longsnp_ids_map = {} 
        
        focal_longsnp_barcode_map = {} # Map from longsnp to barcode ids that that feature has
        
        target_longsnp_barcode_fraction_map = {} # map from target longsnp to fraction of barcodes in set (i.e., Mavg*pj from notes)

        #for species_name in pangenome_species:
        for species_name in desired_speciess:
        
            # Make sure barcodes exist for this species at this timepoint.
            if not barcode_utils.barcodes_exist(species_name, sample_name):
                continue
         
            # Load barcodes      
            allele_barcode_map = barcode_utils.parse_allele_barcode_tuples(species_name, sample_name, corrected=corrected)

            for allele in allele_barcode_map:
        
                if not (allele[-1]=='A' or allele[-1]=='R'):
                    # Not a snp, ignore
                    continue
                    
                snp_str = allele[:-2]
                polarization = allele[-1]
                items = snp_str.split("|")
                contig = items[0]
                location = long(items[1])
                longsnp = (contig, location)
                
                #print longsnp, allele, longsnp_id_map.keys()
                
                if longsnp not in longsnp_id_map: # only looking for linkage within the set of SNVs we have. 
                    continue
                
                longsnp_id = longsnp_id_map[longsnp]
                key = (longsnp_id, polarization)
                    
                if longsnp_id in snp_barcode_timecourse:
                    # one of the focal SNPs we are tracking.
                    
                    # Add in the barcodes
                    for barcode_id, barcode_weight in allele_barcode_map[allele]:
                
                        # don't include barcode if too little coverage
                        # Poor man's error correction
                        if barcode_id not in barcode_depth_map: 
                            continue
                        
                        
                        if key not in focal_longsnp_barcode_map:
                            focal_longsnp_barcode_map[key] = set()
                        # Add barcode to list
                        focal_longsnp_barcode_map[key].add(barcode_id)
                    
                # Only look at snps that have barcodes
                if len(allele_barcode_map[allele])<0.5:
                    continue
            
                # estimated # of fragments per barcode
                total_B = 0

                for barcode_id, barcode_weight in allele_barcode_map[allele]:
                
                    # don't include barcode if too little coverage
                    # Poor man's error correction
                    if barcode_id not in barcode_depth_map: 
                        continue
                    
                    total_B += 1
                        
                    if barcode_id not in barcode_longsnp_ids_map:
                        barcode_longsnp_ids_map[barcode_id] = set()
                        #barcode_depth_map[barcode_id] = 0
                    
                    barcode_longsnp_ids_map[barcode_id].add(key)
                    
                target_longsnp_barcode_fraction_map[key] = total_B*1.0/total_barcodes
        
        
        
        #print barcode_longsnp_ids_map[1185020]
        
        # Make sure we loaded some barcodes...
        if len(barcode_longsnp_ids_map)==0:
            continue
                
        total_longsnps = len(target_longsnp_barcode_fraction_map)
    
        log_corrected_Pstar = log(Pstar/total_longsnps/total_longsnps)
                
        #sys.stderr.write("Done! Loaded %d barcodes\n" % len(barcode_longsnp_ids_map))
        #sys.stderr.write("(other count: %d)\n" % total_barcodes)
        #sys.stderr.write("Checking linkage with %d target snps\n" % total_longsnps)
        
        #sys.stderr.write("Looping through target features...\n")

        num_rejected = 0
        
        for feature in focal_longsnp_barcode_map:
            feature_barcodes = set(focal_longsnp_barcode_map[feature])
            # B = total # barcodes for this feature
            B = len(feature_barcodes)
            
            # Get barcode ids that map to this gene    
            barcode_id_G_map = {}
            for barcode_id in feature_barcodes:
                
                # If it doesn't map to anywhere we know of, 
                # skip it (saves some time later)
                if barcode_id not in barcode_longsnp_ids_map:
                    continue
                     
                barcode_id_G_map[barcode_id] = getG(barcode_depth_map[barcode_id], G_func)
            
            num_barcodes = len(barcode_id_G_map)
            
            # only create entries if there are at least two barcodes
            if num_barcodes<1.5:
                continue
            
            focal_longsnp_id, focal_polarization = feature
                        
            #sys.stderr.write("Processing %d barcodes for %s...\n" % (num_barcodes,focal_longsnp_id))
               
            # Counts number of barcodes that map to different genes
            # within species
            longsnp_id_counter = collections.Counter()
            for barcode_id in feature_barcodes:
                longsnp_id_counter.update( barcode_longsnp_ids_map[barcode_id] )
                    
            # create entry for feature
                    
            snp_barcode_timecourse[focal_longsnp_id]['all'][allele_idx_map[focal_polarization]] += B
            # Now add in stuff from gene counter
            for target_key in longsnp_id_counter:
                
                longsnp_id, polarization = target_key
                
                # Observed of shared barcodes
                S = longsnp_id_counter[target_key]
                
                if S>=min_S:
		    target_barcodes = set(focal_longsnp_barcode_map[target_key])
                    logP = 0
		    target_p = target_longsnp_barcode_fraction_map[target_key]
                    for barcode_id in target_barcodes&feature_barcodes:
                        logP += log(target_p*barcode_id_G_map[barcode_id])
                    if logP < log_corrected_Pstar:
                        
                        num_significant_hits += 1
                        
			#sys.stderr.write("Correlation significant at %s...\n" % logP) 
                        # Surprising! Save it! 
                
                        if longsnp_id not in snp_barcode_timecourse[focal_longsnp_id]['longsnps']:
                            snp_barcode_timecourse[focal_longsnp_id]['longsnps'][longsnp_id] = numpy.array([0,0,0,0])
                        
                        gamete = (focal_polarization, polarization)
                        snp_barcode_timecourse[ focal_longsnp_id]['longsnps'][longsnp_id][gamete_idx_map[gamete]] += S
                
                    else:
                        num_rejected+=1
                            
            # Done!
            
        sys.stderr.write("Done! %d significant hits\n" % num_significant_hits)
        sys.stderr.write("Done! %d rejected\n" % num_rejected)

data_to_save = {'id_longsnp_map': id_longsnp_map, 'allele_idx_map': allele_idx_map, 'gamete_idx_map': gamete_idx_map, 'snp_barcode_timecourse': snp_barcode_timecourse}
          
pickle.dump( data_to_save, open( output_filename, "wb" ) )
sys.stderr.write("Done!\n")

