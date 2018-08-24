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

from numpy.random import shuffle

pangenome_species = parse_midas_data.parse_pangenome_species()
#pangenome_species = ['Bacteroides_vulgatus_57955']

##################
#
# The purpose of this script is to generate a bootstrapped version of the 10x barcodes
# under the "null hypothesis" that genes are unlinked from each other,
# and separately, that SNPs are unlinked from each other. 
#
# It preserves the total number of barcodes (and depth) covered by each item
# It also tries to permute in a way that respects broad classes of depth
#
##################
   

desired_samples = parse_timecourse_data.morteza_samples


for sample_name in desired_samples:

    sys.stderr.write("Processing sample %s...\n" % sample_name)
         
    sys.stderr.write("Loading barcode count map...\n")
    filename = "%s%s/output/all_corrected_barcodes.gz" % (config.barcode_directory, sample_name)
    barcode_file = gzip.GzipFile(filename,"r")
    line = barcode_file.readline() # skip header
    total_barcode_count_map = {}
    for line in barcode_file:
        items = line.split()
        barcode_id = long(items[0])
        barcode_count = long(items[2]) 
        total_barcode_count_map[barcode_id] = barcode_count
    barcode_file.close()
    sys.stderr.write("Done!\n")    
    
    # Used for stratifying barcodes into different coverage bins. Keeps bootstrapped
    # barcodes close-ish to actual values? 
    coverage_bins = numpy.array([0.5, 2.5, 4.5, 6.5, 8.5, 10.5, 50.5, 100.5, 500.5, 1000.5, 1e08])
    coverage_bin_idxs = numpy.arange(0,len(coverage_bins)) # last one is for zero coverage
    sys.stderr.write("Binning count map...\n")
    barcode_bin_idx_map = {}
    for barcode_id in total_barcode_count_map:
        bin_idx = numpy.digitize(total_barcode_count_map[barcode_id], coverage_bins)-1
        barcode_bin_idx_map[barcode_id] = bin_idx
    sys.stderr.write("Done!\n")     
    
    del total_barcode_count_map
    
    # "longgene" is tuple of (species, gene) 
    longgene_id_map = {}
    id_longgene_map = []
    
    # "longsnp" is tuple of (species, snp_allele) 
    longsnp_id_map = {}
    id_longsnp_map = []
    
    gene_master_barcode_list = [[] for bin_idx in coverage_bin_idxs]
    gene_master_allele_list = [[] for bin_idx in coverage_bin_idxs] # (longgene_id, count)
    
    snp_master_barcode_list = [[] for bin_idx in coverage_bin_idxs]
    snp_master_allele_list = [[] for bin_idx in coverage_bin_idxs] # (longsnp_id, count)
        
    # create barcode->species map
    sys.stderr.write("Collating species barcodes...\n")
    for species_name in pangenome_species:
        
        if species_name=='new_species':
            # we don't really need to do this for the new_species things,
            # since we never test linkage from new_gene to new_gene. 
            continue
            
        # Make sure barcodes exist for this timepoint.
        # BG: aside from bugs, shouldn't this always be true? 
        if not barcode_utils.barcodes_exist(species_name, sample_name):
            continue
        
        #sys.stderr.write("%s\n" % species_name)
         
        barcode_filename = "%s%s/output/%s.corrected_barcodes.gz" % (config.barcode_directory, sample_name, species_name)
        barcode_file = gzip.GzipFile(barcode_filename,"r")
        header_line = barcode_file.readline() # header
    
        allele_barcode_map = {}
    
        for line in barcode_file:
            line = line.strip()
            items = line.split("\t")
            allele = items[0].strip()
            allele_barcode_map[allele] = []
        
            if len(items)>1:
        
                barcode_items = items[1].split(",")
                for barcode_item in barcode_items:
                    barcode_subitems = barcode_item.split(":")
                    barcode_id = long(barcode_subitems[0])
                    barcode_weight = long(barcode_subitems[1])
                    
                    allele_barcode_map[allele].append((barcode_id, barcode_weight))
                
        for allele in allele_barcode_map.keys():
        
            if allele.endswith('|A') or allele.endswith('|R'):
                # a SNP allele
                longsnp = (species_name, allele)
                if longsnp not in longsnp_id_map:
                    longsnp_id_map[longsnp] = len(id_longsnp_map)
                    id_longsnp_map.append(longsnp)
                
                longsnp_id = longsnp_id_map[longsnp]
                
                if len(allele_barcode_map[allele])==0:
                    snp_master_allele_list[-1].append((longsnp_id,0))
                else:
                    for barcode_id, barcode_weight in allele_barcode_map[allele]:
                
                        bin_idx = barcode_bin_idx_map[barcode_id]
                    
                        snp_master_barcode_list[bin_idx].append(barcode_id)
                        snp_master_allele_list[bin_idx].append((longsnp_id, barcode_weight))
                
            else:
    
                # form longgene
                longgene = (species_name, allele)
                # update longgene id
                if longgene not in longgene_id_map:
                    longgene_id_map[longgene] = len(id_longgene_map)
                    id_longgene_map.append(longgene)
            
                longgene_id = longgene_id_map[longgene]
    
                if len(allele_barcode_map[allele])==0:
                    gene_master_allele_list[-1].append((longgene_id,0))
                else:
    
                    for barcode_id, barcode_weight in allele_barcode_map[allele]:
                
                        bin_idx = barcode_bin_idx_map[barcode_id]
                    
                        gene_master_barcode_list[bin_idx].append(barcode_id)
                        gene_master_allele_list[bin_idx].append((longgene_id, barcode_weight))
                
    sys.stderr.write("Done!\n")
        
    sys.stderr.write("Shuffling barcodes...\n")
    for bin_idx in coverage_bin_idxs[1:-1]: 
        # don't shuffle the first one, since it's just 1,2-read barcodes
        # don't shuffle the last one, since it's empty barcodes
        #
        # shuffle the barcodes, rather than the alleles, so that alleles stay
        # in species order for output!
        shuffle(gene_master_barcode_list[bin_idx])
        shuffle(snp_master_barcode_list[bin_idx]) 
    sys.stderr.write("Done!\n")               
        
    sys.stderr.write("Writing new species barcodes...\n")
    
    current_gene_idxs = [0 for idx in xrange(0,len(coverage_bin_idxs))]
    current_snp_idxs = [0 for idx in xrange(0,len(coverage_bin_idxs))]
    
    for species_name in pangenome_species:
        
        if species_name=='new_species':
            # we don't really need to do this for the new_species things,
            # since we never test linkage from new_gene to new_gene. 
            continue
            
        # Make sure barcodes exist for this timepoint.
        # BG: aside from bugs, shouldn't this always be true? 
        if not barcode_utils.barcodes_exist(species_name, sample_name):
            continue
        
        #sys.stderr.write("%s\n" % species_name)
        #sys.stderr.write("Genes...\n")
            
        gene_allele_barcode_map = {}
        for bin_idx in coverage_bin_idxs:
            
            if len(gene_master_allele_list[bin_idx])==0:
                continue # nothing to do!
                
            longgene_id = gene_master_allele_list[bin_idx][current_gene_idxs[bin_idx]][0]
            barcode_weight = gene_master_allele_list[bin_idx][current_gene_idxs[bin_idx]][1]
            new_species_name, allele = id_longgene_map[longgene_id] 
            
            num_barcodes=0
            while new_species_name==species_name: # if we're looking at the right species
                # add it to the list!
                
                #num_barcodes+=1
                #if num_barcodes % 1000 == 0:
                #    print num_barcodes
                #    print new_species_name, species_name, longgene_id
                
                if allele not in gene_allele_barcode_map:
                    gene_allele_barcode_map[allele] = []
                
                if len(gene_master_barcode_list[bin_idx])>0:
                    # add barcodes!
                    barcode_id = gene_master_barcode_list[bin_idx][current_gene_idxs[bin_idx]]
                
                    gene_allele_barcode_map[allele].append((barcode_id, barcode_weight))
                    
                # increment idx
                current_gene_idxs[bin_idx] += 1
                if current_gene_idxs[bin_idx]>=len(gene_master_allele_list[bin_idx]):
                    current_gene_idxs[bin_idx]=0
                    new_species_name = 'none'
                else:    
                    longgene_id = gene_master_allele_list[bin_idx][current_gene_idxs[bin_idx]][0]
                    barcode_weight = gene_master_allele_list[bin_idx][current_gene_idxs[bin_idx]][1]
                    new_species_name, allele = id_longgene_map[longgene_id] 
        
        # now output genes
        new_barcode_filename = "%s%s/output/%s.corrected_barcodes.bootstrapped.gz" % (config.barcode_directory, sample_name, species_name)
        new_barcode_file = gzip.GzipFile(new_barcode_filename,"w")
        new_barcode_file.write(header_line) # header
    
        for allele in gene_allele_barcode_map:
            
            barcode_str = ", ".join(["%d:%d" % (barcode_id, barcode_weight) for barcode_id, barcode_weight in gene_allele_barcode_map[allele]])
            output_str = "%s\t%s\n" % (allele, barcode_str)
            new_barcode_file.write(output_str)
                
        # Now do the same thing for snps    
        snp_allele_barcode_map = {}
        snp_allele_barcode_map = {}
        #sys.stderr.write("SNPs...\n")
        for bin_idx in coverage_bin_idxs:
                
            if len(snp_master_allele_list[bin_idx])==0:
                continue # nothing to do!
                
            longsnp_id = snp_master_allele_list[bin_idx][current_snp_idxs[bin_idx]][0]
            barcode_weight = snp_master_allele_list[bin_idx][current_snp_idxs[bin_idx]][1]
            new_species_name, allele = id_longsnp_map[longsnp_id] 
            
            num_barcodes=0
            while new_species_name==species_name: # if we're looking at the right species
                
                #num_barcodes+=1
                #if num_barcodes % 1000 == 0:
                #    print num_barcodes
                #    print new_species_name, species_name, longsnp_id
                
                
                # add it to the list!
                if allele not in snp_allele_barcode_map:
                    snp_allele_barcode_map[allele] = []
                
                if len(snp_master_barcode_list[bin_idx])>0:
                    # add barcodes!
                    barcode_id = snp_master_barcode_list[bin_idx][current_snp_idxs[bin_idx]]
                
                    snp_allele_barcode_map[allele].append((barcode_id, barcode_weight))
                    
                # increment idx
                current_snp_idxs[bin_idx] += 1
                if current_snp_idxs[bin_idx]>=len(snp_master_allele_list[bin_idx]):
                    current_snp_idxs[bin_idx]=0
                    new_species_name = 'none'
                else:    
                    longsnp_id = snp_master_allele_list[bin_idx][current_snp_idxs[bin_idx]][0]
                    barcode_weight = snp_master_allele_list[bin_idx][current_snp_idxs[bin_idx]][1]
                    new_species_name, allele = id_longsnp_map[longsnp_id] 
            
        # now output snps
        for allele in snp_allele_barcode_map:
            
            barcode_str = ", ".join(["%d:%d" % (barcode_id, barcode_weight) for barcode_id, barcode_weight in snp_allele_barcode_map[allele]])
            output_str = "%s\t%s\n" % (allele, barcode_str)
            new_barcode_file.write(output_str)
            
        new_barcode_file.close()
            
    sys.stderr.write("Done!\n")
    
    sys.stderr.write("Purging memory...\n") 
    del longgene_id_map 
    del id_longgene_map 
    del longsnp_id_map 
    del id_longsnp_map
    
    del gene_master_barcode_list 
    del gene_master_allele_list 
    del snp_master_barcode_list 
    del snp_master_allele_list
    sys.stderr.write("Done!\n")
          
sys.stderr.write("Done!\n")

