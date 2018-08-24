import tarfile, bz2, cPickle

# MIDAS STuff
import parse_midas_data
import parse_timecourse_data as parse_sample_data
import parse_timecourse_data
import numpy
from numpy.random import shuffle
import stats_utils
import diversity_utils
import sys

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
parser.add_argument("--species", help="Name of specific species to run code on")
parser.add_argument("-o", "--outdir", help="Where to write output file",metavar="DIR")
parser.add_argument("-Lmax", "--downsample", type=int, help="Where to write output file",default=1e08)
parser.add_argument("-fstar", "--freq-threshold", type=float, help="Frequency has to exceed to be included",default=0.2)
parser.add_argument("--fraction-covered", type=float, help="Fraction of timepoints with sufficient coverage",default=0.5)


args = parser.parse_args()

debug = args.debug
chunk_size = args.chunk_size
species_name = args.species
outdir = args.outdir
max_snps = args.downsample
fstar = args.freq_threshold
fraction_covered = args.fraction_covered

snp_file = outdir+"/"+species_name+".hopach.txt"

snp_samples = parse_timecourse_data.morteza_samples
sample_size = len(snp_samples)
sys.stderr.write("Proceeding with %d temporal samples!\n" % sample_size)

snp_alignment = [] # (construct a # sites x # samples x # bases (4) array) 

final_line_number = 0
while final_line_number >= 0:
    
    sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
    dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_samples=snp_samples, chunk_size=chunk_size,initial_line_number=final_line_number)
    sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
    snp_samples = dummy_samples
        
    for gene_name in allele_counts_map.keys():
        for var_type in allele_counts_map[gene_name].keys():
            locations = allele_counts_map[gene_name][var_type]['locations']
            allele_counts = allele_counts_map[gene_name][var_type]['alleles']
            if len(allele_counts)==0:
                continue
                    
            depths = allele_counts.sum(axis=2)
            freqs = allele_counts[:,:,0]*1.0/(depths+(depths==0))
                
            for snp_idx in xrange(0,len(locations)):
                
                insufficient_coverage = ((depths[snp_idx,:]>0).sum() < fraction_covered*depths.shape[1])
                low_frequency = ((freqs[snp_idx]<fstar).all() or (freqs[snp_idx]>(1-fstar)).all())
                
                if insufficient_coverage or low_frequency:
                    continue
                 
                four_base_counts = [numpy.hstack([allele_counts[snp_idx,sample_idx,:], [0,0]]) for sample_idx in xrange(0,depths.shape[1])]   
                snp_alignment.append( four_base_counts )
         
    
    shuffle(snp_alignment)
    
    if len(snp_alignment) > max_snps:
        snp_alignment = snp_alignment[0:max_snps]   
    
    snp_alignment = numpy.array(snp_alignment)
    
    file = open(snp_file,"w")   
    for snp_idx in xrange(0,snp_alignment.shape[0]):
        
        A = snp_alignment[snp_idx][:,0]
        D = snp_alignment[snp_idx][:,:].sum(axis=1)
        
        fs = A*1.0/(D+(D==0))
        
        file.write("\t".join([str(f) for f in fs]))
        file.write("\n")
    file.close()
    