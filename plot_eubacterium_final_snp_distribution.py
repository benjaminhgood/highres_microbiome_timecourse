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


args = parser.parse_args()

debug = args.debug
chunk_size = args.chunk_size
species_name = 'Eubacterium_eligens_61678'

snp_samples = [parse_timecourse_data.highcoverage_antibiotic, parse_timecourse_data.highcoverage_postantibiotic, parse_timecourse_data.highcoverage_end]
sample_size = len(snp_samples)
sys.stderr.write("Proceeding with %d temporal samples!\n" % sample_size)

final_freqs = []

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
                
                
                insufficient_coverage = (depths[snp_idx,:]>10).sum() < 3
                if insufficient_coverage:
                    continue
                    
                antibiotic_freq = freqs[snp_idx,0]
                postantibiotic_freq = freqs[snp_idx,1]
                final_freq = freqs[snp_idx,2]
                
                if antibiotic_freq>0.8:
                    # re-polarize
                    antibiotic_freq = 1-antibiotic_freq
                    postantibiotic_freq = 1-postantibiotic_freq
                    final_freq = 1-final_freq
                    
                if antibiotic_freq<0.2 and postantibiotic_freq>0.8:
                    final_freqs.append(final_freq)
                    
import pylab
bins = numpy.linspace(0,1,40)
pylab.hist(final_freqs,bins=bins)
pylab.savefig('final_frequency_distribution.pdf',bbox_inches='tight')
 