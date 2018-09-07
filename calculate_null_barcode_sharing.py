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


corrected = True
core_genes_only = True
external_core_genes=False
disable_filter=True
min_S = 1 # or 1
Pstar = 0.1


Dbins = numpy.logspace(0,3,20)-1e-06
Ds = numpy.array(list(Dbins[1:]))
Dbins[-1] = 1e09

#desired_samples = parse_timecourse_data.highcoverage_samples
#desired_samples = [parse_timecourse_data.highcoverage_end]
#desired_samples = parse_timecourse_data.morteza_samples
#desired_samples = parse_timecourse_data.morteza_samples[1:3]
desired_samples = [parse_timecourse_data.highcoverage_antibiotic, parse_timecourse_data.highcoverage_end]
#desired_samples = [parse_timecourse_data.highcoverage_start_2, parse_timecourse_data.highcoverage_antibiotic, parse_timecourse_data.highcoverage_postantibiotic, parse_timecourse_data.highcoverage_end]

pylab.figure(1)
pylab.xlabel('D')
pylab.ylabel('G(D)')
        
    
for sample_idx in xrange(0,len(desired_samples)):
    
    # Process each sample separately
    
    sample_name = desired_samples[sample_idx]

    sys.stderr.write("Processing sample %s...\n" % sample_name)
         
    # map from species_name to numeric species id
    # (used to save memory)
    species_id_map = {pangenome_species[i]:i for i in xrange(0,len(pangenome_species))}
    # reverse map from species id to species_name
    id_species_map = pangenome_species
    
    # map from species name to list of barcodes that map to SNVs in that species
    # organized as map from barcode id -> counts 
    # (to handle barcodes that map to multiple SNVs in a given species)
    species_barcode_map = {}
    
    # map from barcode id to total coverage of barcode, Dmu
    barcode_depth_map = {}
    # map from barcode id to total number of SNVs covered in each species
    # organized as a map from species -> n_s
    barcode_n_map = {}
    
    # create barcode->species map
    sys.stderr.write("Collating species barcodes...\n")
    for species_name in pangenome_species:
        
        species_id = species_id_map[species_name]
        
        species_barcode_map[species_name] = {}
        
        # 'new_species' is set of de-novo assembled genes
        # do not use these for now
        if species_name=='new_species':
            continue
        
        sys.stderr.write("%s...\n" % species_name)
                
        # Make sure barcodes exist for this species at this timepoint.
        if not barcode_utils.barcodes_exist(species_name, sample_name):
            continue
         
        # Load barcodes      
        allele_barcode_map = barcode_utils.parse_allele_barcode_tuples(species_name, sample_name,corrected=corrected)

        for allele in allele_barcode_map:
        
            if not (allele[-1]=='A' or allele[-1]=='R'):
                # Not a snp, ignore
                continue
                    
            # Add in the barcodes
            for barcode_id, barcode_weight in allele_barcode_map[allele]:
                
                if barcode_id not in barcode_depth_map:
                    barcode_depth_map[barcode_id] = 0
                    barcode_n_map[barcode_id] = {}
                    
                barcode_depth_map[barcode_id] += barcode_weight
                
                if species_id not in barcode_n_map[barcode_id]:
                    barcode_n_map[barcode_id][species_id] = 0
                    
                barcode_n_map[barcode_id][species_id] += 1
                
                if barcode_id not in species_barcode_map[species_name]:
                    species_barcode_map[species_name][barcode_id] = 0
                
                species_barcode_map[species_name][barcode_id] += 1
    
    sys.stderr.write("Done!\n")
    
    
    sys.stderr.write("Postprocessing barcodes...\n")
    # Throw out bad barcodes (poor man's error correction)
    filtered_barcode_depth_map = {}
    for barcode_id in barcode_depth_map:
        if barcode_depth_map[barcode_id] > 2.5:
            filtered_barcode_depth_map[barcode_id] = barcode_depth_map[barcode_id]
                
    total_num_barcodes = len(filtered_barcode_depth_map)
    
    # Calculate map from species id to Q_s
    species_id_Q_map = {}
    for species_name in species_barcode_map:
        species_id = species_id_map[species_name]
        species_id_Q_map[species_id] = 0
        
        for barcode_id in species_barcode_map[species_name]:
            
            if barcode_id not in filtered_barcode_depth_map:
                continue
            
            species_id_Q_map[species_id] += species_barcode_map[species_name][barcode_id]*1.0/total_num_barcodes
    
    del barcode_depth_map
    del species_barcode_map
    
    Qtot = 0
    Qdiag = 0
    for species_id in species_id_Q_map:
        Qs = species_id_Q_map[species_id]
        Qtot += Qs
        Qdiag += Qs*Qs
        
    Qvar = Qtot*Qtot-Qdiag
    
    sys.stderr.write("Done!\n")
    
    sys.stderr.write("Estimating G(D)...\n")
    
    binned_Gs = numpy.zeros_like(Ds)
    binned_ns = numpy.zeros_like(Ds)
    binned_vars = numpy.zeros_like(Ds)
    binned_var_counts = numpy.zeros_like(Ds)
    
    for barcode_id in filtered_barcode_depth_map:
        
        D = filtered_barcode_depth_map[barcode_id]
        
        # find the bin(s) that this depth value belongs to
        bin_idxs = ((D<=Dbins[1:])*(D>Dbins[:-1]))
        
        ntot = 0
        ndiag = 0
        for species_id in barcode_n_map[barcode_id]:
            ns = barcode_n_map[barcode_id][species_id] 
            ntot += ns
            ndiag += ns*ns
            
        var = ntot*ntot-ndiag
        binned_vars[bin_idxs] += var
        binned_var_counts[bin_idxs] += 1
        
        for species_id in barcode_n_map[barcode_id]:
            Qs = species_id_Q_map[species_id]
            ns = barcode_n_map[barcode_id][species_id] 
            dg = (ntot-ns)*1.0/(Qtot-Qs)
            binned_Gs[bin_idxs] += ns*dg
            binned_ns[bin_idxs] += ns
            
    
    binned_vars = binned_vars*1.0/(binned_var_counts+(binned_var_counts==0))
    binned_Gs = binned_Gs/(binned_ns+(binned_ns==0))
    

    other_binned_Gs = numpy.sqrt( binned_vars / Qvar )
    other_avg_G = (binned_var_counts * other_binned_Gs).sum()*1.0/binned_var_counts.sum()
    avg_G = (binned_Gs*binned_var_counts).sum()*1.0/binned_var_counts.sum()
    
    sys.stderr.write("Done!\n")
    
    
    line, = pylab.plot(Ds[binned_ns>0], binned_Gs[binned_ns>0],'v-',markersize=3,alpha=0.5,markeredgewidth=0)
    color = pylab.getp(line,'color')
    pylab.loglog(Ds[binned_ns>0], other_binned_Gs[binned_ns>0],'^-',markersize=3,color=color,alpha=0.5,markeredgewidth=0)
    #pylab.loglog(Ds, numpy.ones_like(Ds)*avg_G,'b:')
    #pylab.loglog(Ds, numpy.ones_like(Ds)*other_avg_G,'g:')
    
pylab.savefig('test.pdf',bbox_inches='tight')