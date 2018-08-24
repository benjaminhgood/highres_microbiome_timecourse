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

import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10, ceil, exp
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint
import matplotlib.colors as mcolors

mpl.rcParams['font.size'] = 6
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'


####################################################
#
# Set up Figure (3 panels, arranged in 1x3 grid)
#
####################################################

pylab.figure(1) #,figsize=(5,1.4))
fig = pylab.gcf()

pangenome_species = parse_midas_data.parse_pangenome_species()
#pangenome_species = ['Bacteroides_vulgatus_57955']

# long gene map
# long gene = (species, gene) tuple
# id = int
longgene_id_map = {}
id_longgene_map = []

# Run this algorithm separately for each sample. 
# BG: should we combine across samples? Pros: more data. Cons: could have switching across species between timepoints. 
desired_samples = parse_timecourse_data.morteza_samples
desired_samples = parse_timecourse_data.highcoverage_samples
desired_samples = [parse_timecourse_data.highcoverage_start_1, parse_timecourse_data.highcoverage_start_2] 
#desired_samples = [parse_timecourse_data.highcoverage_end]

coverage_bins = numpy.logspace(0-0.15,log10(2**11)-0.15,12)
#coverage_bins[-1] = 1e08

species_bins = numpy.arange(0,21)-0.5
#species_bins[-1] = 1e08

count_matrix = numpy.zeros((len(coverage_bins)-1,len(species_bins)-1))

num_reads = []
num_species = []

for sample_name in desired_samples:

    sys.stderr.write("Processing sample %s...\n" % sample_name)
    
    sys.stderr.write("Loading depth map...\n")
    barcode_depth_map = barcode_utils.parse_barcode_depth_map(sample_name)
    sys.stderr.write("Done!\n")
    
    # create barcode->species map
    sys.stderr.write("Collating species barcodes...\n")
    # first create intermediate data structure:
    # barcode_id->longgene->count
    barcode_longgene_weight_map = {} 
    for species_name in pangenome_species:
        sys.stderr.write("%s...\n" % species_name)
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
    for barcode_id in barcode_longgene_weight_map:
        weights = numpy.array(barcode_longgene_weight_map[barcode_id].values())
        
        total_weight = weights.sum()
     
        total_depth = barcode_depth_map[barcode_id]
        
        if total_depth > 2.5:
        
            # digitize stuff
            coverage_idx = numpy.digitize(total_depth, coverage_bins)-1
            
            if coverage_idx>=count_matrix.shape[0]:
                coverage_idx = count_matrix.shape[0]-1
            
            weights = weights*1.0/total_weight
            
            effective_num_species = exp( - (weights*numpy.log(weights+(weights==0))).sum() )
            
            species_idx = numpy.digitize(effective_num_species, species_bins)-1
            if species_idx>=count_matrix.shape[1]:
                species_idx = count_matrix.shape[1]-1
            
                
            count_matrix[coverage_idx,species_idx] += 1
            
            num_reads.append(total_weight)
            num_species.append(effective_num_species)

x = numpy.log10(numpy.array(num_reads))
y = numpy.array(num_species)

#print count_matrix

new_coverage_bins = numpy.arange(0, 11)*0.3+0.15
#new_coverage_bins = numpy.arange(0,11)-0.5

pylab.hist2d(x, y, bins=[new_coverage_bins, species_bins], norm=colors.LogNorm()) #, cmap='YlGnBu', cmin=10)
pylab.xlabel('log10(# reads / barcode)')
pylab.ylabel('Effective # species / barcode')
#pylab.semilogx([1],[1],'k.')
pylab.colorbar()
fig.savefig(parse_midas_data.analysis_directory+'barcode_statistics.pdf',bbox_inches='tight')
