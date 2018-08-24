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
from math import log10,ceil
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

pylab.figure(1,figsize=(5,1.4))

fig = pylab.gcf()

outer_grid  = gridspec.GridSpec(1,3, width_ratios=[1,1,1], wspace=0.1)

depth_axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(depth_axis)

depth_axis.set_ylabel('# barcodes >=x',labelpad=2)
depth_axis.set_xlabel('# mapped reads',labelpad=2)

depth_axis.semilogx([1,1e03],[2,2],'k:')
depth_axis.set_ylim([0,1.05])

depth_axis.spines['top'].set_visible(False)
depth_axis.spines['right'].set_visible(False)
depth_axis.get_xaxis().tick_bottom()
depth_axis.get_yaxis().tick_left()

species_axis = plt.Subplot(fig, outer_grid[1])
fig.add_subplot(species_axis)

species_axis.set_xlabel('# species',labelpad=2)
species_axis.set_ylim([0,1.05])
species_axis.set_yticklabels([])

species_axis.spines['top'].set_visible(False)
species_axis.spines['right'].set_visible(False)
species_axis.get_xaxis().tick_bottom()
species_axis.get_yaxis().tick_left()

fraction_axis = plt.Subplot(fig, outer_grid[2])
fig.add_subplot(fraction_axis)

fraction_axis.set_xlabel('max species %',labelpad=2)
fraction_axis.set_xlim([0,1.05])
fraction_axis.set_ylim([0,1.05])
fraction_axis.set_yticklabels([])

fraction_axis.spines['top'].set_visible(False)
fraction_axis.spines['right'].set_visible(False)
fraction_axis.get_xaxis().tick_bottom()
fraction_axis.get_yaxis().tick_left()


pangenome_species = parse_midas_data.parse_pangenome_species()
#pangenome_species = ['Bacteroides_vulgatus_57955']

# long gene map
# long gene = (species, gene) tuple
# id = int
longgene_id_map = {}
id_longgene_map = []

desired_samples = parse_timecourse_data.morteza_samples
#desired_samples = [parse_timecourse_data.highcoverage_start_1]
for sample_name in desired_samples:

    sys.stderr.write("Processing sample %s...\n" % sample_name)
    
    sys.stderr.write("Loading error corrected barcodes...\n")
    barcode_error_map = barcode_utils.parse_barcode_error_correction_map(sample_name)
    sys.stderr.write("Done!\n")
    
    if len(barcode_error_map)==0:
        sys.stderr.write("No barcodes, continuing with next sample!\n")
        continue
    
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
        
            for original_barcode_id, barcode_weight in allele_barcode_map[allele]:
            
                if original_barcode_id in barcode_error_map:
                    barcode_id = barcode_error_map[original_barcode_id]
                else:
                    barcode_id = original_barcode_id
                    
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
    num_species = []
    for barcode_id in barcode_longgene_weight_map:
        weights = numpy.array(barcode_longgene_weight_map[barcode_id].values())
        
        total_weight = weights.sum()
        
        if total_weight > 2.5:
            total_weights.append( total_weight )
            max_weights.append( weights.max() )
            num_species.append( len(weights) )
        
    total_weights = numpy.array( total_weights )
    max_weights = numpy.array( max_weights )
    max_fractions = max_weights*1.0/total_weights
    
    median_total_weights = numpy.median(total_weights)
    
    print ">=3:", len(max_fractions), total_weights.sum(), median_total_weights, numpy.median(max_fractions)
    
    weight_bins = numpy.logspace(0,4,100)
    weight_bins[0] = -1
    weights = weight_bins[1:] 
    
    ns, dummy_bins = numpy.histogram(total_weights, bins=weight_bins)
    cdf = numpy.cumsum(ns)
    survivals = cdf[-1]-cdf
    line, = depth_axis.step(weights, survivals*1.0/survivals[0],'-',where='mid')
    
    depth_axis.set_xlim([1,1e03])
    
    species_bins = numpy.arange(0,101)-0.5
    species_bins[-1] = 1e09
    species_xs = numpy.arange(0,100)
    
    ns, dummy_bins = numpy.histogram(num_species, bins=species_bins)
    cdf = numpy.cumsum(ns)
    survivals = cdf[-1]-cdf
    species_axis.step(species_xs, survivals*1.0/survivals[0],'-',color=pylab.getp(line,'color'),where='mid')
    species_axis.set_xlim([-1,21])
    
    fraction_bins = numpy.linspace(0,1.01,100)
    fraction_bins[0] = -1
    fractions = fraction_bins[1:]
    ns, dummy_bins = numpy.histogram(max_fractions, bins=fraction_bins)
    cdf = numpy.cumsum(ns)
    survivals = cdf[-1]-cdf
    fraction_axis.step(fractions, survivals*1.0/survivals[0],'-',color=pylab.getp(line,'color'),where='mid')
    
    
fig.savefig(parse_midas_data.analysis_directory+'barcode_statistics.pdf',bbox_inches='tight')
