import sys
import pylab
import numpy
import parse_midas_data
import parse_timecourse_data
import stats_utils

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
# Set up Figure (2 panels, arranged in 2x1 grid)
#
####################################################

pylab.figure(1,figsize=(3.42,4))
fig = pylab.gcf()
# make three panels panels
outer_grid  = gridspec.GridSpec(2,1,height_ratios=[1,1], hspace=0.1)

species_axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(species_axis)
species_axis.set_ylabel('# genes $\geq f$')

species_axis.spines['top'].set_visible(False)
species_axis.spines['right'].set_visible(False)
species_axis.get_xaxis().tick_bottom()
species_axis.get_yaxis().tick_left()

species_axis.loglog([1e-06],[0.1],'k.')
species_axis.set_xlim([1e-03,1])
species_axis.set_xticklabels([])


# Genes from pangenome
gene_axis = plt.Subplot(fig, outer_grid[1])
fig.add_subplot(gene_axis)

gene_axis.set_ylabel('# genes $\geq f$')

gene_axis.spines['top'].set_visible(False)
gene_axis.spines['right'].set_visible(False)
gene_axis.get_xaxis().tick_bottom()
gene_axis.get_yaxis().tick_left()

gene_axis.loglog([1e-06],[0.1],'k.')
gene_axis.set_xlim([1e-03,1])
gene_axis.set_xlabel('Effective relative abundance, $f$')
  

sample_time_map = parse_timecourse_data.parse_sample_time_map()
    
species_coverage_matrix, samples, species = parse_midas_data.parse_global_marker_gene_coverages()

species_times, species_time_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, samples)
desired_samples = numpy.array(samples)[species_time_idxs]

species_coverage_matrix = species_coverage_matrix[:,species_time_idxs]

coverage_to_abundance_factor = 1.0/species_coverage_matrix.sum(axis=0)

species_abundance_matrix = species_coverage_matrix * coverage_to_abundance_factor[None,:] 

mean_abundances = []
median_abundances = []
max_abundances = []
min_abundances = []
mean_abundances.extend( species_abundance_matrix.mean(axis=1) ) 
median_abundances.extend( numpy.median(species_abundance_matrix,axis=1) ) 
max_abundances.extend( species_abundance_matrix.max(axis=1) ) 
min_abundances.extend( species_abundance_matrix.min(axis=1) )
    
mean_abundances = numpy.array(mean_abundances)
median_abundances = numpy.array(median_abundances)
max_abundances = numpy.array(max_abundances)
min_abundances = numpy.array(min_abundances)

####
#
# Now do same thing for genes in pangenome
#
####

# get list of pangenome species
pangenome_species = parse_midas_data.parse_pangenome_species()
#pangenome_species = ['Bacteroides_vulgatus_57955', 'new_species']


reference_abundances = {i:[] for i in xrange(0,len(coverage_to_abundance_factor))}
assembly_abundances = {i:[] for i in xrange(0,len(coverage_to_abundance_factor))}


for species_name in pangenome_species:
        
    sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
    gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix =     parse_midas_data.parse_pangenome_data(species_name)
    sys.stderr.write("Done!\n")

    species_times, species_time_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, gene_samples)
    desired_gene_samples = numpy.array(gene_samples)[species_time_idxs]

    gene_depth_matrix = gene_depth_matrix[:,species_time_idxs]

    gene_sample_idx_map = parse_midas_data.calculate_sample_idx_map(desired_gene_samples, desired_samples)
    
    
    idxs = numpy.array([gene_sample_idx_map[i] for i in xrange(0,gene_depth_matrix.shape[1])])
    
    gene_coverage_to_abundance_factor = numpy.array([coverage_to_abundance_factor[gene_sample_idx_map[i]] for i in xrange(0,gene_depth_matrix.shape[1])])
        
    gene_abundances = gene_depth_matrix * gene_coverage_to_abundance_factor

    if species_name!='new_species':
        for i in xrange(0,len(idxs)):
            reference_abundances[idxs[i]].extend( gene_abundances[:,i] )
            
    else:
        for i in xrange(0,len(idxs)):
            assembly_abundances[idxs[i]].extend( gene_abundances[:,i] )
   
for i in sorted(reference_abundances.keys()):
    reference_abundances[i] = numpy.array(reference_abundances[i])
    assembly_abundances[i] = numpy.array(assembly_abundances[i])
    
    
    abundance_bins = numpy.logspace(-6,0,300)
    abundance_bins[0] = -1
    abundances = abundance_bins[1:]

    ns, dummy_bins = numpy.histogram(reference_abundances[i], bins=abundance_bins)
    cdf = numpy.cumsum(ns)
    survivals = cdf[-1]-cdf
    
    line, = species_axis.step(abundances, survivals,'-')
    color = pylab.getp(line,'color')

    ns, dummy_bins = numpy.histogram(assembly_abundances[i], bins=abundance_bins)
    cdf = numpy.cumsum(ns)
    survivals = cdf[-1]-cdf
    
    line, = gene_axis.step(abundances, survivals,'-',color=color)
    

species_axis.set_ylim([1,1e06])
gene_axis.set_ylim([1,1e06])
    

fig.savefig(parse_midas_data.analysis_directory+'gene_temporal_abundances.pdf',bbox_inches='tight')

