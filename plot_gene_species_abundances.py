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
species_axis.set_ylabel('# species $\geq f$')

species_axis.spines['top'].set_visible(False)
species_axis.spines['right'].set_visible(False)
species_axis.get_xaxis().tick_bottom()
species_axis.get_yaxis().tick_left()

species_axis.semilogx([1e-06],[0.1],'k.')
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

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(mean_abundances, min_x=1e-06, max_x=1)
#species_axis.step(xs,ns,'-',color='b',alpha=0.5)

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(median_abundances, min_x=1e-06, max_x=1)
species_axis.step(xs,ns,'-',color='b',label='Median')


xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(max_abundances, min_x=1e-06, max_x=1)
species_axis.step(xs,ns,'-',color='b',alpha=0.5,label='Max')

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(min_abundances, min_x=1e-06, max_x=1)
#species_axis.step(xs,ns,'-',color='g')

species_axis.legend(loc='upper right', frameon=False, fontsize=6, numpoints=1, handlelength=1)

####
#
# Now do same thing for genes in pangenome
#
####

# get list of pangenome species
pangenome_species = parse_midas_data.parse_pangenome_species()
#pangenome_species = ['Bacteroides_vulgatus_57955', 'new_species']

gene_median_abundances = []
gene_mean_abundances = []
gene_max_abundances = []

reference_median_abundances = []
reference_mean_abundances = []
reference_max_abundances = []

assembly_median_abundances = []
assembly_mean_abundances = []
assembly_max_abundances = []


for species_name in pangenome_species:
        
    sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
    gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix =     parse_midas_data.parse_pangenome_data(species_name)
    sys.stderr.write("Done!\n")

    species_times, species_time_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, gene_samples)
    desired_gene_samples = numpy.array(gene_samples)[species_time_idxs]

    print gene_depth_matrix.shape
    print desired_gene_samples
    print species_time_idxs

    gene_depth_matrix = gene_depth_matrix[:,species_time_idxs]

    gene_sample_idx_map = parse_midas_data.calculate_sample_idx_map(desired_gene_samples, desired_samples)
    
    gene_coverage_to_abundance_factor = numpy.array([coverage_to_abundance_factor[gene_sample_idx_map[i]] for i in xrange(0,gene_depth_matrix.shape[1])])
        
    gene_abundances = gene_depth_matrix * gene_coverage_to_abundance_factor


    gene_median_abundances.extend( numpy.median(gene_abundances,axis=1) )
    gene_mean_abundances.extend( gene_abundances.mean(axis=1) )
    gene_max_abundances.extend( gene_abundances.max(axis=1) )
    
    if species_name!='new_species':
        
        reference_median_abundances.extend( numpy.median(gene_abundances,axis=1) )
        reference_mean_abundances.extend( gene_abundances.mean(axis=1) )
        reference_max_abundances.extend( gene_abundances.max(axis=1) )
    
    else:
        
        assembly_median_abundances.extend( numpy.median(gene_abundances,axis=1) )
        assembly_mean_abundances.extend( gene_abundances.mean(axis=1) )
        assembly_max_abundances.extend( gene_abundances.max(axis=1) )
    
    
gene_median_abundances = numpy.array(gene_median_abundances)
gene_mean_abundances = numpy.array(gene_mean_abundances)
gene_max_abundances = numpy.array(gene_max_abundances)


reference_median_abundances = numpy.array(reference_median_abundances)
reference_mean_abundances = numpy.array(reference_mean_abundances)
reference_max_abundances = numpy.array(reference_max_abundances)


assembly_median_abundances = numpy.array(assembly_median_abundances)
assembly_mean_abundances = numpy.array(assembly_mean_abundances)
assembly_max_abundances = numpy.array(assembly_max_abundances)


abundance_bins = numpy.logspace(-6,0,300)
abundance_bins[0] = -1
abundances = abundance_bins[1:]

ns, dummy_bins = numpy.histogram(reference_median_abundances, bins=abundance_bins)
cdf = numpy.cumsum(ns)
survivals = cdf[-1]-cdf
gene_axis.step(abundances, survivals,'-',color='b',label='midas_db')

ns, dummy_bins = numpy.histogram(reference_max_abundances, bins=abundance_bins)
cdf = numpy.cumsum(ns)
survivals = cdf[-1]-cdf
gene_axis.step(abundances, survivals,'-',color='b',alpha=0.5)


ns, dummy_bins = numpy.histogram(assembly_median_abundances, bins=abundance_bins)
cdf = numpy.cumsum(ns)
survivals = cdf[-1]-cdf
gene_axis.step(abundances, survivals,'-',color='r',label='extra')

ns, dummy_bins = numpy.histogram(assembly_max_abundances, bins=abundance_bins)
cdf = numpy.cumsum(ns)
survivals = cdf[-1]-cdf
gene_axis.step(abundances, survivals,'-',color='r',alpha=0.5)


gene_axis.set_ylim([1,survivals[0]])


gene_axis.legend(loc='lower left', frameon=False, fontsize=6, numpoints=1, handlelength=1)

fig.savefig(parse_midas_data.analysis_directory+'species_abundances.pdf',bbox_inches='tight')

