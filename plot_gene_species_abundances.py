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
species_axis.set_xlim([3e-03,3])
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
gene_axis.set_xlim([3e-03,3])
gene_axis.set_xlabel('Effective relative abundance, $f$')
  

sample_time_map = parse_timecourse_data.parse_sample_time_map()
    
species_coverage_matrix, samples, species = parse_midas_data.parse_global_marker_gene_coverages()

species_times, species_time_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, samples)
desired_samples = numpy.array(samples)[species_time_idxs]

species_coverage_matrix = species_coverage_matrix[:,species_time_idxs]

coverage_to_abundance_factor = 1.0/species_coverage_matrix.sum(axis=0)

species_abundance_matrix = species_coverage_matrix * coverage_to_abundance_factor[None,:] 

species_freq_map = {}
for species_idx in xrange(0,len(species)):
    species_name = species[species_idx]
    species_freq_map[species_name] = species_abundance_matrix[species_idx,:]

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

gene_median_abundances = []
gene_mean_abundances = []
gene_max_abundances = []

reference_median_abundances = []
reference_mean_abundances = []
reference_max_abundances = []

assembly_median_abundances = []
assembly_mean_abundances = []
assembly_max_abundances = []

# Figure out how to convert from coverage to abundance in pangenome data
sys.stderr.write("Calculating coverage->abundance conversion factors...\n")
sample_coverage_abundance_map = {} # sample : [total_marker_coverage, total_marker_abundance]
for species_name in pangenome_species:
    if species_name=='new_species':
        continue
    
    # Load data    
    gene_samples, marker_coverages =  parse_midas_data.parse_pangenome_marker_data(species_name)
    
    # Put it in right temporal order
    species_times, species_time_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, gene_samples)
    
    if len(species_times)==0:
        continue

    desired_gene_samples = numpy.array(gene_samples)[species_time_idxs]
    marker_coverages = marker_coverages[species_time_idxs]
    
    gene_sample_idx_map = parse_midas_data.calculate_sample_idx_map(desired_gene_samples, desired_samples)
    gene_coverage_to_abundance_factor = numpy.array([coverage_to_abundance_factor[gene_sample_idx_map[i]] for i in xrange(0,len(desired_gene_samples))])
    species_freqs = numpy.array([species_freq_map[species_name][gene_sample_idx_map[i]] for i in xrange(0,len(desired_gene_samples))])
    
    for sample_idx in xrange(0,len(desired_gene_samples)):
        sample = desired_gene_samples[sample_idx]
        if sample not in sample_coverage_abundance_map:
            sample_coverage_abundance_map[sample] = [0,0]
        
        sample_coverage_abundance_map[sample][0] += marker_coverages[sample_idx]
        sample_coverage_abundance_map[sample][1] += species_freqs[sample_idx]

# Make final conversion factors
abundance_per_coverage_map = {}
for sample in sample_coverage_abundance_map:
    abundance_per_coverage_map[sample] = sample_coverage_abundance_map[sample][1]/sample_coverage_abundance_map[sample][0]       
sys.stderr.write("Done!\n") 

# Now loop through all species 
# For debugging purposes
#pangenome_species = ['Bacteroides_vulgatus_57955', 'new_species']
for species_name in pangenome_species:
        
    sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
    gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix =     parse_midas_data.parse_pangenome_data(species_name)
    sys.stderr.write("Done!\n")

    
    species_times, species_time_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, gene_samples)
    
    if len(species_times)==0:
        continue

    # put them in the right order
    desired_gene_samples = numpy.array(gene_samples)[species_time_idxs]
    gene_depth_matrix = gene_depth_matrix[:,species_time_idxs]
    marker_coverages = marker_coverages[species_time_idxs]
    
    gene_sample_idx_map = parse_midas_data.calculate_sample_idx_map(desired_gene_samples, desired_samples)
    
    gene_coverage_to_abundance_factor = numpy.array([abundance_per_coverage_map[sample] for sample in desired_gene_samples])    
    gene_abundances = gene_depth_matrix * gene_coverage_to_abundance_factor
    
    gene_median_abundances.extend( numpy.median(gene_abundances,axis=1) )
    gene_mean_abundances.extend( gene_abundances.mean(axis=1) )
    gene_max_abundances.extend( gene_abundances.max(axis=1) )
    
    if species_name!='new_species':
        
        marker_gene_abundances = marker_coverages * gene_coverage_to_abundance_factor
    
        #print species_name
        #print marker_gene_abundances
        #print species_freq_map[species_name]

        
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


abundance_bins = numpy.logspace(-5,1,300)
abundances = numpy.array(abundance_bins[1:],copy=True)

abundance_bins[0] = -1
abundance_bins[-1] = 1e+08

ns, dummy_bins = numpy.histogram(reference_median_abundances, bins=abundance_bins)
cdf = numpy.cumsum(ns)
survivals = cdf[-1]-cdf
abundance_cdf = numpy.cumsum(ns*abundances)

#gene_axis.step(abundances, survivals,'-',color='b',label='midas_db')
gene_axis.step(abundances, abundance_cdf/abundance_cdf[-1],'-',color='b',label='midas_db')

ns, dummy_bins = numpy.histogram(reference_max_abundances, bins=abundance_bins)
cdf = numpy.cumsum(ns)
survivals = cdf[-1]-cdf
gene_axis.step(abundances, survivals,'-',color='b',alpha=0.5)


ns, dummy_bins = numpy.histogram(assembly_median_abundances, bins=abundance_bins)
cdf = numpy.cumsum(ns)

abundance_cdf = numpy.cumsum(ns*abundances)

survivals = cdf[-1]-cdf
#gene_axis.step(abundances, survivals,'-',color='r',label='extra')
gene_axis.step(abundances, abundance_cdf/abundance_cdf[-1],'-',color='r',label='extra')
ns, dummy_bins = numpy.histogram(assembly_max_abundances, bins=abundance_bins)
cdf = numpy.cumsum(ns)
survivals = cdf[-1]-cdf
#gene_axis.step(abundances, survivals,'-',color='r',alpha=0.5)


#gene_axis.set_ylim([1,survivals[0]])
gene_axis.set_ylim([1e-03,1])

gene_axis.legend(loc='lower left', frameon=False, fontsize=6, numpoints=1, handlelength=1)

fig.savefig(parse_midas_data.analysis_directory+'species_abundances.pdf',bbox_inches='tight')

