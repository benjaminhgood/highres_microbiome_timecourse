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


sample_time_map = parse_timecourse_data.parse_sample_time_map()
    
species_coverage_matrix, samples, species = parse_midas_data.parse_global_marker_gene_coverages()

species_times, species_time_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, samples)
desired_samples = numpy.array(samples)[species_time_idxs]
species_coverage_matrix = species_coverage_matrix[:,species_time_idxs]
ts = species_times
   
####################################################
#
# Set up Figure (2 panels, arranged in 2x1 grid)
#
####################################################

# Set up figure
fig = plt.figure(figsize=(5, 1.7))

# Set up grids to hold figure panels
outer_grid = gridspec.GridSpec(1, 2, width_ratios=[4,1], wspace=0.1)

trajectory_grid = gridspec.GridSpecFromSubplotSpec(2, 1, height_ratios=[0.3,1],
                subplot_spec=outer_grid[0], hspace=0.05) #, hspace=0.08)

alpha_axis = plt.Subplot(fig, trajectory_grid[0])
fig.add_subplot(alpha_axis)

alpha_axis.set_ylabel('$\\beta$ diversity')
 
alpha_axis.spines['top'].set_visible(False)
alpha_axis.spines['right'].set_visible(False)
alpha_axis.get_xaxis().tick_bottom()
alpha_axis.get_yaxis().tick_left()

freq_axis = plt.Subplot(fig, trajectory_grid[1])
fig.add_subplot(freq_axis)

freq_axis.set_ylabel('Genes')
freq_axis.set_xlabel('Day')
#freq_axis.set_ylim([1e-04,1])
freq_axis.set_ylim([0,1])
freq_axis.set_yticks([])
 
freq_axis.spines['top'].set_visible(False)
freq_axis.spines['right'].set_visible(False)
freq_axis.get_xaxis().tick_bottom()
freq_axis.get_yaxis().tick_left()

# Fill in markers
alpha_axis.plot([parse_timecourse_data.antibiotic_start, parse_timecourse_data.antibiotic_end],[0.31,0.31],'k-')

alpha_axis.plot([parse_timecourse_data.lyme_infection, parse_timecourse_data.antibiotic_start],[0.24,0.24],'k-')

alpha_axis.plot([parse_timecourse_data.hrv_infection-1, parse_timecourse_data.hrv_infection+1], [0.24,0.24], 'k-')

alpha_axis.text((parse_timecourse_data.lyme_infection+parse_timecourse_data.antibiotic_start)/2.0,0.26,'Lyme',ha='center')
alpha_axis.text((parse_timecourse_data.antibiotic_start+parse_timecourse_data.antibiotic_end)/2.0,0.33,'Doxycycline',ha='center')
alpha_axis.text(parse_timecourse_data.hrv_infection,0.26,'HRV',ha='center')

alpha_axis.set_ylim([0,0.4])
alpha_axis.set_yticks([0,0.1,0.2,0.3,0.4])

alpha_axis.set_xlim([ts[0],ts[-1]+1])
freq_axis.set_xlim([ts[0],ts[-1]+1])
alpha_axis.set_xticklabels([])

legend_axis = plt.Subplot(fig, outer_grid[1])
fig.add_subplot(legend_axis)

legend_axis.set_ylim([0,1])
legend_axis.set_xlim([0,1])

legend_axis.spines['top'].set_visible(False)
legend_axis.spines['right'].set_visible(False)
legend_axis.spines['left'].set_visible(False)
legend_axis.spines['bottom'].set_visible(False)

legend_axis.set_xticks([])
legend_axis.set_yticks([])  

###########
#
# Do calculation
#
####

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

rare_frequency_threshold = 1e-02

midas_abundances = numpy.zeros_like(ts)
midas_common_abundances = numpy.zeros_like(midas_abundances)
midas_rare_abundances = numpy.zeros_like(midas_abundances)

denovo_abundances = numpy.zeros_like(midas_abundances)
denovo_common_abundances = numpy.zeros_like(midas_abundances)
denovo_rare_abundances = numpy.zeros_like(midas_abundances)

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
    gene_sample_idxs = numpy.array([gene_sample_idx_map[idx] for idx in xrange(0,len(desired_gene_samples))])
    
    gene_coverage_to_abundance_factor = numpy.array([abundance_per_coverage_map[sample] for sample in desired_gene_samples])    
    gene_abundances = gene_depth_matrix * gene_coverage_to_abundance_factor
    
    total_gene_abundances = gene_abundances.sum(axis=0)
    rare_gene_abundances = (gene_abundances*(gene_abundances<rare_frequency_threshold)).sum(axis=0)
    common_gene_abundances = (gene_abundances*(gene_abundances>=rare_frequency_threshold)).sum(axis=0)
    
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
        
        midas_abundances[gene_sample_idxs] += total_gene_abundances
        midas_common_abundances[gene_sample_idxs] += common_gene_abundances
        midas_rare_abundances[gene_sample_idxs] += rare_gene_abundances
        
    
    else:
        
        assembly_median_abundances.extend( numpy.median(gene_abundances,axis=1) )
        assembly_mean_abundances.extend( gene_abundances.mean(axis=1) )
        assembly_max_abundances.extend( gene_abundances.max(axis=1) )
        
        denovo_abundances[gene_sample_idxs] += total_gene_abundances
        denovo_common_abundances[gene_sample_idxs] += common_gene_abundances
        denovo_rare_abundances[gene_sample_idxs] += rare_gene_abundances
        

total_abundances = denovo_abundances+midas_abundances    


lower = numpy.zeros_like(ts)
upper = denovo_rare_abundances/total_abundances
freq_axis.fill_between(ts,lower,upper,color='#f4a582')

lower = upper+0
upper = upper+denovo_common_abundances/total_abundances
freq_axis.fill_between(ts,lower,upper,color='#ca0020')

lower = upper+0
upper = upper+midas_common_abundances/total_abundances
freq_axis.fill_between(ts,lower,upper,color='#0571b0')

lower = upper+0
upper = upper+midas_rare_abundances/total_abundances
freq_axis.fill_between(ts,lower,upper,color='#92c5de') 

legend_axis.plot([-2,-1],[-2,-1],'s',markersize=6,markeredgewidth=0.0,color='#0571b0', label='Reference')
legend_axis.plot([-2,-1],[-2,-1],'s',markersize=6,markeredgewidth=0.0,color='#92c5de', label='Reference (<1%)')
legend_axis.plot([-2,-1],[-2,-1],'s',markersize=6,markeredgewidth=0.0,color='#ca0020', label='Assembled')
legend_axis.plot([-2,-1],[-2,-1],'s',markersize=6,markeredgewidth=0.0,color='#f4a582', label='Assembled (<1%)')

legend_axis.legend(loc='lower right',frameon=False,fontsize=5,numpoints=1,ncol=1,handlelength=1)   

fig.savefig(parse_midas_data.analysis_directory+'pangenome_abundances.pdf',bbox_inches='tight')

