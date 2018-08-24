import matplotlib  
matplotlib.use('Agg') 
import parse_midas_data
import parse_timecourse_data
import bacterial_phylogeny_utils
import pylab
import sys
import numpy
import hmp_utils
from math import log10, fabs, log
from scipy.stats import gaussian_kde

import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint

mpl.rcParams['font.size'] = 7
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

species_abundance_distribution_map = hmp_utils.parse_species_abundance_distributions()

species_coverage_matrix, samples, species = parse_midas_data.parse_global_marker_gene_coverages()
sample_time_map = parse_timecourse_data.parse_sample_time_map()
ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, samples)
samples = numpy.array(samples)[sample_idxs]
species_coverage_matrix = species_coverage_matrix[:,sample_idxs]
total_coverage = species_coverage_matrix.sum(axis=0)
species_freq_matrix = numpy.clip(species_coverage_matrix*1.0/total_coverage,0, 2)    

initial_idx, antibiotic_idx, final_idx = parse_timecourse_data.get_initial_antibiotic_final_idxs(samples)

sample_list = list(samples)

initial_idx_1 = sample_list.index(parse_timecourse_data.highcoverage_start_1)
initial_idx_2 = sample_list.index(parse_timecourse_data.highcoverage_start_2)
antibiotic_idx = sample_list.index(parse_timecourse_data.highcoverage_antibiotic)
final_idx = sample_list.index(parse_timecourse_data.highcoverage_end)


initial_species_freqs = numpy.fmax(species_freq_matrix[:,initial_idx_1], species_freq_matrix[:,initial_idx_2])


antibiotic_species_freqs = species_freq_matrix[:,antibiotic_idx]
final_species_freqs = species_freq_matrix[:,final_idx]

antibiotic_fold_changes = numpy.clip(antibiotic_species_freqs/(initial_species_freqs+(initial_species_freqs<1e-07)),1.3e-02,1e02/1.3)
final_fold_changes = numpy.clip(final_species_freqs/(initial_species_freqs+(initial_species_freqs<1e-07)),1.3e-02,1e02/1.3)

good_idxs = initial_species_freqs>=1e-04
species = numpy.array(species)[good_idxs]
initial_species_freqs = initial_species_freqs[good_idxs]
antibiotic_fold_changes = antibiotic_fold_changes[good_idxs]
final_fold_changes = final_fold_changes[good_idxs]

species_freq_trajectories = numpy.clip(species_freq_matrix[good_idxs,:],1e-05,1)
species_freq_trajectories = [species_freq_trajectories[idx,:] for idx in xrange(0,species_freq_trajectories.shape[0])]

# sort everything by descending order of XXX
initial_species_freqs, species, antibiotic_fold_changes, final_fold_changes, species_freq_trajectories = (numpy.array(x) for x in zip(*sorted(zip(initial_species_freqs, species, antibiotic_fold_changes, final_fold_changes, species_freq_trajectories), key=lambda pair: (-pair[2],-pair[3]), reverse=True)))

#initial_species_freqs, species, antibiotic_fold_changes, final_fold_changes = (numpy.array(x) for x in zip(*sorted(zip(initial_species_freqs, species, antibiotic_fold_changes, final_fold_changes), key=lambda pair: (-pair[0],-pair[2]), reverse=True)))

pylab.figure(1,figsize=(7,4))
fig = pylab.gcf()
outer_grid  = gridspec.GridSpec(4,1,height_ratios=[1,1,1,2],hspace=0.1)

abundance_axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(abundance_axis)
abundance_axis.set_ylabel('Abundance')

foldchange_axis = plt.Subplot(fig, outer_grid[1])
fig.add_subplot(foldchange_axis)
foldchange_axis.set_ylabel('$f(t)/f(0)$')

ptr_axis = plt.Subplot(fig, outer_grid[2])
fig.add_subplot(ptr_axis)
ptr_axis.set_ylabel('PTR')

fst_axis = plt.Subplot(fig, outer_grid[3])
fig.add_subplot(fst_axis)
ptr_axis.set_ylabel('$F_2$')

xticks = []
xticklabels = []

for species_idx in xrange(0,len(species)):
    species_name = species[species_idx]

    xticks.append(species_idx)
    xticklabels.append(species_name)
    
    # First plot abundances
    f0 = initial_species_freqs[species_idx]
    f_trajectory = species_freq_trajectories[species_idx,:]
    
    hmp_fs = species_abundance_distribution_map[species_name]
        
    abundance_axis.semilogy(numpy.ones_like(f_trajectory)*species_idx, f_trajectory,'ko',markersize=1,zorder=1)
    abundance_axis.semilogy([species_idx],[f0],'bs',markersize=2,zorder=1,markeredgewidth=0)   
    
    good_idxs = (hmp_fs>1e-04)
    if good_idxs.sum()>2:
    
        hmp_fs = hmp_fs[good_idxs]
        log_hmp_fs = numpy.log(hmp_fs)
    
        
        kernel = gaussian_kde(log_hmp_fs)
    
        theory_fs = numpy.logspace(-4,0,100)
        theory_log_fs = numpy.log(theory_fs)
        theory_pdf = kernel(theory_log_fs)
        theory_pdf = theory_pdf / theory_pdf.max() * 0.45
    
        abundance_axis.fill_betweenx(theory_fs,species_idx-theory_pdf, species_idx+theory_pdf,linewidth=0,facecolor='0.7',zorder=0) 
    
    
    ra = antibiotic_fold_changes[species_idx]
    rf = final_fold_changes[species_idx]
    
    if rf > ra:
        symbol = '^'
    else:
        symbol = 'v'
        
    #pylab.plot([species_idx],[ra],'bo',markersize=2)
    foldchange_axis.plot([species_idx,species_idx],[ra,rf],'b-')
    foldchange_axis.plot([species_idx], [rf],symbol,color='b',markersize=2,markeredgewidth=0)

abundance_axis.set_ylim([3e-04,1])
abundance_axis.set_xlim([-2,len(species)+1])
abundance_axis.set_xticklabels([])

foldchange_axis.semilogy([-2,len(species)+1],[1,1],'k:')
foldchange_axis.set_xlim([-2,len(species)+1])
foldchange_axis.set_ylim([1e-02,1e02])
foldchange_axis.fill_between([-2,len(species)+1],[1e-03,1e-03],[2e-02,2e-02],color='0.8')
foldchange_axis.set_xticklabels([])

ptr_axis.set_ylim([0.9,1.9])
ptr_axis.set_xlim([-2,len(species)+1])

fst_axis.set_ylim([1e-06,3e-02])
fst_axis.set_xlim([-2,len(species)+1])
fst_axis.set_xticks(xticks)
fst_axis.set_xticklabels(xticklabels,rotation=90)
fst_axis.tick_params(axis='x', labelsize=5)

fig = pylab.gcf()
fig.savefig('%s/antibiotic_effects.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
#change_fig.savefig('%s/species_freq_change.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
#focal_fig.savefig('%s/species_focal_freq_timecourse.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
