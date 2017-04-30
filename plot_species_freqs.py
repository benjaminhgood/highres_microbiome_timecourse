import matplotlib  
matplotlib.use('Agg') 
import parse_midas_data
import parse_timecourse_data
import pylab
import sys
import numpy

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

min_marker_coverage = 20
min_prevalence=5

bacteroides_color = '#084594'
alistipes_color = '#B10026'
rest_color = '0.7'


# Set up figure
fig = plt.figure(figsize=(14, 2.5))

# Set up grids to hold figure panels
outer_grid = gridspec.GridSpec(1, 2, width_ratios=[4,2], wspace=0.1)

main_axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(main_axis)

main_axis.set_xlabel('Days')
main_axis.set_ylabel('Frequency of species')
main_axis.set_ylim([1e-04,1])
 
main_axis.spines['top'].set_visible(False)
main_axis.spines['right'].set_visible(False)
main_axis.get_xaxis().tick_bottom()
main_axis.get_yaxis().tick_left()


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


good_species_list = []
prevalences = []
    
species_coverage_matrix, samples, species = parse_midas_data.parse_global_marker_gene_coverages()

sample_time_map = parse_timecourse_data.parse_sample_time_map()

ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, samples)

for i in xrange(0,len(species)):
        
    species_coverages = species_coverage_matrix[i,:]
    total_coverage = species_coverage_matrix.sum(axis=0)
    species_freqs = species_coverages*1.0/total_coverage
    prevalence = (species_coverages>=min_marker_coverage).sum()
    
    if (species_coverages[sample_idxs]>=20).sum() >= 2:
        
        if species[i].startswith('Bacteroides'):
            linewidth=1.5
            alpha=1
            zorder=1
        else:
            linewidth=0.5
            alpha=0.5
            zorder=0
        
        line, = main_axis.semilogy(ts, species_freqs[sample_idxs],'.-',markersize=3,linewidth=linewidth,alpha=alpha,zorder=zorder)
        colorVal = pylab.getp(line,'color')
        legend_axis.plot([-2,-1],[-2,-1],'.-',markersize=3,markeredgewidth=0.0, label=species[i],alpha=alpha,linewidth=linewidth)
    
        
    if prevalence >= min_prevalence:
        good_species_list.append(species[i])
        prevalences.append(prevalence)
        
legend_axis.legend(loc='upper center',frameon=False,fontsize=5,numpoints=1,ncol=2,handlelength=1)   
main_axis.set_xlim([0,ts.max()+1])
fig.savefig('%s/species_freq_timecourse.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')