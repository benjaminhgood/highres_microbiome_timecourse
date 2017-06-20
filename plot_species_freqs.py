import matplotlib  
matplotlib.use('Agg') 
import parse_midas_data
import parse_timecourse_data
import pylab
import sys
import numpy
from math import log10, fabs, log

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

species_coverage_matrix, samples, species = parse_midas_data.parse_global_marker_gene_coverages()

sample_time_map = parse_timecourse_data.parse_sample_time_map()

ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, samples)

species_coverage_matrix = species_coverage_matrix[:,sample_idxs]

total_coverage = species_coverage_matrix.sum(axis=0)
species_freq_matrix = species_coverage_matrix*1.0/total_coverage    

shannon_diversity = -1*(species_freq_matrix*numpy.log(species_freq_matrix+(species_freq_matrix==0))).sum(axis=0)


# Set up figure
fig = plt.figure(figsize=(14, 4))

# Set up grids to hold figure panels
outer_grid = gridspec.GridSpec(1, 2, width_ratios=[4,2], wspace=0.1)

trajectory_grid = gridspec.GridSpecFromSubplotSpec(2, 1, height_ratios=[1,1],
                subplot_spec=outer_grid[0], hspace=0.1) #, hspace=0.08)

freq_axis = plt.Subplot(fig, trajectory_grid[0])
fig.add_subplot(freq_axis)

freq_axis.set_ylabel('Species abundance (marker gene)')
freq_axis.set_ylim([1e-04,1])
 
freq_axis.spines['top'].set_visible(False)
freq_axis.spines['right'].set_visible(False)
freq_axis.get_xaxis().tick_bottom()
freq_axis.get_yaxis().tick_left()

alpha_axis = plt.Subplot(fig, trajectory_grid[1])
fig.add_subplot(alpha_axis)

alpha_axis.set_xlabel('Days')
alpha_axis.set_ylabel('Shannon diversity')
 
alpha_axis.spines['top'].set_visible(False)
alpha_axis.spines['right'].set_visible(False)
alpha_axis.get_xaxis().tick_bottom()
alpha_axis.get_yaxis().tick_left()

freq_axis.fill_between([parse_timecourse_data.antibiotic_start, parse_timecourse_data.antibiotic_end],[1e-04,1e-04],[1,1],color=parse_timecourse_data.antibiotics_color,linewidth=0)
freq_axis.fill_between([parse_timecourse_data.lyme_infection, parse_timecourse_data.antibiotic_start],[1e-04,1e-04],[1,1],color=parse_timecourse_data.lyme_color,linewidth=0)
freq_axis.plot([parse_timecourse_data.hrv_infection, parse_timecourse_data.hrv_infection], [1e-04,1], 'k-',linewidth=0.25,zorder=-1)

alpha_axis.fill_between([parse_timecourse_data.antibiotic_start, parse_timecourse_data.antibiotic_end],[2,2],[4,4],color=parse_timecourse_data.antibiotics_color,linewidth=0)
alpha_axis.fill_between([parse_timecourse_data.lyme_infection, parse_timecourse_data.antibiotic_start],[2,2],[4,4],color=parse_timecourse_data.lyme_color,linewidth=0)
alpha_axis.plot([parse_timecourse_data.hrv_infection, parse_timecourse_data.hrv_infection], [2,4], 'k-',linewidth=0.25,zorder=-1)
       

alpha_axis.set_ylim([2,4])

alpha_axis.set_xlim([0,160])
freq_axis.set_xlim([0,160])
freq_axis.set_xticklabels([])


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


# Set up freq change figure
change_fig = plt.figure(figsize=(6, 2.5))

special_t_idx = 11

# Set up grids to hold figure panels
change_outer_grid = gridspec.GridSpec(1, 2, width_ratios=[2,1], wspace=0.1)

change_axis = plt.Subplot(change_fig, change_outer_grid[0])
change_fig.add_subplot(change_axis)

change_axis.set_xlabel('Abundance (day %d)' % ts[0])
change_axis.set_ylabel('Abundance (day %d)' % ts[special_t_idx])

change_axis.set_ylim([1e-04,1])
change_axis.set_xlim([1e-04,1])


change_axis.loglog([1e-04,1],[1e-04,1],'-',color='0.7')
 
change_axis.spines['top'].set_visible(False)
change_axis.spines['right'].set_visible(False)
change_axis.get_xaxis().tick_bottom()
change_axis.get_yaxis().tick_left()

change_legend_axis = plt.Subplot(change_fig, change_outer_grid[1])
change_fig.add_subplot(change_legend_axis)

change_legend_axis.set_ylim([0,1])
change_legend_axis.set_xlim([0,1])

change_legend_axis.spines['top'].set_visible(False)
change_legend_axis.spines['right'].set_visible(False)
change_legend_axis.spines['left'].set_visible(False)
change_legend_axis.spines['bottom'].set_visible(False)

change_legend_axis.set_xticks([])
change_legend_axis.set_yticks([]) 

# Set up figure
focal_fig = plt.figure(figsize=(14, 2.3))

# Set up grids to hold figure panels
focal_outer_grid = gridspec.GridSpec(1, 2, width_ratios=[4,2], wspace=0.1)



focal_freq_axis = plt.Subplot(focal_fig, focal_outer_grid[0])
focal_fig.add_subplot(focal_freq_axis)

focal_freq_axis.set_ylabel('Species abundance (marker gene)')
focal_freq_axis.set_xlabel('Days')

focal_freq_axis.set_ylim([1e-04,1])
focal_freq_axis.set_xlim([0,160])
 
focal_freq_axis.spines['top'].set_visible(False)
focal_freq_axis.spines['right'].set_visible(False)
focal_freq_axis.get_xaxis().tick_bottom()
focal_freq_axis.get_yaxis().tick_left()

focal_legend_axis = plt.Subplot(focal_fig, focal_outer_grid[1])
focal_fig.add_subplot(focal_legend_axis)

focal_legend_axis.set_ylim([0,1])
focal_legend_axis.set_xlim([0,1])

focal_legend_axis.spines['top'].set_visible(False)
focal_legend_axis.spines['right'].set_visible(False)
focal_legend_axis.spines['left'].set_visible(False)
focal_legend_axis.spines['bottom'].set_visible(False)

focal_legend_axis.set_xticks([])
focal_legend_axis.set_yticks([])  

focal_freq_axis.fill_between([parse_timecourse_data.antibiotic_start, parse_timecourse_data.antibiotic_end],[1e-04,1e-04],[1,1],color=parse_timecourse_data.antibiotics_color,linewidth=0)
focal_freq_axis.fill_between([parse_timecourse_data.lyme_infection, parse_timecourse_data.antibiotic_start],[1e-04,1e-04],[1,1],color=parse_timecourse_data.lyme_color,linewidth=0)
focal_freq_axis.plot([parse_timecourse_data.hrv_infection, parse_timecourse_data.hrv_infection], [1e-04,1], 'k-',linewidth=0.25,zorder=-1)


good_species_list = []
prevalences = []
    

for i in xrange(0,len(species)):
        
    species_coverages = species_coverage_matrix[i,:]
    species_freqs = species_freq_matrix[i,:]
    clipped_species_freqs = numpy.clip(species_freqs,1e-04,1)
    prevalence = (species_coverages>=min_marker_coverage).sum()
    
    if (species_coverages>=20).sum() >= 2:
        
        if species[i].startswith('Bacteroides'):
            linewidth=1.5
            alpha=1
            zorder=4
            symbol = 's'
        else:
            linewidth=0.5
            alpha=0.5
            zorder=1
            symbol = 'o'
        
        line, = freq_axis.semilogy(ts, species_freqs,'.-',markersize=3,linewidth=linewidth,alpha=alpha,zorder=zorder)
        colorVal = pylab.getp(line,'color')
        legend_axis.plot([-2,-1],[-2,-1],'.-',markersize=3,markeredgewidth=0.0, label=species[i],alpha=alpha,linewidth=linewidth)
        
        
        # if underwent factor of 10 change
        if fabs(log10(clipped_species_freqs[0]/clipped_species_freqs[special_t_idx]))>0.9:
            
            line, = change_axis.plot([clipped_species_freqs[0]],[clipped_species_freqs[special_t_idx]],symbol,markersize=4,markeredgewidth=0,alpha=alpha)
            colorVal = pylab.getp(line,'color')
            change_legend_axis.plot([-2,-1],[-2,-1],symbol,markersize=3,markeredgewidth=0.0, label=species[i],alpha=alpha,linewidth=linewidth)
            
            focal_freq_axis.semilogy(ts, clipped_species_freqs,'.-',markersize=3,linewidth=linewidth,alpha=alpha,zorder=zorder)
            focal_legend_axis.plot([-2,-1],[-2,-1],'.-',markersize=3,markeredgewidth=0.0, label=species[i],alpha=alpha,linewidth=linewidth)
            
        else:
            change_axis.plot([clipped_species_freqs[0]],[clipped_species_freqs[special_t_idx]],symbol,markersize=2,color='0.7',markeredgewidth=0)
            
        
    if prevalence >= min_prevalence:
        good_species_list.append(species[i])
        prevalences.append(prevalence)
        
    

alpha_axis.plot(ts, shannon_diversity,'b.-',linewidth=1.5)        
legend_axis.legend(loc='upper center',frameon=False,fontsize=5,numpoints=1,ncol=2,handlelength=1)   
change_legend_axis.legend(loc='upper center',frameon=False,fontsize=5,numpoints=1,ncol=1,handlelength=1)   
focal_legend_axis.legend(loc='upper center',frameon=False,fontsize=5,numpoints=1,ncol=1,handlelength=1)   



fig.savefig('%s/species_freq_timecourse.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
change_fig.savefig('%s/species_freq_change.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
focal_fig.savefig('%s/species_focal_freq_timecourse.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
