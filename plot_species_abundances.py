import matplotlib  
matplotlib.use('Agg') 
import parse_midas_data
import parse_timecourse_data
import bacterial_phylogeny_utils
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

#ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, parse_timecourse_data.highcoverage_samples)

species_coverage_matrix = species_coverage_matrix[:,sample_idxs]

total_coverage = species_coverage_matrix.sum(axis=0)
species_freq_matrix = numpy.clip(species_coverage_matrix*1.0/total_coverage,0, 2)    


shannon_diversity = -1*(species_freq_matrix*numpy.log(species_freq_matrix+(species_freq_matrix==0))).sum(axis=0)

js_divergence_matrix = []
# Calculate JSDs
for t0 in xrange(0,len(shannon_diversity)):
    js_divergences = []
    for t in xrange(0,len(shannon_diversity)):
    
        initial_freqs = (species_freq_matrix[:,t0])
        current_freqs = species_freq_matrix[:,t]
        mixture_freqs = initial_freqs*0.5+current_freqs*0.5
    
        good_idxs = (mixture_freqs>0)
    
        initial_freqs = initial_freqs[good_idxs]
        current_freqs = current_freqs[good_idxs]
        mixture_freqs = mixture_freqs[good_idxs]
    
        js_divergence = 0.5*(current_freqs*numpy.log((current_freqs+(current_freqs==0))/mixture_freqs)).sum() + 0.5*(initial_freqs*numpy.log((initial_freqs+(initial_freqs==0))/mixture_freqs)).sum()
        js_divergences.append(js_divergence)
    js_divergence_matrix.append(js_divergences)

js_divergence_matrix = numpy.array(js_divergence_matrix)

# Set up figure
fig = plt.figure(figsize=(14, 2.7))

# Set up grids to hold figure panels
outer_grid = gridspec.GridSpec(1, 2, width_ratios=[4,2], wspace=0.1)

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

freq_axis.set_ylabel('Species')
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




good_species_list = []
prevalences = []
    
freq_threshold = 5e-02

family_freqs = {}
genus_family_map = bacterial_phylogeny_utils.get_genus_family_map()

phylum_freqs = {}
genus_phylum_map = bacterial_phylogeny_utils.get_genus_phylum_map()


display_freqs = []
display_idxs = []

rest_freqs = numpy.zeros_like(species_freq_matrix[0,:])

for i in xrange(0,len(species)):
        
    species_coverages = species_coverage_matrix[i,:]
    species_freqs = species_freq_matrix[i,:]
    species_name = species[i]
    
    family_name = bacterial_phylogeny_utils.get_family_name(species_name, genus_family_map)
    if family_name not in family_freqs:
        family_freqs[family_name] = numpy.zeros_like(rest_freqs)        
    family_freqs[family_name] += species_freqs
    
    phylum_name = bacterial_phylogeny_utils.get_phylum_name(species_name, genus_phylum_map)
    if phylum_name not in phylum_freqs:
        phylum_freqs[phylum_name] = numpy.zeros_like(rest_freqs)        
    phylum_freqs[phylum_name] += species_freqs
    
    
    if (species_freqs>freq_threshold).any():
        # it's a display freq!
        display_freqs.append(species_freqs)
        display_idxs.append(i)
    
        
          
    else:
        rest_freqs += species_freqs

# sort families in ascending order of frequency
sorted_families = list(family_freqs.keys())
sorted_families = list(sorted(sorted_families, key = lambda f: family_freqs[f][0]))
family_order_map = {}
for idx in xrange(0,len(sorted_families)):
    family_order_map[sorted_families[idx]] = idx
    if (family_freqs[sorted_families[idx]]>freq_threshold).any():
        print sorted_families[idx]

# sort families in ascending order of frequency
sorted_phyla = list(phylum_freqs.keys())
sorted_phyla = list(sorted(sorted_phyla, key = lambda f: phylum_freqs[f][0]))
phylum_freq_matrix = numpy.array([phylum_freqs[phylum] for phylum in sorted_phyla])

phylum_order_map = {}
for idx in xrange(0,len(sorted_phyla)):
    phylum_order_map[sorted_phyla[idx]] = idx
    if (phylum_freqs[sorted_phyla[idx]]>freq_threshold).any():
        print sorted_phyla[idx]

#print family_order_map

# sort species in descending order of family, then in descending order of abundance        
sorted_idxs = range(0,len(display_idxs))
#sorted_idxs = list(sorted(sorted_idxs, key = lambda idx: (family_order_map[bacterial_phylogeny_utils.get_family_name(species[idx])], (display_freqs[idx][8])/2.0)))
sorted_idxs = list(sorted(sorted_idxs, key = lambda idx: (family_order_map[bacterial_phylogeny_utils.get_family_name(species[idx])], (display_freqs[idx][0])/2.0)))



lower = numpy.ones_like(rest_freqs)

for idx in reversed(sorted_idxs):
    species_name = species[idx]
    family_name = bacterial_phylogeny_utils.get_family_name(species_name, genus_family_map)
    
    print family_name, species_name
     
    upper = lower
    lower = upper-display_freqs[idx]
    line, = freq_axis.plot(ts,-1*numpy.ones_like(ts),'-')
    colorVal = pylab.getp(line,'color')
    legend_axis.plot([-2,-1],[-2,-1],'s',markersize=3,markeredgewidth=0.0, label=species[idx])
        
    freq_axis.fill_between(ts, lower,upper,color=colorVal)
    
upper = lower
lower = numpy.zeros_like(rest_freqs)
freq_axis.fill_between(ts, lower, upper, color='0.7')
legend_axis.plot([-2,-1],[-2,-1],'s',markersize=3,markeredgewidth=0.0, color='0.7',label='Other')
     
    
alpha_axis.plot(ts, shannon_diversity,'b.-',linewidth=1.5)   
alpha_axis.plot(ts[js_divergence_matrix[0]>0], js_divergence_matrix[0][js_divergence_matrix[0]>0],'r.-',linewidth=1.5)  
alpha_axis.plot(ts[js_divergence_matrix[3]>0], js_divergence_matrix[3][js_divergence_matrix[3]>0],'r.-',linewidth=1.5,alpha=0.5)  


legend_axis.legend(loc='upper center',frameon=False,fontsize=5,numpoints=1,ncol=2,handlelength=1)   



fig.savefig('%s/species_freq_timecourse.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')

### Now do same thing at phylum level

shannon_diversity = -1*(phylum_freq_matrix*numpy.log(phylum_freq_matrix+(phylum_freq_matrix==0))).sum(axis=0)

js_divergence_matrix = []
# Calculate JSDs
for t0 in xrange(0,len(shannon_diversity)):
    js_divergences = []
    for t in xrange(0,len(shannon_diversity)):
    
        initial_freqs = (phylum_freq_matrix[:,t0])
        current_freqs = phylum_freq_matrix[:,t]
        mixture_freqs = initial_freqs*0.5+current_freqs*0.5
    
        good_idxs = (mixture_freqs>0)
    
        initial_freqs = initial_freqs[good_idxs]
        current_freqs = current_freqs[good_idxs]
        mixture_freqs = mixture_freqs[good_idxs]
    
        js_divergence = 0.5*(current_freqs*numpy.log((current_freqs+(current_freqs==0))/mixture_freqs)).sum() + 0.5*(initial_freqs*numpy.log((initial_freqs+(initial_freqs==0))/mixture_freqs)).sum()
        js_divergences.append(js_divergence)
    js_divergence_matrix.append(js_divergences)

js_divergence_matrix = numpy.array(js_divergence_matrix)


# Set up figure
fig = plt.figure(figsize=(14, 2.7))

# Set up grids to hold figure panels
outer_grid = gridspec.GridSpec(1, 2, width_ratios=[4,2], wspace=0.1)

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

freq_axis.set_ylabel('Species')
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


#change_fig.savefig('%s/species_freq_change.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
#focal_fig.savefig('%s/species_focal_freq_timecourse.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')

# sort species in descending order of family, then in descending order of abundance        
sorted_idxs = range(0,len(sorted_phyla))


lower = numpy.ones_like(rest_freqs)
other_freqs = numpy.zeros_like(rest_freqs)

for idx in reversed(sorted_idxs):
    phylum_name = sorted_phyla[idx]
     
    
    if phylum_freq_matrix[idx].max()<1e-02:
        other_freqs+=phylum_freq_matrix[idx]
        continue
     
    upper = lower
    lower = upper-phylum_freq_matrix[idx]
    line, = freq_axis.plot(ts,-1*numpy.ones_like(ts),'-')
    colorVal = pylab.getp(line,'color')
    legend_axis.plot([-2,-1],[-2,-1],'s',markersize=3,markeredgewidth=0.0, label=phylum_name)
        
    freq_axis.fill_between(ts, lower,upper,color=colorVal)
    
upper = lower
lower = numpy.zeros_like(rest_freqs)
freq_axis.fill_between(ts, lower, upper, color='0.7')
legend_axis.plot([-2,-1],[-2,-1],'s',markersize=3,markeredgewidth=0.0, color='0.7',label='Other')
     

legend_axis.legend(loc='upper center',frameon=False,fontsize=5,numpoints=1,ncol=2,handlelength=1)   



fig.savefig('%s/phylum_freq_timecourse.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')

