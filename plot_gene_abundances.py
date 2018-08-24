import parse_midas_data
import parse_timecourse_data
import numpy
import pylab
import matplotlib 
import parse_humann2_data 
import sys


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


samples = numpy.array(parse_timecourse_data.morteza_samples)
sample_time_map = parse_timecourse_data.parse_sample_time_map()

ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, samples)

samples = samples[sample_idxs]


sys.stderr.write("Loading humann2 genefamily data...\n")
uniref_names, freq_matrix = parse_humann2_data.parse_gene_family_abundance_matrix(parse_timecourse_data.morteza_samples, include_unmapped=False, include_unknown=False)
sys.stderr.write("Done!\n")
kegg_modules, module_freq_matrix = parse_humann2_data.parse_kegg_module_abundance_matrix(samples, level=3, uniref_abundance_data=(uniref_names, freq_matrix))


# Calculate alpha beta diversity for gene families

shannon_diversity = -1*(freq_matrix*numpy.log(freq_matrix+(freq_matrix==0))).sum(axis=0)

js_divergence_matrix = []
# Calculate JSDs
for t0 in xrange(0,len(shannon_diversity)):
    js_divergences = []
    for t in xrange(0,len(shannon_diversity)):
    
        initial_freqs = (freq_matrix[:,t0])
        current_freqs = freq_matrix[:,t]
        mixture_freqs = initial_freqs*0.5+current_freqs*0.5
    
        good_idxs = (mixture_freqs>0)
    
        initial_freqs = initial_freqs[good_idxs]
        current_freqs = current_freqs[good_idxs]
        mixture_freqs = mixture_freqs[good_idxs]
    
        js_divergence = 0.5*(current_freqs*numpy.log((current_freqs+(current_freqs==0))/mixture_freqs)).sum() + 0.5*(initial_freqs*numpy.log((initial_freqs+(initial_freqs==0))/mixture_freqs)).sum()
        js_divergences.append(js_divergence)
    js_divergence_matrix.append(js_divergences)

js_divergence_matrix = numpy.array(js_divergence_matrix)

print js_divergence_matrix[0,:]
print js_divergence_matrix[3,:]

# Do it again for modules

module_shannon_diversity = -1*(module_freq_matrix*numpy.log(module_freq_matrix+(module_freq_matrix==0))).sum(axis=0)

module_js_divergence_matrix = []
# Calculate JSDs
for t0 in xrange(0,len(module_shannon_diversity)):
    js_divergences = []
    for t in xrange(0,len(module_shannon_diversity)):
    
        initial_freqs = (module_freq_matrix[:,t0])
        current_freqs = module_freq_matrix[:,t]
        mixture_freqs = initial_freqs*0.5+current_freqs*0.5
    
        good_idxs = (mixture_freqs>0)
    
        initial_freqs = initial_freqs[good_idxs]
        current_freqs = current_freqs[good_idxs]
        mixture_freqs = mixture_freqs[good_idxs]
    
        js_divergence = 0.5*(current_freqs*numpy.log((current_freqs+(current_freqs==0))/mixture_freqs)).sum() + 0.5*(initial_freqs*numpy.log((initial_freqs+(initial_freqs==0))/mixture_freqs)).sum()
        js_divergences.append(js_divergence)
    module_js_divergence_matrix.append(js_divergences)

module_js_divergence_matrix = numpy.array(module_js_divergence_matrix)

print module_js_divergence_matrix[0,:]
print module_js_divergence_matrix[3,:]

# Set up figure
fig = plt.figure(figsize=(14, 2.7))

# Set up grids to hold figure panels
outer_grid = gridspec.GridSpec(1, 2, width_ratios=[4,2], wspace=0.1)

trajectory_grid = gridspec.GridSpecFromSubplotSpec(3, 1, height_ratios=[0.3,0.3,1],
                subplot_spec=outer_grid[0], hspace=0.05) #, hspace=0.08)

beta_axis = plt.Subplot(fig, trajectory_grid[0])
fig.add_subplot(beta_axis)

beta_axis.set_ylabel('$\\beta$ diversity')
 
beta_axis.spines['top'].set_visible(False)
beta_axis.spines['right'].set_visible(False)
beta_axis.get_xaxis().tick_bottom()
beta_axis.get_yaxis().tick_left()

beta2_axis = plt.Subplot(fig, trajectory_grid[1])
fig.add_subplot(beta2_axis)

beta2_axis.set_ylabel('$\\beta$ diversity')
 
beta2_axis.spines['top'].set_visible(False)
beta2_axis.spines['right'].set_visible(False)
beta2_axis.get_xaxis().tick_bottom()
beta2_axis.get_yaxis().tick_left()

freq_axis = plt.Subplot(fig, trajectory_grid[2])
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
beta_axis.plot([parse_timecourse_data.antibiotic_start, parse_timecourse_data.antibiotic_end],[0.26,0.26],'k-')

beta_axis.plot([parse_timecourse_data.lyme_infection, parse_timecourse_data.antibiotic_start],[0.24,0.24],'k-')

beta_axis.plot([parse_timecourse_data.hrv_infection-1, parse_timecourse_data.hrv_infection+1], [0.24,0.24], 'k-')

beta_axis.text((parse_timecourse_data.lyme_infection+parse_timecourse_data.antibiotic_start)/2.0,0.26,'Lyme',ha='center')
beta_axis.text((parse_timecourse_data.antibiotic_start+parse_timecourse_data.antibiotic_end)/2.0,0.29,'Doxycycline',ha='center')
beta_axis.text(parse_timecourse_data.hrv_infection,0.26,'HRV',ha='center')



beta_axis.set_ylim([0,0.3])
beta_axis.set_yticks([0,0.1,0.2,0.3])

beta_axis.set_xlim([ts[0],ts[-1]+1])
beta_axis.set_xticklabels([])

#beta2_axis.set_ylim([0,0.3])
#beta2_axis.set_yticks([0,0.1,0.2,0.3])

beta2_axis.set_xlim([ts[0],ts[-1]+1])
beta2_axis.set_xticklabels([])


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

# Set up figure 2 (rest pathways)
fig2 = plt.figure(figsize=(14, 2.7))

outer_grid = gridspec.GridSpec(1, 2, width_ratios=[4,2], wspace=0.1)


rest_freq_axis = plt.Subplot(fig2, outer_grid[0])
fig2.add_subplot(rest_freq_axis)

rest_freq_axis.set_ylabel('Genes')
rest_freq_axis.set_xlabel('Day')
#freq_axis.set_ylim([1e-04,1])
rest_freq_axis.set_ylim([0,1])
#rest_freq_axis.set_yticks([])
 
rest_freq_axis.spines['top'].set_visible(False)
rest_freq_axis.spines['right'].set_visible(False)
rest_freq_axis.get_xaxis().tick_bottom()
rest_freq_axis.get_yaxis().tick_left()

rest_legend_axis = plt.Subplot(fig2, outer_grid[1])
fig2.add_subplot(rest_legend_axis)

rest_legend_axis.set_ylim([0,1])
rest_legend_axis.set_xlim([0,1])

rest_legend_axis.spines['top'].set_visible(False)
rest_legend_axis.spines['right'].set_visible(False)
rest_legend_axis.spines['left'].set_visible(False)
rest_legend_axis.spines['bottom'].set_visible(False)

rest_legend_axis.set_xticks([])
rest_legend_axis.set_yticks([])  


#beta_axis.plot(ts, shannon_diversity,'b.-',linewidth=1.5)   
beta_axis.plot(ts[js_divergence_matrix[0]>0], js_divergence_matrix[0][js_divergence_matrix[0]>0],'r.-',linewidth=1.5)  
beta_axis.plot(ts[js_divergence_matrix[3]>0], js_divergence_matrix[3][js_divergence_matrix[3]>0],'r.-',linewidth=1.5,alpha=0.5)  

#beta2_axis.plot(ts, shannon_diversity,'b.-',linewidth=1.5)   
beta2_axis.plot(ts[module_js_divergence_matrix[0]>0], module_js_divergence_matrix[0][module_js_divergence_matrix[0]>0],'r.-',linewidth=1.5)  
beta2_axis.plot(ts[module_js_divergence_matrix[3]>0], module_js_divergence_matrix[3][module_js_divergence_matrix[3]>0],'r.-',linewidth=1.5,alpha=0.5)  




######
#
# Plot module abundances themselves
#
#####
good_species_list = []
prevalences = []
    
freq_threshold = 2e-02
smaller_freq_threshold = 2e-03

display_freqs = []
display_idxs = []

rest_display_freqs = []
rest_display_idxs = []

rest_freqs = numpy.zeros_like(module_freq_matrix[0,:])
rest_rest_freqs = numpy.zeros_like(rest_freqs)

for i in xrange(0,len(kegg_modules)):
        
    module_freqs = module_freq_matrix[i,:]
    module_name = kegg_modules[i]
    if (module_freqs>=freq_threshold).any():
        # it's a display freq!
        display_freqs.append(module_freqs)
        display_idxs.append(i)
    
    else:
        rest_freqs += module_freqs
        
        if (module_freqs>=smaller_freq_threshold).any():
            rest_display_idxs.append(i)
            rest_display_freqs.append(module_freqs)
        else:
            rest_rest_freqs += module_freqs

# sort species in descending order of family, then in descending order of abundance        
sorted_idxs = range(0,len(display_idxs))
sorted_idxs = list(sorted(sorted_idxs, key = lambda idx: ((display_freqs[idx][0])/2.0)))

lower = numpy.ones_like(rest_freqs)

for idx in reversed(sorted_idxs):
    module_name = kegg_modules[display_idxs[idx]]
    print module_name  
    upper = lower
    lower = upper-display_freqs[idx]
    line, = freq_axis.plot(ts,-1*numpy.ones_like(ts),'-')
    colorVal = pylab.getp(line,'color')
    legend_axis.plot([-2,-1],[-2,-1],'s',markersize=3,markeredgewidth=0.0, label=module_name)
        
    freq_axis.fill_between(ts, lower,upper,color=colorVal)
    
upper = lower
lower = numpy.zeros_like(rest_freqs)
freq_axis.fill_between(ts, lower, upper, color='0.7')
legend_axis.plot([-2,-1],[-2,-1],'s',markersize=3,markeredgewidth=0.0, color='0.7',label='Other')
     
legend_axis.legend(loc='upper center',frameon=False,fontsize=5,numpoints=1,ncol=2,handlelength=1)   

freq_axis.set_xlim([ts[0],ts[-1]+1])

####
# Now do same thing for zoomed in on "rest" guys
#
###

# sort species in descending order of family, then in descending order of abundance        
rest_sorted_idxs = range(0,len(rest_display_idxs))
rest_sorted_idxs = list(sorted(rest_sorted_idxs, key = lambda idx: ((rest_display_freqs[idx][0])/2.0)))

lower = rest_freqs

for idx in reversed(rest_sorted_idxs):
    module_name = kegg_modules[rest_display_idxs[idx]]
    print module_name 
    upper = lower
    lower = upper-rest_display_freqs[idx]
    line, = rest_freq_axis.plot(ts,-1*numpy.ones_like(ts),'-')
    colorVal = pylab.getp(line,'color')
    rest_legend_axis.plot([-2,-1],[-2,-1],'s',markersize=3,markeredgewidth=0.0, label=module_name)
        
    rest_freq_axis.fill_between(ts, lower,upper,color=colorVal)
    
upper = lower
lower = numpy.zeros_like(rest_rest_freqs)
rest_freq_axis.fill_between(ts, lower, upper, color='0.7')
rest_legend_axis.plot([-2,-1],[-2,-1],'s',markersize=3,markeredgewidth=0.0, color='0.7',label='Other')
     
rest_legend_axis.legend(loc='upper center',frameon=False,fontsize=5,numpoints=1,ncol=2,handlelength=1)   

rest_freq_axis.set_xlim([ts[0],ts[-1]+1])
rest_freq_axis.set_ylim([0,rest_freqs.max()])

fig.savefig('%s/gene_freq_timecourse.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
fig2.savefig('%s/gene_rest_freq_timecourse.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
