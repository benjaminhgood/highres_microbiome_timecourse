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
#samples = numpy.array([parse_timecourse_data.highcoverage_start_1, parse_timecourse_data.highcoverage_start_2, parse_timecourse_data.highcoverage_antibiotic])
sample_time_map = parse_timecourse_data.parse_sample_time_map()

ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, samples)

samples = samples[sample_idxs]

sys.stderr.write("Loading humann2 genefamily data...\n")
uniref_names, freq_matrix = parse_humann2_data.parse_gene_family_abundance_matrix(samples, include_unmapped=False, include_unknown=False)
sys.stderr.write("Done!\n")

# Calculate depth distribution


pylab.figure(1,figsize=(3.42,2))

abundance_bins = numpy.logspace(-8,-3,300)
abundances = numpy.array(abundance_bins[1:],copy=True)
abundance_bins[0] = -1
abundance_bins[-1] = 1e+08

abundance_threshold = 1e-07

for sample_idx in xrange(0,len(samples)):

    freqs = freq_matrix[:,sample_idx]

    ns, dummy_bins = numpy.histogram(freqs, bins=abundance_bins)
    cdf = numpy.cumsum(ns)

    abundance_cdf = numpy.cumsum(ns*abundances)
    survivals = cdf[-1]-cdf


    pylab.step(abundances, abundance_cdf/abundance_cdf[-1],'-')

pylab.semilogx([1],[1e-03],'k.')
pylab.xlim([1e-08,1e-03])
pylab.xlabel('Gene family abundance, $f$')
pylab.ylabel('% contribution from families $\leq f$') 
fig = pylab.gcf()
fig.savefig(parse_midas_data.analysis_directory+'rare_gene_abundances.pdf',bbox_inches='tight')

common_freq_matrix = freq_matrix*(freq_matrix>=abundance_threshold)
rare_freq_matrix = freq_matrix*(freq_matrix<abundance_threshold)

doubled_freq_matrix = numpy.hstack([common_freq_matrix, rare_freq_matrix])

sys.stderr.write("Calculating KEGG modules...\n")
kegg_modules, doubled_module_freq_matrix = parse_humann2_data.parse_kegg_module_abundance_matrix(samples, level=3, uniref_abundance_data=(uniref_names, doubled_freq_matrix))

common_module_freq_matrix = doubled_module_freq_matrix[:,:len(samples)]
rare_module_freq_matrix = doubled_module_freq_matrix[:,len(samples):]

#print common_module_freq_matrix

#print common_module_freq_matrix.shape
#print rare_module_freq_matrix.shape


# Set up figure
fig2 = plt.figure(figsize=(14, 4))

# Set up grids to hold figure panels
outer_grid = gridspec.GridSpec(1, 2, width_ratios=[4,2], wspace=0.1)

trajectory_grid = gridspec.GridSpecFromSubplotSpec(2, 1, height_ratios=[1,1],
                subplot_spec=outer_grid[0], hspace=0.05) #, hspace=0.08)


common_freq_axis = plt.Subplot(fig2, trajectory_grid[0])
fig2.add_subplot(common_freq_axis)

common_freq_axis.set_ylabel('Common gene families')
#common_freq_axis.set_xlabel('Day')
#freq_axis.set_ylim([1e-04,1])
common_freq_axis.set_ylim([0,1])
common_freq_axis.set_yticks([])

 
common_freq_axis.spines['top'].set_visible(False)
common_freq_axis.spines['right'].set_visible(False)
common_freq_axis.get_xaxis().tick_bottom()
common_freq_axis.get_yaxis().tick_left()

rare_freq_axis = plt.Subplot(fig2, trajectory_grid[1])
fig2.add_subplot(rare_freq_axis)

rare_freq_axis.set_ylabel('Rare gene families')
rare_freq_axis.set_xlabel('Day')
#freq_axis.set_ylim([1e-04,1])
rare_freq_axis.set_ylim([0,1])
rare_freq_axis.set_yticks([])

 
rare_freq_axis.spines['top'].set_visible(False)
rare_freq_axis.spines['right'].set_visible(False)
rare_freq_axis.get_xaxis().tick_bottom()
rare_freq_axis.get_yaxis().tick_left()

legend_axis = plt.Subplot(fig2, outer_grid[1])
fig2.add_subplot(legend_axis)

legend_axis.set_ylim([0,1])
legend_axis.set_xlim([0,1])

legend_axis.spines['top'].set_visible(False)
legend_axis.spines['right'].set_visible(False)
legend_axis.spines['left'].set_visible(False)
legend_axis.spines['bottom'].set_visible(False)

legend_axis.set_xticks([])
legend_axis.set_yticks([])  


######
#
# Plot module abundances themselves
#
#####
    
freq_threshold = 2e-02

common_display_freqs = []
rare_display_freqs = []
display_idxs = []

common_rest_freqs = numpy.zeros_like(common_module_freq_matrix[0,:])
rare_rest_freqs = numpy.zeros_like(rare_module_freq_matrix[0,:])

for i in xrange(0,len(kegg_modules)):
    module_name = kegg_modules[i]    
    common_freqs = common_module_freq_matrix[i,:]
    rare_freqs = rare_module_freq_matrix[i,:]
    if (common_freqs>=freq_threshold).any() or (rare_freqs>=freq_threshold).any():
        # it's a display freq!
        common_display_freqs.append(common_freqs)
        rare_display_freqs.append(rare_freqs)
        display_idxs.append(i)
    
    else:
        common_rest_freqs += common_freqs
        rare_rest_freqs += rare_freqs
        
       
# sort species in descending order of family, then in descending order of abundance        
sorted_idxs = range(0,len(display_idxs))
sorted_idxs = list(sorted(sorted_idxs, key = lambda idx: (common_display_freqs[idx][0], rare_display_freqs[idx][0])))

common_lower = numpy.ones_like(common_rest_freqs)
rare_lower = numpy.ones_like(rare_rest_freqs)

for idx in reversed(sorted_idxs):
    module_name = kegg_modules[display_idxs[idx]]
    
    common_upper = common_lower
    common_lower = common_upper-common_display_freqs[idx]
    
    rare_upper = rare_lower
    rare_lower = rare_upper-rare_display_freqs[idx]
    
    line, = common_freq_axis.plot(ts,-1*numpy.ones_like(ts),'-')
    colorVal = pylab.getp(line,'color')
    legend_axis.plot([-2,-1],[-2,-1],'s',markersize=3,markeredgewidth=0.0, label=module_name)
        
    common_freq_axis.fill_between(ts, common_lower,common_upper,color=colorVal)
    rare_freq_axis.fill_between(ts, rare_lower, rare_upper,color=colorVal)
    
common_upper = common_lower
common_lower = numpy.zeros_like(common_rest_freqs)
common_freq_axis.fill_between(ts, common_lower, common_upper, color='0.7')

rare_upper = rare_lower
rare_lower = numpy.zeros_like(rare_rest_freqs)
rare_freq_axis.fill_between(ts, rare_lower, rare_upper, color='0.7')

legend_axis.plot([-2,-1],[-2,-1],'s',markersize=3,markeredgewidth=0.0, color='0.7',label='Other')
     
legend_axis.legend(loc='upper center',frameon=False,fontsize=5,numpoints=1,ncol=2,handlelength=1)   

common_freq_axis.set_xlim([ts[0],ts[-1]+1])
common_freq_axis.set_xticklabels([])
rare_freq_axis.set_xlim([ts[0],ts[-1]+1])
fig2.savefig('%s/gene_rare_freq_timecourse.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
