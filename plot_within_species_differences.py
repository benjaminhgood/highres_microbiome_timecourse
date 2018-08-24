import numpy
import pylab
import diversity_utils
import calculate_substitution_rates
import parse_timecourse_data
import parse_midas_data
import sys

from numpy.random import normal

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
ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, parse_timecourse_data.morteza_samples)
all_samples = numpy.array(parse_timecourse_data.morteza_samples)[sample_idxs]


species_names = []
x_values = []
y_values = [] 
y_errors = []

other_x_values = []
other_y_values = []
other_y_errors = []

third_x_values = []
third_y_values = []
third_y_errors = []

pi_values = []    
  
compressed_fst_data = []  
  
good_species_list = parse_midas_data.parse_good_species_list()
for species_name in good_species_list:
    
    sys.stderr.write("Calculating Fst values for %s...\n" % species_name)           
    substitution_rate_map = calculate_substitution_rates.load_substitution_rate_map(species_name)
    
    desired_samples = numpy.array(parse_timecourse_data.morteza_samples)

    samples, pi_matrix, opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, 'all_pi', allowed_samples=desired_samples)

    samples, f2_difference_matrix, f2_opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, 'all_f2', allowed_samples=desired_samples)

    dummy_samples, f2_variance_matrix, dummy_opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, 'all_f2variance', allowed_samples=desired_samples)

    if len(samples)<2:
        continue
    
    f2_matrix = f2_difference_matrix/f2_opportunity_matrix
    std_f2_matrix = numpy.sqrt(f2_variance_matrix)/f2_opportunity_matrix
    
    initial_idxs = parse_timecourse_data.get_initial_idxs(samples)
    antibiotic_idxs = parse_timecourse_data.get_antibiotic_idxs(samples)
    final_idxs = parse_timecourse_data.get_final_idxs(samples)
    
    if (len(initial_idxs)>1 and len(antibiotic_idxs)>0 and len(final_idxs)>0):
        
        # calculate compressed fst measures
        
        initial_fst_distribution = []
        for i in xrange(0,len(initial_idxs)):
            for j in xrange(i+1,len(initial_idxs)):
                initial_fst_distribution.append( f2_matrix[initial_idxs[i],initial_idxs[j]]) 
                
        max_initial_fst = max(initial_fst_distribution)
                    
        initial_antibiotic_fst_distribution = []
        for i in xrange(0,len(initial_idxs)):
            chunk_initial_antibiotic_fst_distribution = []
            for j in xrange(0,len(antibiotic_idxs)):
                chunk_initial_antibiotic_fst_distribution.append( f2_matrix[initial_idxs[i],antibiotic_idxs[j]])
            initial_antibiotic_fst_distribution.append( max(chunk_initial_antibiotic_fst_distribution) )
        max_initial_antibiotic_fst = max(initial_antibiotic_fst_distribution)
        
        initial_final_fst_distribution = []
        for i in xrange(0,len(initial_idxs)):
            chunk_initial_final_fst_distribution = []
            for j in xrange(0,len(final_idxs)):
                chunk_initial_final_fst_distribution.append( f2_matrix[initial_idxs[i],final_idxs[j]])
            initial_final_fst_distribution.append(max(chunk_initial_final_fst_distribution))
            
        max_initial_final_fst = max(initial_final_fst_distribution) 
        
        antibiotic_final_fst_distribution = []
        for i in xrange(0,len(antibiotic_idxs)):
            chunk_antibiotic_final_fst_distribution = []
            for j in xrange(0,len(final_idxs)):
                chunk_antibiotic_final_fst_distribution.append( f2_matrix[antibiotic_idxs[i],final_idxs[j]])
            antibiotic_final_fst_distribution.append( max(chunk_antibiotic_final_fst_distribution) )
        max_antibiotic_final_fst = max(antibiotic_final_fst_distribution)
        
        initial_fst_distribution = numpy.array(initial_fst_distribution)
        initial_antibiotic_fst_distribution = numpy.array(initial_antibiotic_fst_distribution)
        
        initial_final_fst_distribution = numpy.array(initial_final_fst_distribution)
        antibiotic_final_fst_distribution = numpy.array(antibiotic_final_fst_distribution)
            
        compressed_fst_data.append((species_name, max_initial_fst, max_initial_antibiotic_fst, max_initial_final_fst, max_antibiotic_final_fst, initial_fst_distribution, initial_antibiotic_fst_distribution, initial_final_fst_distribution, antibiotic_final_fst_distribution))
          
    ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, samples)

    initial_idx, antibiotic_idx, final_idx = parse_timecourse_data.get_initial_antibiotic_final_idxs(samples)

    # don't plot if can't compare at the "hallmarks"
    if (initial_idx==-1):
        continue

    # Make sure we have initial, antibiotic, and end    
                 
    pi12_matrix = pi_matrix*1.0/opportunity_matrix
    self_pis = numpy.diag(pi12_matrix)
    pi1_matrix = self_pis[:,None]
    pi2_matrix = self_pis[None,:]
    
    
    pi12s = pi12_matrix[initial_idx]
    f2s = f2_matrix[initial_idx]
    std_f2s = std_f2_matrix[initial_idx]
    
    good_idxs = (ts!=ts[initial_idx])
    
    pi12s = pi12s[good_idxs]
    f2s = f2s[good_idxs]
    std_f2s = std_f2s[good_idxs]
    self_pis = self_pis[good_idxs]
    new_ts = ts[good_idxs]
    
    if (f2s[antibiotic_idx]>(2*std_f2s[antibiotic_idx])): #.any():
        print "*", species_name
    else:
        #continue
        print species_name        
    print "F2(initial->i): ", f2s
    #print f2_matrix[initial_idx,antibiotic_idx], f2_matrix[antibiotic_idx,final_idx]  
    
    species_names.append(species_name)
    x_values.append(new_ts)
    y_values.append(f2s)
    y_errors.append(std_f2s)
    pi_values.append(self_pis)

    f2s = f2_matrix[0]
    std_f2s = std_f2_matrix[0]
    
    good_idxs = (ts!=ts[0])
    
    f2s = f2s[good_idxs]
    std_f2s = std_f2s[good_idxs]
    new_ts = ts[good_idxs]

    other_x_values.append(new_ts)
    other_y_values.append(f2s)
    other_y_errors.append(std_f2s)
    
    f2s = f2_matrix[1]
    std_f2s = std_f2_matrix[1]
    
    good_idxs = (ts!=ts[1])
    
    f2s = f2s[good_idxs]
    std_f2s = std_f2s[good_idxs]
    new_ts = ts[good_idxs]

    third_x_values.append(new_ts)
    third_y_values.append(f2s)
    third_y_errors.append(std_f2s)
    
 
# Calculate antibiotic effects on species abundance
# copy/pasted from plot_antibiotic_effects.py
 
species_coverage_matrix, samples, species = parse_midas_data.parse_global_marker_gene_coverages()
sample_time_map = parse_timecourse_data.parse_sample_time_map()
ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, samples)
samples = numpy.array(samples)[sample_idxs]
species_coverage_matrix = species_coverage_matrix[:,sample_idxs]
total_coverage = species_coverage_matrix.sum(axis=0)
species_freq_matrix = numpy.clip(species_coverage_matrix*1.0/total_coverage,0, 2)    

species_freq_map = {species[species_idx]: species_freq_matrix[species_idx,:] for species_idx in xrange(0,len(species))}

species_freq_matrix = []
for species_idx in xrange(0,len(compressed_fst_data)):
    species_name = compressed_fst_data[species_idx][0]
    
    species_freq_matrix.append( species_freq_map[species_name] )
    
species_freq_matrix = numpy.array(species_freq_matrix)

    
initial_idx, antibiotic_idx, final_idx = parse_timecourse_data.get_initial_antibiotic_final_idxs(samples)

initial_species_freqs = species_freq_matrix[:,initial_idx]
antibiotic_species_freqs = species_freq_matrix[:,antibiotic_idx]
final_species_freqs = species_freq_matrix[:,final_idx]

antibiotic_fold_changes = numpy.clip(antibiotic_species_freqs/(initial_species_freqs+(initial_species_freqs<1e-07)),1.3e-02,1e02/1.3)
final_fold_changes = numpy.clip(final_species_freqs/(initial_species_freqs+(initial_species_freqs<1e-07)),1.3e-02,1e02/1.3)

good_idxs = initial_species_freqs>=1e-05
initial_species_freqs = initial_species_freqs[good_idxs]
antibiotic_fold_changes = antibiotic_fold_changes[good_idxs]
final_fold_changes = final_fold_changes[good_idxs]

# sort everything by descending order of XXX
compressed_fst_data, initial_species_freqs, antibiotic_fold_changes, final_fold_changes = zip(*sorted(zip(compressed_fst_data, initial_species_freqs, antibiotic_fold_changes, final_fold_changes), key=lambda pair: (-pair[2],-pair[3]), reverse=True))

####
#
# Now plot everything
#
###
    
pylab.figure(1,figsize=(3.42,2*len(species_names)))
fig = pylab.gcf()

outer_grid  = gridspec.GridSpec(len(species_names),1, height_ratios=[1]*len(species_names), hspace=0.25)

for species_idx in xrange(0,len(species_names)):

    axis = plt.Subplot(fig, outer_grid[species_idx])
    fig.add_subplot(axis)
    axis.set_title(species_names[species_idx])
    axis.plot(x_values[species_idx], y_values[species_idx],'b.-')
    axis.plot(other_x_values[species_idx], other_y_values[species_idx],'g.-',alpha=0.5)
    axis.plot(third_x_values[species_idx], third_y_values[species_idx],'r.-',alpha=0.5)
    
    axis.set_ylim([0, max([y_values[species_idx].max(),other_y_values[species_idx].max(),third_y_values[species_idx].max()])*1.1])
    
    for data_idx in xrange(0,len(x_values[species_idx])):
        axis.plot([x_values[species_idx][data_idx], x_values[species_idx][data_idx]],[y_values[species_idx][data_idx]-2*y_errors[species_idx][data_idx],y_values[species_idx][data_idx]+2*y_errors[species_idx][data_idx]],'b-')
        
    for data_idx in xrange(0,len(other_x_values[species_idx])):
        axis.plot([other_x_values[species_idx][data_idx], other_x_values[species_idx][data_idx]],[other_y_values[species_idx][data_idx]-2*other_y_errors[species_idx][data_idx],other_y_values[species_idx][data_idx]+2*other_y_errors[species_idx][data_idx]],'g-',alpha=0.5)
    
    for data_idx in xrange(0,len(third_x_values[species_idx])):
        axis.plot([third_x_values[species_idx][data_idx], third_x_values[species_idx][data_idx]],[third_y_values[species_idx][data_idx]-2*third_y_errors[species_idx][data_idx], third_y_values[species_idx][data_idx]+2*third_y_errors[species_idx][data_idx]],'g-',alpha=0.5)
        
fig.savefig('%s/within_species_differences.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')


pylab.figure(2,figsize=(3.42,2*len(species_names)))
fig = pylab.gcf()

outer_grid  = gridspec.GridSpec(len(species_names),1, height_ratios=[1]*len(species_names), hspace=0.25)

for species_idx in xrange(0,len(species_names)):

    axis = plt.Subplot(fig, outer_grid[species_idx])
    fig.add_subplot(axis)
    axis.set_title(species_names[species_idx])
    axis.plot(x_values[species_idx], pi_values[species_idx],'.-')
    
fig.savefig('%s/within_species_diversities.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')

pylab.figure(3,figsize=(12,3))
fig = pylab.gcf()
outer_grid  = gridspec.GridSpec(2,1,height_ratios=[1,1.75],hspace=0.1)

abundance_axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(abundance_axis)
abundance_axis.set_ylabel('$f(t)/f_0$')
f2_axis = plt.Subplot(fig, outer_grid[1])
fig.add_subplot(f2_axis)
f2_axis.set_ylabel('$F_2$')
f2_axis.set_xlabel('Species (arbitrary order)')

xticks = []
xticklabels = []
  
for species_idx in xrange(0,len(compressed_fst_data)):
    
    species_name, control_fst, antibiotic_fst, final_fst, other_fst, initial_fst_distribution, initial_antibiotic_fst_distribution, initial_final_fst_distribution, antibiotic_final_fst_distribution = compressed_fst_data[species_idx]
    
    
    #f2_axis.semilogy([species_idx],[control_fst],'k_',markersize=3)
    #f2_axis.plot([species_idx,species_idx],[antibiotic_fst,final_fst],'k-')
    
    print control_fst, final_fst, antibiotic_fst
    
    #if final_fst > antibiotic_fst*1.1:
    #    f2_axis.plot([species_idx],[final_fst],'k^',markersize=3,markeredgewidth=0)
        #print antibiotic_fst, final_fst,'^',final_fst > antibiotic_fst,final_fst, antibiotic_fst
        #print antibiotic_fst>9e-05, final_fst<9e-05
        #print type(antibiotic_fst), type(final_fst)
    #elif antibiotic_fst>final_fst*1.1:
    #    f2_axis.plot([species_idx],[final_fst],'kv',markersize=3,markeredgewidth=0)   
    #else:
    #    pass
     
    f2_axis.semilogy( [species_idx-0.2]*len(initial_fst_distribution)+normal(0,1,size=len(initial_fst_distribution))*0.03,initial_fst_distribution,'.',color='0.7',markersize=3,alpha=0.5,markeredgewidth=0)
    
    f2_axis.plot([species_idx-0.2], [initial_fst_distribution.mean()],'_',markersize=3,color='0.7')
    
    
    
    f2_axis.plot( [species_idx]*len(initial_antibiotic_fst_distribution)+normal(0,1,size=len(initial_antibiotic_fst_distribution))*0.03,initial_antibiotic_fst_distribution,'.',color='r',markersize=3,alpha=0.5,markeredgewidth=0)
    
    f2_axis.plot([species_idx], [initial_antibiotic_fst_distribution.mean()],'_',markersize=3,color='r')
    
    
    f2_axis.plot( [species_idx+0.2]*len(initial_final_fst_distribution)+normal(0,1,size=len(initial_final_fst_distribution))*0.03,initial_final_fst_distribution,'.',color='b',markersize=3,alpha=0.5,markeredgewidth=0)
    
    f2_axis.plot([species_idx+0.2], [initial_final_fst_distribution.mean()],'_',markersize=3,color='b')
    
        
    ra = antibiotic_fold_changes[species_idx]
    rf = final_fold_changes[species_idx]
    
    if rf > ra:
        symbol = '^'
    else:
        symbol = 'v'
        
    #pylab.plot([species_idx],[ra],'bo',markersize=2)
    abundance_axis.plot([species_idx,species_idx],[ra,rf],'k-')
    abundance_axis.plot( [species_idx],[rf],symbol,color='k',markersize=3,markeredgewidth=0)
    
    xticks.append(species_idx)
    xticklabels.append(species_name)
    
abundance_axis.semilogy([-2,len(compressed_fst_data)+1],[1,1],'k:')

abundance_axis.set_xlim([-2,len(compressed_fst_data)+1])
f2_axis.set_xlim([-2,len(compressed_fst_data)+1])

abundance_axis.set_ylim([1e-02,1e02])

abundance_axis.fill_between([-2,len(compressed_fst_data)+1],[1e-03,1e-03],[2e-02,2e-02],color='0.8')


f2_axis.set_xlabel('Species')
f2_axis.set_xticks(xticks)
f2_axis.set_xticklabels(xticklabels,rotation=90)
abundance_axis.set_xticks([])


fig.savefig('%s/within_species_resiliency.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')

