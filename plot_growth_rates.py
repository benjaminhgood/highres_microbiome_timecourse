import pylab
import numpy
import sys
import parse_midas_data
import parse_timecourse_data
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl

from math import log10,floor,ceil

species_coverage_matrix, samples, midas_species_names = parse_midas_data.parse_global_marker_gene_coverages()

sample_time_map = parse_timecourse_data.parse_sample_time_map()

ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, samples)

species_coverage_matrix = species_coverage_matrix[:,sample_idxs]

total_coverage = species_coverage_matrix.sum(axis=0)
species_freq_matrix = numpy.clip(species_coverage_matrix*1.0/total_coverage,0, 2)    
species_freq_ts = ts

midas_simple_species_names = []
for species_name in midas_species_names:
    
    name_items = species_name.split("_")
    midas_simple_species_names.append("%s %s" % (name_items[0], name_items[1]))


file = open("ptrs5.csv","r")
header = file.readline()
sample_strs = header.split(",")

samples = []
for sample_str in sample_strs:
    if sample_str.strip()=="":
        continue
        
    sample_substrs = sample_str.split("_")
    
    samples.append(sample_substrs[0])
    
sample_time_map = parse_timecourse_data.parse_sample_time_map()

sample_times = numpy.array([sample_time_map[sample] for sample in samples])

data = {}

for line in file:
    items = line.split(",")
    species_name = items[0].strip()
    
    ptr_strs = items[1:]
    #print species_name, ptr_strs
    ptrs = []
    
    for ptr_str in ptr_strs:
        
        if ptr_str.strip()=="":
            ptrs.append(-1)
        else:
            ptrs.append(float(ptr_str))
    
    ptrs = numpy.array(ptrs)
    
    new_sample_times, ptrs = (numpy.array(t) for t in zip(*sorted(zip(sample_times, ptrs))))
    
    data[species_name] = (new_sample_times[ptrs>0], ptrs[ptrs>0])

ptr_species_names = set(data.keys())

good_data = {}
for idx in xrange(0,len(midas_species_names)):
    
    species_name = midas_species_names[idx]
    simple_species_name = midas_simple_species_names[idx]
    
    if simple_species_name in data:
        
        good_data[species_name] = (species_freq_ts, species_freq_matrix[idx,:], data[simple_species_name][0], data[simple_species_name][1])



species_list = good_data.keys()
num_axes = len(species_list)

antibiotic_t = sample_time_map[parse_timecourse_data.highcoverage_antibiotic]

mpl.rcParams['font.size'] = 6.0
mpl.rcParams['lines.linewidth'] = 0.25
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

fig = plt.figure(figsize=(7, num_axes*1.7))

# Set up grids to hold figure panels
outer_grid = gridspec.GridSpec(num_axes, 2, wspace=0.3, hspace=0.1)

for idx in xrange(0,num_axes):
    
    
    species_name = species_list[idx]
    
    name_items = species_name.split("_")

    pretty_name = "%s %s" % (name_items[0], name_items[1])
    
    freq_ts, freqs, ptr_ts, ptrs = good_data[species_name]
     
    ptr_axis = plt.Subplot(fig, outer_grid[idx,0])
    fig.add_subplot(ptr_axis)
   
    ptr_axis.plot(ptr_ts, ptrs, 'b.-')
    ptr_axis.plot([antibiotic_t, antibiotic_t], [0,2],'k:')   
    ptr_axis.set_ylabel(pretty_name)
    ptr_axis.set_ylim([0.8,1.6])
    ptr_axis.set_xlim([0,160])
    
    freq_axis = plt.Subplot(fig, outer_grid[idx, 1])
    fig.add_subplot(freq_axis)
   
    max_freq = 10**(ceil(log10(freqs.max())))
    min_freq = 10**(floor(log10(freqs[freqs>0].min())))
    
   
    freq_axis.semilogy(freq_ts, freqs, 'r.-')
    freq_axis.plot([antibiotic_t, antibiotic_t], [1e-06,2],'k:')   
    freq_axis.set_ylim([min_freq, max_freq])
    #freq_axis.set_ylim([0,freqs.max()*1.1])
    
    freq_axis.set_xlim([0,160])
    
    
    if idx==0:
        ptr_axis.set_title('PTR')
        freq_axis.set_title('Relative abundance')
        
    if idx==(len(species_list)-1):
        ptr_axis.set_xlabel('Day')
        freq_axis.set_xlabel('Day')
    else:
        freq_axis.set_xticklabels([])
        ptr_axis.set_xticklabels([])
    
        


pylab.savefig('ptr.pdf',bbox_inches='tight')
            
