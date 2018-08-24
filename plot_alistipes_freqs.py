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

mpl.rcParams['font.size'] = 6
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'


pylab.figure(1,figsize=(5,1))


species_coverage_matrix, samples, species = parse_midas_data.parse_global_marker_gene_coverages()
sample_time_map = parse_timecourse_data.parse_sample_time_map()
ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, samples)
samples = numpy.array(samples)[sample_idxs]
species_coverage_matrix = species_coverage_matrix[:,sample_idxs]
total_coverage = species_coverage_matrix.sum(axis=0)
species_freq_matrix = numpy.clip(species_coverage_matrix*1.0/total_coverage,0, 2)    

for species_idx in xrange(0,len(species)):
    
    species_name = species[species_idx]
    
    if not species_name.startswith('Alistipes'):
        continue
     
    short_name = species_name.split("_")[1]
     
    pylab.semilogy(ts, species_freq_matrix[species_idx,:],'.-',markersize=3,label=short_name)

pylab.title('Alistipes species abundances')
pylab.legend(frameon=False,numpoints=1,ncol=5)  
pylab.xlabel('Days')
pylab.ylabel('Abundance')
pylab.ylim([1e-03,1])
fig = pylab.gcf()
fig.savefig('%s/alistipes_freqs.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
#change_fig.savefig('%s/species_freq_change.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
#focal_fig.savefig('%s/species_focal_freq_timecourse.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
