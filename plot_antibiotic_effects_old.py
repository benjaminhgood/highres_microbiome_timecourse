import matplotlib  
matplotlib.use('Agg') 
import parse_midas_data
import parse_timecourse_data
import bacterial_phylogeny_utils
import pylab
import sys
import numpy
from math import log10, fabs, log

species_coverage_matrix, samples, species = parse_midas_data.parse_global_marker_gene_coverages()
sample_time_map = parse_timecourse_data.parse_sample_time_map()
ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, samples)
samples = numpy.array(samples)[sample_idxs]
species_coverage_matrix = species_coverage_matrix[:,sample_idxs]
total_coverage = species_coverage_matrix.sum(axis=0)
species_freq_matrix = numpy.clip(species_coverage_matrix*1.0/total_coverage,0, 2)    

initial_idx, antibiotic_idx, final_idx = parse_timecourse_data.get_initial_antibiotic_final_idxs(samples)

initial_species_freqs = species_freq_matrix[:,initial_idx]
antibiotic_species_freqs = species_freq_matrix[:,antibiotic_idx]
final_species_freqs = species_freq_matrix[:,final_idx]

for species_idx in xrange(0,len(species)):
    species_name = species[species_idx]
    
    f0 = initial_species_freqs[species_idx]
    fa = antibiotic_species_freqs[species_idx]
    ff = final_species_freqs[species_idx]
    
    if numpy.max([f0,fa,ff]) < 1e-04:
        continue
    
    f0 = max([f0,4e-07])
    fa = max([fa,1e-06])
    ff = max([ff,1e-06])
    
    if fa>f0:
        symbol = '^'
    else:
        symbol = 'v' 
        
    line, = pylab.plot([ff],[f0],'o')
    color = pylab.getp(line,'color')
    pylab.plot([ff,ff],[f0,fa],'-',color=color)
    pylab.plot([ff],[fa],symbol,color=color)

    if species_name.startswith('Alistipes'):
        print species_name, f0,fa,ff

pylab.loglog([1e-05,1],[1e-05,1],'k:')
pylab.xlim([3e-07,1])
pylab.ylim([1e-06,1])
fig = pylab.gcf()
fig.savefig('%s/antibiotic_effects.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
#change_fig.savefig('%s/species_freq_change.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
#focal_fig.savefig('%s/species_focal_freq_timecourse.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
