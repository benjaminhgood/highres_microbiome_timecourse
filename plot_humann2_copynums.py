import parse_midas_data
import parse_timecourse_data
import parse_humann2_data
import sys
import numpy

desired_samples = numpy.array(parse_timecourse_data.morteza_samples)
sample_time_map = parse_timecourse_data.parse_sample_time_map()

ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, desired_samples)

desired_samples = list(desired_samples[sample_idxs])
start_idx = desired_samples.index(parse_timecourse_data.highcoverage_start_1)
start_idx_2 = desired_samples.index(parse_timecourse_data.highcoverage_start_2)
antibiotic_idx = desired_samples.index(parse_timecourse_data.highcoverage_antibiotic)

desired_species_names = set(['g__Alistipes.s__Alistipes_onderdonkii'])

sys.stderr.write("Loading humann2 genefamily data...\n")
uniref_names, freq_matrix = parse_humann2_data.parse_species_specific_gene_family_abundance_matrix(desired_samples, desired_species_names, include_unknown=False)
sys.stderr.write("Done!\n")

copynum_matrix = freq_matrix/numpy.median(freq_matrix, axis=0)[None,:]
abundances = freq_matrix[:,0]
print "Med = ", numpy.median(abundances)
print "nonzero med =", numpy.median(abundances[abundances>0])
print "Avg = ", abundances.mean()
print "max = ", abundances.max()
print "min = ", abundances.min()
print "nonzero min =", abundances[abundances>0].min()
import pylab


pylab.figure(figsize=(10,2))
for gene_idx in xrange(0,copynum_matrix.shape[0]):

    cs = copynum_matrix[gene_idx,:]
    uniref_name = uniref_names[gene_idx]
    
    if (cs[start_idx:start_idx_2].max()<0.1)*(cs.max()>0.5):
        print uniref_name, cs[start_idx], cs[antibiotic_idx]
        pylab.semilogy(ts,cs,'-',linewidth=0.5,zorder=1)
    else:
        pylab.semilogy(ts,cs,'-',color='0.7',alpha=0.5,zorder=0)

pylab.ylim([1e-02,10])
pylab.savefig('test.png',bbox_inches='tight')
