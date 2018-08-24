import numpy
import pylab
import config
import parse_midas_data
import parse_timecourse_data as parse_sample_data
# get all pathways from initial sample
import pathway_utils
import matplotlib.colors as colors

samples, pathways, coverage_matrix = pathway_utils.load_pathways(parse_sample_data.morteza_samples)

coverage_matrix /= coverage_matrix.sum(axis=0)

use_pathways = numpy.zeros(len(pathways))
pathways = numpy.array(pathways)

for sample_idx in xrange(0,len(samples)):
    
    coverages = coverage_matrix[:,sample_idx]
    print len(coverages)
    print len(pathways)
    pathway_idxs = range(0,len(pathways))
    pathway_idxs = sorted(pathway_idxs, key = lambda idx: coverages[idx],reverse=True)
    use_pathways[pathway_idxs[:30]] = 1

pathways = pathways[use_pathways>0.5]
coverage_matrix = coverage_matrix[use_pathways>0.5,:]

for idx in xrange(0,len(pathways)):
    print pathways[idx], coverage_matrix[idx,:]
print len(pathways)
pylab.pcolor(coverage_matrix*100) #,norm=colors.LogNorm())
pylab.xlabel('Samples')
pylab.ylabel('Most abundant pathways')
pylab.colorbar(label='Relative abundance (%)')

pylab.savefig('%s/pathway_abundance_timecourse.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')