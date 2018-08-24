import config
import numpy

def load_pathways(samples):

    pathway_coverages = {}
    
    for sample_idx in xrange(0,len(samples)):
        sample_name = samples[sample_idx]
        file = open("%s%s.join_pathabundance.tsv" % (config.humann2_directory, sample_name),"r")

        pathways = []
        coverages = []

        file.readline()
        for line in file:
            items = line.split("\t") 
            pathway_name = items[0]
            pathway_rpk = float(items[1])
            if "|" in pathway_name:
                continue # species resolved version
            else:
                if pathway_name not in pathway_coverages:
                    pathway_coverages[pathway_name] = numpy.zeros(len(samples))*1.0
                pathway_coverages[pathway_name][sample_idx] = pathway_rpk
                
        file.close()
        
        
    pathways = list(sorted(pathway_coverages))
    coverage_matrix = []
    final_pathways = []
    for pathway_name in pathways:
        if pathway_name == 'UNMAPPED':
            continue
        if pathway_name == 'UNINTEGRATED':
            continue
        final_pathways.append(pathway_name)
        coverage_matrix.append(pathway_coverages[pathway_name])
        
    coverage_matrix = numpy.array(coverage_matrix)
    final_pathways = numpy.array(final_pathways)
    
    return samples, final_pathways, coverage_matrix
    