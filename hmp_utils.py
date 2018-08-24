import os
import bz2
import numpy

data_directory = os.path.expanduser("~/ben_nandita_hmp_data/")

def parse_global_marker_gene_coverages():

    file = bz2.BZ2File("%sspecies/coverage.txt.bz2" %  (data_directory),"r")
    line = file.readline() # header
    samples = line.split()[1:]
    species = []
    species_coverage_matrix = []
    for line in file:
        items = line.split()
        species_name = items[0]
        #print items
        coverages = numpy.array([float(item) for item in items[1:]])
        
        species.append(species_name)
        species_coverage_matrix.append(coverages)
    
    file.close()    
    species, species_coverage_matrix = zip(*sorted(zip(species, species_coverage_matrix), key=lambda pair: pair[1].sum(), reverse=True))
    
    species_coverage_matrix = numpy.array(species_coverage_matrix)
    return species_coverage_matrix, samples, species
    
def parse_species_abundance_distributions(fmin=1e-04):
    
    species_coverage_matrix, samples, species = parse_global_marker_gene_coverages()
    
    total_coverage = species_coverage_matrix.sum(axis=0)
    
    good_idxs = (total_coverage>1)
    
    species_coverage_matrix = species_coverage_matrix[:,good_idxs]
    
    freq_matrix = species_coverage_matrix*1.0/species_coverage_matrix.sum(axis=0)[None,:]
    
    species_abundance_distribution_map = {}
    
    for species_idx in xrange(0,len(species)):
        fs = freq_matrix[species_idx]
        species_abundance_distribution_map[species[species_idx]] = fs
        
    return species_abundance_distribution_map
        