import sys
import pylab
import numpy
import parse_midas_data
import parse_timecourse_data
import stats_utils
  

sample_time_map = parse_timecourse_data.parse_sample_time_map()
    
species_coverage_matrix, samples, species = parse_midas_data.parse_global_marker_gene_coverages()

species_times, species_time_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, samples)
desired_samples = numpy.array(samples)[species_time_idxs]

species_coverage_matrix = species_coverage_matrix[:,species_time_idxs]

coverage_to_abundance_factor = 1.0/species_coverage_matrix.sum(axis=0)

species_abundance_matrix = species_coverage_matrix * coverage_to_abundance_factor[None,:] 


####
#
# Now do same thing for genes in pangenome
#
####

for species_name in ['new_species']:
        
    sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
    gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix =     parse_midas_data.parse_pangenome_data(species_name)
    sys.stderr.write("Done!\n")

    species_times, species_time_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, gene_samples)
    
    if len(species_times)==0:
        continue

    desired_gene_samples = numpy.array(gene_samples)[species_time_idxs]

    gene_depth_matrix = gene_depth_matrix[:,species_time_idxs]

    gene_sample_idx_map = parse_midas_data.calculate_sample_idx_map(desired_gene_samples, desired_samples)
    
    gene_coverage_to_abundance_factor = numpy.array([coverage_to_abundance_factor[gene_sample_idx_map[i]] for i in xrange(0,gene_depth_matrix.shape[1])])
        
    gene_abundances = gene_depth_matrix * gene_coverage_to_abundance_factor

    sys.stderr.write(", ".join(desired_gene_samples)+"\n")

    sys.stderr.write(", ".join([str(1.0/c) for c in gene_coverage_to_abundance_factor])+"\n")

    initial_idxs = numpy.logical_or(desired_gene_samples==parse_timecourse_data.highcoverage_start_1, desired_gene_samples==parse_timecourse_data.highcoverage_start_2)
    hrv_idxs = desired_gene_samples==parse_timecourse_data.highcoverage_hrv
    lyme_idxs = desired_gene_samples==parse_timecourse_data.highcoverage_lyme
    antibiotic_idxs = desired_gene_samples==parse_timecourse_data.highcoverage_antibiotic
    final_idxs = desired_gene_samples==parse_timecourse_data.highcoverage_end
    
    initial_abundances = gene_abundances[:,initial_idxs].mean(axis=1)
    hrv_abundances = gene_abundances[:,hrv_idxs].mean(axis=1)
    lyme_abundances = gene_abundances[:,lyme_idxs].mean(axis=1)
    antibiotic_abundances = gene_abundances[:,antibiotic_idxs].mean(axis=1)
    final_abundances = gene_abundances[:,final_idxs].mean(axis=1)
    
    hrv_gene_idxs = numpy.logical_and(initial_abundances < 0.1*hrv_abundances, hrv_abundances>1e-03)
     
    lyme_gene_idxs = numpy.logical_and(initial_abundances < 0.1*lyme_abundances, lyme_abundances>1e-03)
    
    antibiotic_gene_idxs = numpy.logical_and(initial_abundances < 0.1*antibiotic_abundances, antibiotic_abundances>1e-03)
    
    final_gene_idxs = numpy.logical_and(initial_abundances < 0.1*final_abundances, final_abundances>1e-03)
    
    #any_gene_idxs = final_gene_idxs
    any_gene_idxs = antibiotic_gene_idxs
    #any_gene_idxs = (hrv_gene_idxs+lyme_gene_idxs+antibiotic_gene_idxs+final_gene_idxs)>0
    
    sys.stderr.write("HRV genes: %d\n" % hrv_gene_idxs.sum())
    sys.stderr.write("Lyme genes: %d\n" % lyme_gene_idxs.sum())
    sys.stderr.write("Antibiotic genes: %d\n" % antibiotic_gene_idxs.sum())
    sys.stderr.write("Final genes: %d\n" % final_gene_idxs.sum())
    sys.stderr.write("Any genes: %d\n" % any_gene_idxs.sum())
    
    desired_gene_names = numpy.array(gene_names)[any_gene_idxs]
    
    for gene_name in desired_gene_names:
        print gene_name