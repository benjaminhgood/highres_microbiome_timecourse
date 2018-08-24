import numpy
import sys
from math import log10
from numpy.random import binomial, random_sample
import parse_midas_data
import parse_timecourse_data
import timecourse_utils
import barcode_utils

corrected = False
species_name = 'Bacteroides_vulgatus_57955'
sample_time_map = parse_timecourse_data.parse_sample_time_map()
desired_samples = parse_timecourse_data.morteza_samples

# Load gene coverage information for species_name
sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name,allowed_samples=desired_samples)
sys.stderr.write("Done!\n")  
marker_coverage_times, marker_coverage_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, gene_samples)

marker_coverages = marker_coverages[marker_coverage_idxs]
gene_depth_matrix = gene_depth_matrix[:,marker_coverage_idxs]
gene_reads_matrix = gene_reads_matrix[:,marker_coverage_idxs]
gene_names = numpy.array(gene_names)
gene_copynum_matrix = gene_depth_matrix*1.0/(marker_coverages+(marker_coverages==0))

rare_genes = gene_names[(gene_copynum_matrix.max(axis=1)<0.1)]

sys.stderr.write("%d total genes, %d putatively rare\n" % (len(gene_names), len(rare_genes)))

gene_count_map = {gene_name: 0 for gene_name in rare_genes}

sys.stderr.write("Loading barcodes...\n")
for sample_name in desired_samples:
    sys.stderr.write("%s\n" % sample_name)
    
    # Make sure barcodes exist for this timepoint.
    if not barcode_utils.barcodes_exist(species_name, sample_name):
        continue
         
    # Load barcodes      
    allele_barcode_map = barcode_utils.parse_allele_barcode_tuples(species_name, sample_name,corrected=corrected)

    for gene_name in gene_count_map:
        if gene_name in allele_barcode_map:
            gene_count_map[gene_name] += len(allele_barcode_map[gene_name])
sys.stderr.write("Done!\n")           
for gene_name in sorted(gene_count_map):
    if gene_count_map[gene_name] > 100:
        print gene_name