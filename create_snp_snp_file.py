import parse_midas_data
import parse_timecourse_data
import timecourse_utils
import core_gene_utils
import sys

################################################################################
#
# Standard header to read in argument information
#
################################################################################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('species', type=str, help="species")
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
args = parser.parse_args()

debug = args.debug
chunk_size = args.chunk_size
species_prefix = args.species
################################################################################


min_freq = 0.2
min_coverage = 10
desired_samples = parse_timecourse_data.morteza_samples

pangenome_species = parse_midas_data.parse_pangenome_species()
for species_name in pangenome_species:
    if species_name.startswith(species_prefix):
        break

sample_time_map = parse_timecourse_data.parse_sample_time_map()        
sys.stderr.write("Loading core genes...\n")
core_genes = core_gene_utils.parse_core_genes(species_name,external_filtering=False) # Don't use HMP data in assaying "core"-nes)
non_shared_genes = core_gene_utils.parse_non_shared_reference_genes(species_name) # Use HMP filter to remove shared genes
core_genes = core_genes & non_shared_genes
sys.stderr.write("Done! %d core genes \n" % (len(core_genes)))

print ">%s" % species_name
num_tracked=0
final_line_number = 0
while final_line_number >= 0:
    
    sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
    samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_variant_types=set(['1D','2D','3D','4D']),chunk_size=chunk_size,allowed_samples=desired_samples, initial_line_number=final_line_number, allowed_genes=core_genes)
    sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
    
    sample_ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, samples)

    # Calculate fixation matrix
    sys.stderr.write("Calculating allele freqs...\n")
    chunk_alts, chunk_depths, chunk_snp_infos = timecourse_utils.calculate_read_count_matrix(allele_counts_map, passed_sites_map, allowed_variant_types=set(['1D','2D','3D','4D']))    
    sys.stderr.write("Done!\n")
    
    chunk_alts = chunk_alts[:,sample_idxs]
    chunk_depths = chunk_depths[:,sample_idxs]
    
    for snp_idx in xrange(0,len(chunk_alts)):
        
        chromosome, location, gene_name, variant_type = chunk_snp_infos[snp_idx]
        
        if gene_name not in core_genes:
            continue
        
        good_idxs = chunk_depths[snp_idx,:]>min_coverage
        
        if (good_idxs).sum() < 10:
            continue
        
        As = chunk_alts[snp_idx, good_idxs]  
        Ds = chunk_depths[snp_idx,good_idxs]
        fs = As*1.0/Ds
        
        avg_f = fs.mean()
        
        if avg_f>0.5:
            # polarize
            fs = 1-fs
            
        if (fs<min_freq).all():
            continue
        
        num_tracked+=1    
        # A good SNV to track!
        print "%s|%d" % (chromosome, location)  
        
sys.stderr.write("Done! Tracking %d SNVs\n" % num_tracked) 