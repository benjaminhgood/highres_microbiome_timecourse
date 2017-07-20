import parse_timecourse_data
import barcode_utils

species_name = "Bacteroides_vulgatus_57955"

################################################################################
#
# Standard header to read in argument information
#
################################################################################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
parser.add_argument('--other-species', type=str, help='Run the script for a different species')

args = parser.parse_args()

debug = args.debug
chunk_size = args.chunk_size
other_species = args.other_species

if other_species:
    species_name = other_species
    other_species_str = "_%s" % species_name
else:
    other_species_str = ""
 

target_snp_groups = {}

 
sample_name =  parse_timecourse_data.highcoverage_end


target_snp_filename = "Bacteroides_vulgatus_57955_final_snps.txt"
target_snp_file = open(target_snp_filename,"r")
target_snp_file.readline() # header
for line in target_snp_file:
    items = line.split("\t")
    if items[0].strip() == species_name:
        
        contig = items[1].strip()
        position = long(items[2])
        allele = items[3].strip()
        group = long(items[4])
        
        
        
        desired_alleles.append(gene_name)
        
allele_barcode_map = barcode_utils.parse_allele_barcodes(species_name, sample_name)
barcode_allele_map = barcode_utils.calculate_barcode_allele_map(allele_barcode_map)
num_shared_barcodes = barcode_utils.calculate_num_shared_barcodes(allele_barcode_map, barcode_allele_map, desired_alleles)

for allele in desired_alleles:
    linked_genes = []
    linked_snps = {}
    
    for other_allele, n in sorted(num_shared_barcodes[allele].iteritems(), key=lambda (k,v): (v,k), reverse=True):
        
        if other_allele[-1] in set(['A','R']):
            # A snp
            short_other_allele = other_allele[:-2]
            if short_other_allele not in linked_snps:
                linked_snps[short_other_allele] = [0, 0]
            
            if other_allele[-1]=='A':
                linked_snps[short_other_allele][0] += n
                linked_snps[short_other_allele][1] += n
            else:
                linked_snps[short_other_allele][0] += n
        
        else:
           # A gene
           if n>1:
               linked_genes.append((other_allele, n))
    
           
    print allele
    print ", ".join(["%s:%d" % (other_allele, n) for other_allele, n in linked_genes])
    snp_strs = []
    for snp_allele, allele_counts in sorted(linked_snps.iteritems(), key=lambda (k,v): (v[0],k), reverse=True):
        if allele_counts[0]>1:
            snp_strs.append("%s:%d|%d" % (snp_allele, allele_counts[0], allele_counts[1]))
    print ", ".join(snp_strs)
    
    