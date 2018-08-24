import numpy
import sys

snp_output_file = open(sys.argv[1],"r")
desired_panel_idx = long(sys.argv[2])

feature_species = None
feature_alleles = []
snp_output_file.readline()
for line in snp_output_file:
    items = line.split()
    species = items[0]
    allele = items[1]
    panel_idx = long(items[2])
    
    if panel_idx==desired_panel_idx:
        feature_species = species
        feature_alleles.append(allele)
        
feature_antialleles = []
for allele in feature_alleles:
    
    items = allele.split("|")
    contig = items[0]
    location = items[1]
    polarization = items[2]
    if polarization=='A':
        anti_polarization='R'
    else:
        anti_polarization='A'
        
    anti_allele = "|".join((contig,location,anti_polarization))
    feature_antialleles.append(anti_allele)
    
print "\t".join([feature_species]+feature_alleles)
print "\t".join([feature_species]+feature_antialleles)

    
    

