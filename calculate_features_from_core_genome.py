import sys
import parse_midas_data
import core_gene_utils

#####
#
# Load precalculated fixation events
#
####
target_species=sys.argv[1] 
pangenome_species = parse_midas_data.parse_pangenome_species()

for species_name in pangenome_species:
    
    if not species_name.startswith(target_species):
            continue
            

    core_genes = core_gene_utils.parse_core_genes(species_name,external_filtering=False) # Don't use HMP data in assaying "core"-nes
    shared_pangenome_genes = core_gene_utils.parse_shared_genes(species_name,external_filtering=True) # Use HMP data to weed out potentially shared genes
    core_genes = core_genes - shared_pangenome_genes
    
    for gene_name in sorted(core_genes):
        print species_name, gene_name
        