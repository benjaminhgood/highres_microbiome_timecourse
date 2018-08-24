import sys

#####
#
# Load precalculated fixation events
#
####
genes_only = True
snps_only = False
target_species=sys.argv[2] #'Alistipes_senegalensis_58364'
upper_limit = 1e03

file = open(sys.argv[1],"r")
for line in file:
    items = line.split("\t")
    species = items[0]
    change_strs = items[1:]
    snp_changes = []
    gene_changes = []
    
    if target_species!=None:
        if not species.startswith(target_species):
            continue
    
    for change_str in change_strs:
        subitems = change_str[1:-1].split(",")
        if len(subitems)>2:
            # snp change
            gene_name = subitems[0][1:-1]
            contig_name = subitems[1].strip()[1:-1]
            location = long(subitems[2])
            
            snp_changes.append("%s|%d|A" % (contig_name, location))
            snp_changes.append("%s|%d|R" % (contig_name, location))
            
        else:
            gene_name = subitems[0][1:-1]
            gene_changes.append(gene_name)
            
    if not genes_only:
        # can output SNPs
        if len(snp_changes)<upper_limit:
            sys.stderr.write("%s\n" % species)
            for snp_change in snp_changes:
                print species, snp_change
   
    if not snps_only:
        for gene_change in gene_changes:
            print species, gene_change
    