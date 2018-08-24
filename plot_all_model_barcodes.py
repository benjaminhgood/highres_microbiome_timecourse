import sys
import cPickle
import numpy

#####
#
# Load precalculated fixation events
#
####

# create output filename
pickle_filename = sys.argv[1]

d = cPickle.load( open( pickle_filename,"rb") )
for feature in d:
    
    focal_species = d[feature]['focal_species'] 
    focal_alleles = [a for a in d[feature]['focal_alleles']]
    allele_str = ", ".join(focal_alleles)
        
    all_barcodes = d[feature]['all']
    
    print "Focal feature: %s, %s" % (focal_species, allele_str)
    print "Total barcodes: %s" % (", ".join([str(n) for n in all_barcodes]))
    
    other_species_list = d[feature]['species'].keys()
    other_species_max_barcodes = []
    for other_species_name in other_species_list:
        max_barcodes = 0
        for other_gene_name in sorted(d[feature]['species'][other_species_name]):
                    
            barcodes, pvalues = d[feature]['species'][other_species_name][other_gene_name]
            total_barcodes = barcodes.sum()
            if total_barcodes>max_barcodes:
                max_barcodes=total_barcodes
        
        other_species_max_barcodes.append(max_barcodes)
    
    # sort the species list
    if len(other_species_list)>0:
        other_species_max_barcodes, other_species_list = (list(x) for x in zip(*sorted(zip(other_species_max_barcodes, other_species_list), key=lambda pair: pair[0], reverse=True)))
    
    # only look at the top two species
    for other_species_name in other_species_list[0:2]:
        if other_species_name==focal_species:
            print "%s->self" % (focal_species)
        else:
            print "%s->%s" % (focal_species, other_species_name)
        print "Linked genes:"
        for other_gene_name in sorted(d[feature]['species'][other_species_name]):
                      
            barcodes, pvalues = d[feature]['species'][other_species_name][other_gene_name]
            
            print "%s (%d): %s" % (other_gene_name, barcodes.sum(), ", ".join(["%d (%g)" % (n,p) for n,ntot,p in zip(barcodes, all_barcodes, pvalues)]))    
            
        print "|"
    
    print "-"        
    print ""
    print ""
                