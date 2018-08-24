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
#species_name='Alistipes_onderdonkii_55464'
#species_name='Barnesiella_intestinihominis_62208'
#species_name='Alistipes_senegalensis_58364'
species_name='Faecalibacterium_prausnitzii_61481'
#species_name='Eubacterium_eligens_61678'
#species_name='Bacteroides_vulgatus_57955'
#species_name='Bacteroides_caccae_53434'
#species_name='Bacteroides_coprocola_61586'

d = cPickle.load( open( pickle_filename,"rb") )
for allele_name in d:
        
    all_barcodes = d[allele_name]['all']
    
    
    if 'species' not in d[allele_name]:
        continue
    
    print "Target gene: %s, %s " % (species_name, allele_name)
    print "Target gene total barcodes: %s" % (", ".join([str(n) for n in all_barcodes]))
    
    other_species_list = d[allele_name]['species'].keys()
    other_species_max_barcodes = []
    for other_species_name in other_species_list:
        max_barcodes = 0
        for other_gene_name in sorted(d[allele_name]['species'][other_species_name]):
                
            if other_gene_name=='any':
                continue
                
            barcodes = d[allele_name]['species'][other_species_name][other_gene_name]
            total_barcodes = barcodes.sum()
            if total_barcodes>max_barcodes:
                max_barcodes=total_barcodes
        
        other_species_max_barcodes.append(max_barcodes)
    
    # sort the species list
    other_species_max_barcodes, other_species_list = (list(x) for x in zip(*sorted(zip(other_species_max_barcodes, other_species_list), key=lambda pair: pair[0], reverse=True)))
    
    for other_species_name in other_species_list[0:2]:
        
        species_barcodes = d[allele_name]['species'][other_species_name]['any']
        
        num_total_genes = 0
        num_good_genes = 0
        
        
            
        for other_gene_name in sorted(d[allele_name]['species'][other_species_name]):
                
            num_total_genes += 1
            if other_gene_name=='any':
                continue
                
            barcodes = d[allele_name]['species'][other_species_name][other_gene_name]
            
            if True:
                num_good_genes += 1
                    
        if num_good_genes>0: 
            
            if other_species_name==species_name:
                print "%s->self:  " % species_name 
            else:
                print "%s->%s: " % (species_name, other_species_name) 
        
                
            
            print "Total barcodes shared w/ species: %s" % (", ".join([str(n) for n in species_barcodes]))
            
            # Found a species with a significant hit.
            # Print out all the genes
                
            for other_gene_name in sorted(d[allele_name]['species'][other_species_name]):
                    
                if other_gene_name=='any':
                    continue
                    
                other_gene_barcodes = d[allele_name]['species'][other_species_name][other_gene_name]
                    
                if other_gene_barcodes.sum()>2.5:

                    #print "%s (%d): %s" % (other_gene_name, other_gene_barcodes.sum(), ", ".join(["%g" % (long((n*1.0/(ntot+(ntot==0))*(ntot>=5))*100)/100.0) for n,ntot in zip(other_gene_barcodes, all_barcodes)])) 
                    print "Linked genes:"   
                    print "%s (%d): %s" % (other_gene_name, other_gene_barcodes.sum(), ", ".join(["%d" % (n) for n,ntot in zip(other_gene_barcodes, all_barcodes)]))    
            
        print "|"
    
print "-"        
print ""
print ""
                