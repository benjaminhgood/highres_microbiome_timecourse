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
species_name = 'Bacteroides_vulgatus_57955'
print "Rare genes for %s" % species_name
print "-"
print "|"
for allele_name in sorted(d):
        
        
    all_barcodes = d[allele_name]['all']
    
    print "%s: %s" % (allele_name, ", ".join([str(n) for n in all_barcodes]))
    
    if 'species' not in d[allele_name]:
        continue
    
    for other_species_name in d[allele_name]['species']:
        
            species_barcodes = d[allele_name]['species'][other_species_name]['any']
        
            num_total_genes = 0
            num_good_genes = 0
        
            
            for other_gene_name in sorted(d[allele_name]['species'][other_species_name]):
                
                if other_gene_name=='any':
                    continue
                
                num_total_genes += 1
                
                barcodes = d[allele_name]['species'][other_species_name][other_gene_name]
            
                if ((barcodes>=4)*(barcodes>0.049*all_barcodes)).any():
                #if ((all_barcodes>=5)*(barcodes>=6)*(barcodes>=0.5*all_barcodes)).any():
                #if barcodes.sum() >= 10:
                #if (((barcodes>=3).sum()>=2) and (barcodes.sum()>=10)):
                    num_good_genes += 1
                    #print barcodes*1.0/all_barcodes, species_barcodes*1.0/all_barcodes
            
            
                
            if num_good_genes>0: 
            
                if other_species_name==species_name:
                    print "%s->self:  " % species_name, 
                else:
                    print "%s->%s: " % (species_name, other_species_name), 
        
                
            
                print "%s" % (", ".join([str(n) for n in species_barcodes]))
            
                print 
                print "%d, %d" % (num_good_genes, num_total_genes)
            
                
                # Found a species with a significant hit.
                # Print out all the genes
                
                for other_gene_name in sorted(d[allele_name]['species'][other_species_name]):
                    
                    if other_gene_name=='any':
                        continue
                    
                    other_gene_barcodes = d[allele_name]['species'][other_species_name][other_gene_name]
                    
                    if ((other_gene_barcodes>=4)*(other_gene_barcodes>0.049*all_barcodes)).any():

                        print "%s (%d): %s" % (other_gene_name, other_gene_barcodes.sum(), ", ".join(["%g" % (long((n*1.0/(ntot+(ntot==0))*(ntot>=5))*100)/100.0) for n,ntot in zip(other_gene_barcodes, all_barcodes)]))    
        
    print "|"
    
print "-"        
                