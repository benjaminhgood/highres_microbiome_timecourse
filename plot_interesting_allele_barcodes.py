import sys
import cPickle
import numpy

#####
#
# Load precalculated fixation events
#
####

# create output filename
within_host_filename = sys.argv[1]
pickle_filename = sys.argv[2]

# First load interesting alleles
file = open(within_host_filename,"r")
interesting_alleles = {}
for line in file:
    items = line.split("\t")
    species = items[0]
    change_strs = items[1:]
    
    for change_str in change_strs:
        subitems = change_str[1:-1].split(",")
        if len(subitems)>2:
            # snp change
            gene_name = subitems[0][1:-1]
            contig_name = subitems[1].strip()[1:-1]
            location = long(subitems[2])
            
            if species not in interesting_alleles:
                interesting_alleles[species] = []
            
            interesting_alleles[species].append("%s|%d|A" % (contig_name, location))
            interesting_alleles[species].append("%s|%d|R" % (contig_name, location))
            
        else:
            gene_name = subitems[0][1:-1]
            if species not in interesting_alleles:
                interesting_alleles[species] = [] 
            interesting_alleles[species].append(gene_name)
file.close()


d = cPickle.load( open( pickle_filename,"rb") )
for species_name in interesting_alleles.keys():
    print "Target species %s:" % species_name
    print "-"
    print "|"
    for allele_name in interesting_alleles[species_name]:
        
        if allele_name not in d:
            continue
        
        all_barcodes = d[allele_name]['all']
    
        print "Target gene %s: %s" % (allele_name, ", ".join([str(n) for n in all_barcodes]))
    
        if 'species' not in d[allele_name]:
            continue
    
        for other_species_name in d[allele_name]['species']:
        
            species_barcodes = d[allele_name]['species'][other_species_name]['any']
        
            num_total_genes = 0
            num_good_genes = 0
        
            
            for other_gene_name in sorted(d[allele_name]['species'][other_species_name]):
                
                num_total_genes += 1
                if other_gene_name=='any':
                    continue
                
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
        
                
            
                print "# shared barcodes at species level: %s" % (", ".join([str(n) for n in species_barcodes]))
            
                #print "%d, %d" % (num_good_genes, num_total_genes)
                print "Linked genes:"
                
                # Found a species with a significant hit.
                # Print out all the genes
                
                for other_gene_name in sorted(d[allele_name]['species'][other_species_name]):
                    
                    if other_gene_name=='any':
                        continue
                    
                    other_gene_barcodes = d[allele_name]['species'][other_species_name][other_gene_name]
                    
                    if (other_gene_barcodes.sum()>=10) and ((other_gene_barcodes>=4).sum()>=2) and ((other_gene_barcodes>0.049*all_barcodes).sum()>=1):

                        print "%s (%d): %s" % (other_gene_name, other_gene_barcodes.sum(), ", ".join(["%g" % (long((n*1.0/(ntot+(ntot==0))*(ntot>=5))*100)/100.0) for n,ntot in zip(other_gene_barcodes, all_barcodes)]))    
            
        print "|"
    
    print "-"        
    print ""
    print ""
                