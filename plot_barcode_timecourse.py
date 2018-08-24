import cPickle

#filenames = ["antibiotic_genes_to_track.p"]
filenames = ["final_genes_to_track.p"]

for filename in filenames:
    print filename
    d = cPickle.load( open( filename,"rb") )

    num_total_genes = 0
    num_good_genes = 0

    linked_species_map = {}

    for gene_name in d:
    
        num_total_genes += 1
        all_barcodes = d[gene_name]['all']
    
        for species_name in d[gene_name]['species']:
        
            species_barcodes = d[gene_name]['species'][species_name]['any']
        
            #print species_barcodes
            one_good_gene = False
        
            for other_gene_name in d[gene_name]['species'][species_name]:
            
                if other_gene_name=='any':
                    continue
                
                barcodes = d[gene_name]['species'][species_name][other_gene_name]
            
                if ((all_barcodes>=5)*(barcodes>=5)*(barcodes>=0.1*all_barcodes)).any():
                #if ((all_barcodes>=5)*(barcodes>=6)*(barcodes>=0.5*all_barcodes)).any():
                    one_good_gene = True
                    #print barcodes*1.0/all_barcodes, species_barcodes*1.0/all_barcodes
                
            if one_good_gene: # or filename.startswith('fake'):
                num_good_genes += 1
                if species_name not in linked_species_map:
                    linked_species_map[species_name] = 0
   
                #if species_name.startswith('Alistipes_onderdonkii'):
                #if species_name.startswith('Bacteroides_caccae'):
                if species_name.startswith('Bacteroides_vulgatus'):
                    print gene_name
                    print d[gene_name]['all']
                    for other_species_name in d[gene_name]['species']:
                        print  other_species_name, d[gene_name]['species'][other_species_name]
                
                linked_species_map[species_name] += 1
            
    print num_total_genes, num_good_genes
    print linked_species_map