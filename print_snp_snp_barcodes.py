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
for focal_longsnp in d:
    
    focal_species, focal_snp_str = focal_longsnp
      
    
    if len(d[focal_longsnp]['longsnps'])==0:
        continue
    
    print "Focal SNP: %s, %s" % (focal_species, focal_snp_str)
    print "Total barcodes:"
    for allele in d[focal_longsnp]['all']:
        all_barcodes = d[focal_longsnp]['all'][allele]
        print "%s (%d): %s" % (allele, all_barcodes.sum(), ", ".join([str(n) for n in all_barcodes]))
    
    other_longsnp_list = []
    other_total_barcodes = []
    other_total_gametes = []
    
    for other_longsnp in d[focal_longsnp]['longsnps']:
        total_barcodes = 0
        for gamete in d[focal_longsnp]['longsnps'][other_longsnp]:
            total_barcodes += d[focal_longsnp]['longsnps'][other_longsnp][gamete][0].sum()
        
        other_longsnp_list.append(other_longsnp)    
        other_total_barcodes.append(total_barcodes)
        other_total_gametes.append(len(d[focal_longsnp]['longsnps'][other_longsnp]))
        
    dummy1, dummy2, sorted_longsnp_list = zip(*sorted(zip(other_total_gametes, other_total_barcodes, other_longsnp_list), key=lambda pair: (pair[1],pair[0]), reverse=True))

    for other_longsnp in sorted_longsnp_list:
        if focal_longsnp[1]==other_longsnp[1]:
            # mapping to self
            # we don't need this unless it's interesting
            if len(d[focal_longsnp]['longsnps'][other_longsnp])>2:
                print "*** %s->self" % focal_longsnp[1]
                for gamete in d[focal_longsnp]['longsnps'][other_longsnp]:
                
                    barcodes, pvalues = d[focal_longsnp]['longsnps'][other_longsnp][gamete]
            
                    print "%s%s (%d): %s" % (gamete[0], gamete[1], barcodes.sum(), ", ".join(["%d (%g)" % (n,p) for n,p in zip(barcodes, pvalues)])) 
                print "|"   
                break
            continue
                
        print "%s->%s" % (focal_longsnp[1], other_longsnp[1])
        for gamete in d[focal_longsnp]['longsnps'][other_longsnp]:
                
            barcodes, pvalues = d[focal_longsnp]['longsnps'][other_longsnp][gamete]
            
            print "%s%s (%d): %s" % (gamete[0], gamete[1], barcodes.sum(), ", ".join(["%d (%g)" % (n,p) for n,p in zip(barcodes, pvalues)]))    
          
        
        print "|" 
    
    print ""
        