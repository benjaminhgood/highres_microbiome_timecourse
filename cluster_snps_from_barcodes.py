import sys
import cPickle
import numpy
from math import fabs
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import squareform
from sklearn import metrics 


#####
#
# Load precalculated fixation events
#
####


min_distance = 1000

distance_bins = numpy.array([0,100]+[500*i for i in xrange(1,21)]+[2e04,1e05,1e07])
distances = distance_bins[1:]
distance_counts = numpy.zeros_like(distances)*1.0
distance_fractions = numpy.zeros_like(distance_counts)
distance_barcodes = numpy.zeros_like(distance_counts)
all_distance_counts = numpy.zeros_like(distances)*1.0

# create output filename
pickle_filename = sys.argv[1]

blacklisted_longsnp_ids = set()
fourgamete_pairs = set()

sys.stderr.write("Loading pickle from disk...\n")
d = cPickle.load( open( pickle_filename,"rb") )
sys.stderr.write("Done!\n")

id_longsnp_map = d['id_longsnp_map']
allele_idx_map = d['allele_idx_map']
gamete_idx_map = d['gamete_idx_map']
d = d['snp_barcode_timecourse']

# Bad combos:
# AA and AR
# RA and RR
# AA and RA
# RR and RR

# Good combos:
# AA and RR
# AR and RA

# 4 choose 2 = 6 


# mapping for distance matrix
idx_id_map = []
id_idx_map = {}

dmax = 1e14
d0 = 1e04

test_id = d.keys()[0]
#print id_longsnp_map[test_id]

# First figure out which snps to remove
for focal_longsnp_id in d:
      
    #if id_longsnp_map[focal_longsnp_id]==("NC_009614", 263271):
    #    print d[focal_longsnp_id]
              
    # Can't look if there are no linked SNVs
    if len(d[focal_longsnp_id]['longsnps'])==0:
        continue
    
    if focal_longsnp_id not in d[focal_longsnp_id]['longsnps']:
        continue
            
    self_gametes = d[focal_longsnp_id]['longsnps'][focal_longsnp_id]
    
    self_barcodes = self_gametes.sum()
    
    self_bad_barcodes = self_barcodes - self_gametes[gamete_idx_map[('A','A')]] - self_gametes[gamete_idx_map[('R','R')]]
    
    if self_bad_barcodes>0.5:
        continue

    location = id_longsnp_map[focal_longsnp_id][1]    
    
    if location > 1000000:
        continue
        
    idx = len(idx_id_map)
    idx_id_map.append( focal_longsnp_id )
    id_idx_map[focal_longsnp_id] = idx
 
num_snps = len(idx_id_map)
#print num_snps

distance_matrix = numpy.ones((num_snps,num_snps))*d0

for idx in xrange(0,num_snps):
    focal_longsnp_id = idx_id_map[idx]

    # Can't look if there are no linked SNVs
    if len(d[focal_longsnp_id]['longsnps'])==0:
        continue
        
    location = id_longsnp_map[focal_longsnp_id][1]    
    
    alleles = d[focal_longsnp_id]['all']
        
    total_alleles = alleles.sum()    
        
    self_gametes = d[focal_longsnp_id]['longsnps'][focal_longsnp_id]
    
    for other_longsnp_id in d[focal_longsnp_id]['longsnps']:
    
        if other_longsnp_id==focal_longsnp_id:
            continue # don't look at mapping to itself
    
        if other_longsnp_id not in id_idx_map:
            continue
            
        other_idx = id_idx_map[other_longsnp_id]
                
        other_location = id_longsnp_map[other_longsnp_id][1]
        distance = fabs(other_location-location) 
        
        # could have distance condition here...
            
        gametes = d[focal_longsnp_id]['longsnps'][other_longsnp_id]
    
        total_barcodes = gametes.sum()
        
        fraction = total_barcodes*1.0/total_alleles
        
        gamete_freqs = gametes*1.0/total_barcodes
        
        # minimum number of barcodes
        if total_barcodes < 10:
            continue
        
        gamete_idxs = (gametes>0.5)
        num_gametes = gamete_idxs.sum()
        stringent_gamete_idxs = ((gametes>2.5)*(gamete_freqs>0.05))
        stringent_num_gametes = stringent_gamete_idxs.sum()
    
        
        if num_gametes >= 2.5:
            distance = dmax
            
        elif num_gametes==2:
        
            # First check whether there are bad combinations of gametes
            if not ((gamete_idxs[gamete_idx_map[('A','A')]] and gamete_idxs[gamete_idx_map[('R','R')]]) or (gamete_idxs[gamete_idx_map[('A','R')]] and gamete_idxs[gamete_idx_map[('R','A')]])):
                # bad!
                distance = dmax
            else:
                
                if stringent_num_gametes==2:
                    # good to go!
                    distance = d0/2 - total_barcodes
                else:
                    distance = d0 #- total_barcodes
        
        else: # (num_gametes==1)
            
            distance = d0 #- total_barcodes
            
        distance_matrix[idx,other_idx] = distance
        distance_matrix[other_idx,idx] = distance
        

# Now do clustering        
distance_matrix = (distance_matrix+distance_matrix.T)/2.0
numpy.fill_diagonal(distance_matrix,0)
Y = squareform(distance_matrix)

sys.stderr.write("Clustering...")
Z =  linkage(Y, method='complete')
sys.stderr.write("Done!\n")

sys.stderr.write("Forming flat clusters...\n")
nodes = fcluster(Z, d0-1e-09, criterion="distance")
sys.stderr.write("Done!\n")

sys.stderr.write("Writing output...\n")
num_realized_clusters = len(set(nodes))
cluster_idx_map = {}
for idx in xrange(0,len(nodes)):
    
    cluster_label = nodes[idx]
    if cluster_label not in cluster_idx_map:
        cluster_idx_map[cluster_label] = []
    cluster_idx_map[cluster_label].append(idx)

sorted_cluster_labels = sorted(cluster_idx_map.keys(), key=lambda x: len(cluster_idx_map[x]), reverse=True)

sys.stderr.write("%d clusters\n" % len(sorted_cluster_labels))
for cluster_label in sorted_cluster_labels:

    if len(cluster_idx_map[cluster_label])>1.5:
        sys.stderr.write("%d\n" % len(cluster_idx_map[cluster_label]))
    
        record_strs = ["%s|%d" % id_longsnp_map[idx_id_map[idx]] for idx in cluster_idx_map[cluster_label]]
        print " ".join(record_strs)
    
sys.stderr.write("Done!\n")