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

num_gametes = [0,0,0,0,0] # 1,2,3,4

clusters = []

blacklisted_longsnps = set()

sys.stderr.write("Loading pickle from disk...\n")
d = cPickle.load( open( pickle_filename,"rb") )
sys.stderr.write("Done!\n")

for focal_longsnp in d:
    
    focal_species, focal_snp_str = focal_longsnp
      
    # Can't look if there are no linked SNVs
    if len(d[focal_longsnp]['longsnps'])==0:
        continue
    
    if focal_longsnp not in d[focal_longsnp]['longsnps']:
        continue
        
        
    gametes = d[focal_longsnp]['longsnps'][focal_longsnp]
    
    if ((('A','R') in gametes) and (('A','A') in gametes)) or ((('R','R') in gametes) and (('R','A') in gametes)):
        # likely a duplication or something weird
        num_gametes[0]+=1
        blacklisted_longsnps.add(focal_longsnp)
        continue
    
    for other_longsnp in d[focal_longsnp]['longsnps']:
    
        if other_longsnp==focal_longsnp:
            continue # don't look at mapping to itself
        
        if other_longsnp in blacklisted_longsnps:
            continue # don't look at bad ones!
        
        #print other_longsnp
            
        gametes = d[focal_longsnp]['longsnps'][other_longsnp]
    
        if len(gametes)==1:
            num_gametes[1]+=1
        elif len(gametes)==3:
            num_gametes[3]+=1
        elif len(gametes)==4:
            num_gametes[4]+=1
        elif len(gametes)==2:
            
            # could be two gametes, or secretly three
            if ((('A','R') in gametes) and (('A','A') in gametes)) or ((('R','R') in gametes) and (('R','A') in gametes)):
                # secretly three
                num_gametes[3]+=1
            else:
                num_gametes[2]+=1
                
                # this is a good one:
                found_cluster=False
                for cluster in clusters:
                    if (focal_longsnp in cluster) or (other_longsnp in cluster):
                        cluster.add(focal_longsnp)
                        cluster.add(other_longsnp)
                        found_cluster=True
                        break
                if not found_cluster:
                    clusters.append(set([focal_longsnp, other_longsnp]))
                    

# now go through and merge clusters
cluster_merge=True
while cluster_merge:
    cluster_merge=False
    for i in xrange(0,len(clusters)):
        for j in xrange(i+1,len(clusters)):
            if len(clusters[i] & clusters[j])>0:
                clusters[i].update(clusters[j])
                clusters[j] = set()   
                cluster_merge=True
                
new_clusters = []
for cluster in clusters:
    if len(cluster)>0:
        new_clusters.append(cluster)

# Sort by length
clusters = list(sorted(new_clusters,key=lambda x: len(x),reverse=True))        

#print num_gametes
#print len(clusters), "clusters:"
#for cluster in clusters:
#    print len(cluster)

for cluster in clusters:
    record_strs = []
    for longsnp in cluster:
        record_strs.append(longsnp[1])
    print " ".join(record_strs)