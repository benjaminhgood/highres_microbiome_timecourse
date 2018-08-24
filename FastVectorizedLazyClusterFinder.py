import numpy
import sys
import cPickle
from scipy.stats import chi2
from math import log

filename = sys.argv[1]

snp_alignment = cPickle.load( open( filename, "rb" ) )

def calculate_vectorized_distances_between(avg_fs,Ds,ns,f2,D2,n2):
    
    # make vectorized versions
    f1s = numpy.array(avg_fs)
    total_Ds = numpy.array(Ds)
    n1s = numpy.array(ns)
    
    D1s = total_Ds*1.0/n1s[:,None,None] 
    
    if len(f1s)==0:
        return numpy.array([1e300]), numpy.array([False])
    
    #print "Entered function"
    good_idxs = (D1s>0)*(D2[None,:,:]>0)
    dfs = (good_idxs).sum(axis=2)
    
     
    Cs = numpy.power(2*(dfs+1e-300),-0.5)
    
    # Calculate average
    n3s = n1s+n2
    D3s = D1s+D2[None,:,:]
    
    f3s = (f1s*D1s+(f2*D2)[None,:,:])/(D3s+(D3s==0))
    # other polarization
    f3_primes = (f1s*D1s+((1-f2)*D2)[None,:,:])/(D3s+(D3s==0))
    
    
    variances = 2*f3s*(1-f3s)/(D3s+(D3s==0))
    variances_prime = 2*f3_primes*(1-f3_primes)/(D3s+(D3s==0))
    
    raw_distance = ((numpy.square(f2[None,:,:]-f1s)/(variances+(variances==0)))*good_idxs).sum(axis=2)
    distances = Cs*(raw_distance-dfs-2*numpy.log(n3s)[:,None])
    
    raw_distance = ((numpy.square(f1s-(1-f2[None,:,:]))/(variances_prime+(variances_prime==0)))*good_idxs).sum(axis=2)
    distance_primes = Cs*( raw_distance - dfs - 2*numpy.log(n3s)[:,None] )
        
    flips = (distance_primes<distances)
    distances = numpy.fmin(distances, distance_primes)
        
    return distances, flips

# create initial cluster_A, cluster_D matrices
# one for each trajectory

cluster_As = []
cluster_Ds = []
cluster_sizes = []

max_num_snps = snp_alignment.shape[1]
#max_num_snps = 1000
dmax = 9
dmax = 5

sys.stderr.write("Clustering %d snps...\n" % max_num_snps)

for snp_idx in xrange(0,max_num_snps):
    Ds = snp_alignment[:,snp_idx,:].sum(axis=1)
    As = snp_alignment[:,snp_idx,0]
    As = numpy.reshape(As, (1,len(As)))
    Ds = numpy.reshape(Ds, (1,len(Ds)))
    
    cluster_As.append(As)
    cluster_Ds.append(Ds)
    cluster_sizes.append(1)


# returns clustered versions of these things (don't save which ones are in which)
def calculate_cluster_trajectories(cluster_As, cluster_Ds, cluster_sizes):

    # calculate avg fs for singletons
    cluster_avg_fs = []
    for cluster_idx in xrange(0,len(cluster_As)):
        cluster_avg_fs.append( cluster_As[cluster_idx]*1.0 / (cluster_Ds[cluster_idx] + (cluster_Ds[cluster_idx]==0)) )

    # Do everything twice
    for overall_iteration in xrange(0,2):

        finished_cluster_As = []
        finished_cluster_Ds = []
        finished_cluster_avg_fs = []
        finished_cluster_sizes = []

        recently_merged_As = []
        recently_merged_Ds = []
        recently_merged_sizes = []
        recently_merged_avg_fs = []
        recently_merged_distances = []

        iteration = 0
        while (len(cluster_As)+len(recently_merged_As))>0:
    
            iteration+=1
            if iteration%100==0:
            
                sys.stderr.write("Iteration %d...\n" % (iteration))
    
            # Get target trajectory
            # start w/ non-recently merged oens
            if len(cluster_As)>0:
                target_cluster_As = cluster_As.pop()
                target_cluster_Ds = cluster_Ds.pop()
                target_cluster_avg_fs = cluster_avg_fs.pop()
                target_cluster_size = cluster_sizes.pop()
                target_cluster_last_merge = -1e300
            else:
                target_cluster_As = recently_merged_As.pop()
                target_cluster_Ds = recently_merged_Ds.pop()
                target_cluster_avg_fs = recently_merged_avg_fs.pop()
                target_cluster_size = recently_merged_sizes.pop()
                target_cluster_last_merge = recently_merged_distances.pop()
    
            # Look for distance between target and recently merged clusters
            
            best_pair = None
            best_distance = 1e300
        
            distances, flips = calculate_vectorized_distances_between( recently_merged_avg_fs, recently_merged_Ds, recently_merged_sizes, target_cluster_avg_fs, target_cluster_Ds, target_cluster_size)
        
            recent_i = distances.argmin()
            recent_flip = flips[recent_i] 
            best_recent_distance = distances[recent_i]
            best_recent_pair = (recent_i,recent_flip)
        
            if len(recently_merged_avg_fs)>0:
                last_recent_distance = recently_merged_distances[recent_i]
            else:
                last_recent_distance = -1e300
        
            if (best_recent_distance<2) or (best_recent_distance <= last_recent_distance) or (best_recent_distance <= target_cluster_last_merge) or (len(cluster_As)<1):
                
                # pre-emptively merge without comparing rest of trajectories!
                i = -1
                flip = False
                best_distance = 1e300
                best_pair = (i,flip)
            
            else:
                
                # calculate distances with rest of guys!
        
            
                best_pair = None
                best_distance = 1e300
        
                distances, flips = calculate_vectorized_distances_between( cluster_avg_fs, cluster_Ds, cluster_sizes, target_cluster_avg_fs, target_cluster_Ds, target_cluster_size)
        
                i = distances.argmin()
                flip = flips[i] 
                best_distance = distances[i]
                best_pair = (i,flip)
        
            if best_distance < best_recent_distance:
                best_best_distance = best_distance
            else:
                best_best_distance = best_recent_distance 
            
            if best_best_distance > dmax:
                # if bad distance add it to "finished" bin
                finished_cluster_As.append(target_cluster_As)
                finished_cluster_Ds.append(target_cluster_Ds)
                finished_cluster_sizes.append(target_cluster_size)
                finished_cluster_avg_fs.append( target_cluster_avg_fs )
                # Now go to next record
                continue
            
            # otherwise, merge target and best match
            
            if best_distance < best_recent_distance:
                # merge w/ record in cluster_XXX
                other_cluster_As = cluster_As.pop(i)
                other_cluster_Ds = cluster_Ds.pop(i)
                other_cluster_avg_fs = cluster_avg_fs.pop(i)
                other_cluster_size = cluster_sizes.pop(i)
            
                if flip:
                    target_cluster_As = target_cluster_Ds - target_cluster_As
                    
            
                merged_As = target_cluster_As + other_cluster_As
                merged_Ds = target_cluster_Ds + other_cluster_Ds
                merged_size = target_cluster_size + other_cluster_size
                merged_distance = best_distance
                merged_avg_fs = merged_As*1.0/(merged_Ds+(merged_Ds==0))
                
                recently_merged_As.append(merged_As)
                recently_merged_Ds.append(merged_Ds)
                recently_merged_sizes.append(merged_size)
                recently_merged_distances.append( best_distance )
                recently_merged_avg_fs.append(merged_avg_fs)
            
            else:
                
                # merge w/ record in recently merged cluster
                # pre-emptive merge!
                
                if recent_flip:
                    target_cluster_As = target_cluster_Ds - target_cluster_As
                    
                
                recently_merged_As[recent_i] += target_cluster_As
                recently_merged_Ds[recent_i] += target_cluster_Ds
                recently_merged_sizes[recent_i] += target_cluster_size
                recently_merged_avg_fs[recent_i] = recently_merged_As[recent_i]*1.0/(recently_merged_Ds[recent_i]+(recently_merged_Ds[recent_i]==0))
                recently_merged_distances[recent_i] = max([best_recent_distance, target_cluster_last_merge, last_recent_distance]) 
                
            
        # end of loop.. everything should be in finished
        cluster_As = finished_cluster_As
        cluster_Ds = finished_cluster_Ds
        cluster_sizes = finished_cluster_sizes
        cluster_avg_fs = finished_cluster_avg_fs
        
    return cluster_As, cluster_Ds, cluster_sizes

cluster_As, cluster_Ds, cluster_sizes = calculate_cluster_trajectories(cluster_As, cluster_Ds, cluster_sizes)

# calculate averages
cluster_As = numpy.array(cluster_As)
cluster_Ds = numpy.array(cluster_Ds)
cluster_sizes = numpy.array(cluster_sizes)

cluster_avg_fs = cluster_As*1.0/(cluster_Ds+(cluster_Ds==0))
 
sys.stderr.write("Done! Found %d clusters...\n" % (cluster_sizes>1.5).sum())    
    
import pylab
pylab.figure(1,figsize=(7,4))
for cluster_idx in xrange(0,len(cluster_sizes)):
    
    avg_fs = cluster_avg_fs[cluster_idx][0,:]
    avg_Ds = (cluster_Ds[cluster_idx]*1.0 / cluster_sizes[cluster_idx])[0,:]
    avg_good_idxs = (avg_Ds>1)
    
    avg_ts = numpy.arange(0,len(avg_fs))
    avg_ts = avg_ts[avg_good_idxs]
    avg_fs = avg_fs[avg_good_idxs]
    avg_Ds = avg_Ds[avg_good_idxs]
    
    if len(avg_Ds)<1:
        continue
        
    flip = (avg_fs[0]>0.5)
    if flip:
        avg_fs = 1-avg_fs
            
    if cluster_sizes[cluster_idx]>1:
        
        pylab.plot(avg_ts, avg_fs,'-',linewidth=0.5)

pylab.ylim([-0.05,1.05])
pylab.xlabel('Timepoint')
pylab.ylabel('Allele frequency')
pylab.savefig('clustering_example.png',bbox_inches='tight',dpi=300)
#pylab.show()

# 


