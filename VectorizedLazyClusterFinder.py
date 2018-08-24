import numpy
import sys
import cPickle
from scipy.stats import chi2
from math import log

filename = sys.argv[1]

snp_alignment = cPickle.load( open( filename, "rb" ) )

def calculate_vectorized_distances_between(f1s,D1s,n1s,f2,D2,n2):
    
    if len(f1s)==0:
        return numpy.array([1e300]), numpy.array([False])
    
    #print "Entered function"
    good_idxs = (D1s>0)*(D2[None,:]>0)
    dfs = (good_idxs).sum(axis=1)
    
    Cs = (1.0/dfs/2.0)**(0.5)
    
    
    # Calculate average
    n3s = n1s+n2
    D3s = D1s+D2[None,:]
    f3s = (f1s*D1s+(f2*D2)[None,:])/(D3s+(D3s==0))
    # other polarization
    f3_primes = (f1s*D1s+((1-f2)*D2)[None,:])/(D3s+(D3s==0))
    
    variances = 2*f3s*(1-f3s)/(D3s+(D3s==0))
    variances_prime = 2*f3_primes*(1-f3_primes)/(D3s+(D3s==0))
    
    #print "After variances"
    
    
    #scaled_mse = (numpy.square(f2-f1)/(variances+(variances==0)))[good_idxs].sum()

    #scaled_mse_prime = (numpy.square(f2-(1-f1))/(variances_prime+(variances_prime==0)))[good_idxs].sum()

    #print "After mse"
    

    #print scaled_mse, scaled_mse_prime


    #distance = -1*log(chi2.sf(scaled_mse, df)+1e-300)-log(n3)
    #distance_prime = -1*log(chi2.sf(scaled_mse_prime, df)+1e-300)-log(n3)
    
    #print "After distances"
    
    #return distance, False
    
    #print Cs.shape
    #print variances.shape
    #print good_idxs.shape
    #print dfs.shape
    
    distances = Cs*(((numpy.square(f2[None,:]-f1s)/(variances+(variances==0)))*good_idxs).sum(axis=1)-dfs-2*numpy.log(n3s))
    
    distance_primes = Cs*(((numpy.square(f1s-(1-f2[None,:]))/(variances_prime+(variances_prime==0)))*good_idxs).sum(axis=1) - dfs - 2*numpy.log(n3s) )
    
    #distances = numpy.clip(distances,0,1e300)
    #distance_primes = numpy.clip(distance_primes,0,1e300)
    
    
    flips = (distance_primes<distances)
    distances = numpy.fmin(distances, distance_primes)
        
    return distances, flips

# create initial cluster_A, cluster_D matrices
# one for each trajectory

cluster_As = []
cluster_Ds = []

#max_num_snps = snp_alignment.shape[1]
max_num_snps = 10000
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

clustering_operations = []

# Do everything twice
for overall_iteration in xrange(0,2):

    finished_cluster_As = []
    finished_cluster_Ds = []

    recently_merged_As = []
    recently_merged_Ds = []
    recently_merged_distances = []

    iteration = 0
    while (len(cluster_As)+len(recently_merged_distances))>0:
    
        iteration+=1
        if iteration%100==0:
            
            sys.stderr.write("Iteration %d...\n" % (iteration))
    
        # Get target trajectory
        # start w/ non-recently merged oens
        if len(cluster_As)>0:
            target_cluster_As = cluster_As.pop(0)
            target_cluster_Ds = cluster_Ds.pop(0)
        else:
            target_cluster_As = recently_merged_As.pop(0)
            target_cluster_Ds = recently_merged_Ds.pop(0)
            recently_merged_distances.pop(0)
    
        target_cluster_avg_Ds = target_cluster_Ds.mean(axis=0) 
    
        total_Ds = target_cluster_Ds.sum(axis=0)
        target_cluster_avg_fs = target_cluster_As.sum(axis=0)*1.0 / (total_Ds+(total_Ds==0)) 
        target_cluster_size = target_cluster_As.shape[0]
    
        # first calculate cluster averages
        
        recently_merged_avg_Ds = []
        recently_merged_avg_fs = []
        recently_merged_sizes = []

        for cluster_idx in xrange(0,len(recently_merged_As)):
        
            recently_merged_avg_Ds.append( recently_merged_Ds[cluster_idx].mean(axis=0) )
    
            total_Ds = recently_merged_Ds[cluster_idx].sum(axis=0)
            recently_merged_avg_fs.append( recently_merged_As[cluster_idx].sum(axis=0)*1.0 / (total_Ds+(total_Ds==0)) )
        
            recently_merged_sizes.append( recently_merged_As[cluster_idx].shape[0] )
        #sys.stderr.write("Done!\n")
    
        #sys.stderr.write("Calculating cluster distances...\n")
    
        recently_merged_avg_fs = numpy.array(recently_merged_avg_fs)
        recently_merged_avg_Ds = numpy.array(recently_merged_avg_Ds)
        recently_merged_sizes = numpy.array(recently_merged_sizes)
    
        best_pair = None
        best_distance = 1e300
        
        distances, flips = calculate_vectorized_distances_between( recently_merged_avg_fs, recently_merged_avg_Ds, recently_merged_sizes, target_cluster_avg_fs, target_cluster_avg_Ds, target_cluster_size)
        
        recent_i = distances.argmin()
        recent_flip = flips[recent_i] 
        best_recent_distance = distances[recent_i]
        best_recent_pair = (recent_i,recent_flip)
        
        if len(recently_merged_avg_fs)>0:
            last_recent_distance = recently_merged_distances[recent_i]
        else:
            last_recent_distance = -1
        
        if (best_recent_distance <= last_recent_distance) or (best_recent_distance<2):
            
            # pre-emptive merge!
            other_cluster_As = recently_merged_As.pop(recent_i)
            other_cluster_Ds = recently_merged_Ds.pop(recent_i)
            
            merged_distance = recently_merged_distances.pop(recent_i)
            merged_flip = recent_flip
            
        else:
        
            # Now we have to compare it to the rest of the trajectories!
            # first calculate cluster averages
        
            cluster_avg_Ds = []
            cluster_avg_fs = []
            cluster_sizes = []

            for cluster_idx in xrange(0,len(cluster_As)):
        
                cluster_avg_Ds.append( cluster_Ds[cluster_idx].mean(axis=0) )
    
                total_Ds = cluster_Ds[cluster_idx].sum(axis=0)
                cluster_avg_fs.append( cluster_As[cluster_idx].sum(axis=0)*1.0 / (total_Ds+(total_Ds==0)) )
        
                cluster_sizes.append( cluster_As[cluster_idx].shape[0] )
        
    
            cluster_avg_fs = numpy.array(cluster_avg_fs)
            cluster_avg_Ds = numpy.array(cluster_avg_Ds)
            cluster_sizes = numpy.array(cluster_sizes)
    
            best_pair = None
            best_distance = 1e300
        
            distances, flips = calculate_vectorized_distances_between( cluster_avg_fs, cluster_avg_Ds, cluster_sizes, target_cluster_avg_fs, target_cluster_avg_Ds, target_cluster_size)
        
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
                # Now go to next record
                continue
            
            # otherwise, merge and move to end!
            
            if best_distance < best_recent_distance:
                # merge w/ record in cluster_XXX
                other_cluster_As = cluster_As.pop(i)
                other_cluster_Ds = cluster_Ds.pop(i)
            
                merged_distance = best_distance
                merged_flip = flip
            
            else:
                # merge w/ record in recently merged cluster
                # pre-emptive merge!
                other_cluster_As = recently_merged_As.pop(recent_i)
                other_cluster_Ds = recently_merged_Ds.pop(recent_i)
            
                merged_distance = best_recent_distance
                merged_flip = recent_flip
            
            
        # Now actually Pull off the merge
            
        if merged_flip:
            target_cluster_As = target_cluster_Ds - target_cluster_As
            
        recently_merged_As.append( numpy.vstack( [other_cluster_As, target_cluster_As] ) )  
        recently_merged_Ds.append( numpy.vstack( [other_cluster_Ds, target_cluster_Ds] ) )
        recently_merged_distances.append( merged_distance )
            
    # end of loop.. everything should be in finished
    cluster_As = finished_cluster_As
    cluster_Ds = finished_cluster_Ds
    
# calculate one last set of averages    
cluster_avg_Ds = []
cluster_avg_fs = []
cluster_sizes = []

# first calculate cluster averages
#sys.stderr.write("Calculating cluster averages...\n")
for cluster_idx in xrange(0,len(cluster_As)):
        
    cluster_avg_Ds.append( cluster_Ds[cluster_idx].mean(axis=0) )
    
    total_Ds = cluster_Ds[cluster_idx].sum(axis=0)
    cluster_avg_fs.append( cluster_As[cluster_idx].sum(axis=0)*1.0 / (total_Ds+(total_Ds==0)) )
        
    cluster_sizes.append( cluster_As[cluster_idx].shape[0] )
 
sys.stderr.write("Done! Found %d clusters...\n" % len(cluster_sizes))    
    
import pylab
pylab.figure(1,figsize=(7,4))
for cluster_idx in xrange(0,len(cluster_sizes)):
    
    avg_fs = cluster_avg_fs[cluster_idx]
    avg_Ds = cluster_avg_Ds[cluster_idx]
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
    
    if avg_fs.max()<0.1:
        continue
    
    
    color = None    
    for snp_idx in xrange(0,len(cluster_As[cluster_idx])):
        
        As = cluster_As[cluster_idx][snp_idx]
        Ds = cluster_Ds[cluster_idx][snp_idx]
        
        if flip:
            As = Ds-As
        
        fs = As*1.0/(Ds+(Ds==0))
        ts = numpy.arange(0,len(fs))
        good_idxs = (Ds>0)
            
        fs = fs[good_idxs]
        ts = ts[good_idxs]
        
        if cluster_sizes[cluster_idx]==1:
            pylab.plot(ts, fs,'-',color='0.7',alpha=0.5,zorder=0,linewidth=0.5)
        elif color==None:
            line, = pylab.plot(ts, fs,'-',linewidth=0.5,alpha=0.5)
            color = pylab.getp(line,'color')
        else:
            pass
            #pylab.plot(ts, fs,'-',color=color,alpha=0.5,linewidth=0.5)
            
    if cluster_sizes[cluster_idx]>1:
        
        pylab.plot(avg_ts, avg_fs,'o-',color=color,markeredgecolor='k',linewidth=2)

pylab.ylim([-0.05,1.05])
pylab.xlabel('Timepoint')
pylab.ylabel('Allele frequency')
pylab.savefig('clustering_example.png',bbox_inches='tight',dpi=300)
#pylab.show()

# 


