import numpy
import sys
import cPickle
from scipy.stats import chi2
from math import log

filename = sys.argv[1]

snp_alignment = cPickle.load( open( filename, "rb" ) )

def calculate_distance_between(f1,D1,n1,f2,D2,n2):
    
    #print "Entered function"
    good_idxs = (D1>0)*(D2>0)
    df = (good_idxs).sum()
    
    if df<0.5:
        return 1e03,False
    
    # Calculate average
    n3 = n1+n2
    D3 = D1+D2
    f3 = (f1*D1+f2*D2)/(D3+(D3==0))
    
    # other polarization
    f3_prime = (f1*D1+(1-f2)*D2)/(D3+(D3==0))
    
    # Estimated variances
    
    #print "After average"
    
    # (factor of two from paired end reads?)
    variances = 2*f3*(1-f3)/(D3+(D3==0))
    variances_prime = 2*f3_prime*(1-f3_prime)/(D3+(D3==0))
    
    #print "After variances"
    
    
    scaled_mse = (numpy.square(f2-f1)/(variances+(variances==0)))[good_idxs].sum()

    scaled_mse_prime = (numpy.square(f2-(1-f1))/(variances_prime+(variances_prime==0)))[good_idxs].sum()

    #print "After mse"
    

    #print scaled_mse, scaled_mse_prime


    distance = -1*log(chi2.sf(scaled_mse, df)+1e-300)-log(n3)
    distance_prime = -1*log(chi2.sf(scaled_mse_prime, df)+1e-300)-log(n3)
    
    #print "After distances"
    
    #return distance, False
    
    if distance_prime < distance:
        distance = distance_prime
        flip = True
    else:
        flip = False
    
    #print "Returning"
    
        
    return distance, flip

def calculate_best_distance(cluster_avg_fs, cluster_avg_Ds, cluster_sizes):
    best_pair = None
    best_distance = 1e300
    num_clusters = len(cluster_avg_fs)
    
    #print len(cluster_sizes), len(cluster_avg_fs), len(cluster_avg_Ds)    
    
    for i in xrange(0,num_clusters):
        for j in xrange(i+1, num_clusters):
            
            distance, flip = calculate_distance_between( cluster_avg_fs[i], cluster_avg_Ds[i], cluster_sizes[i], cluster_avg_fs[j], cluster_avg_Ds[j], cluster_sizes[j] )
            
            if distance < best_distance:
                best_pair = (i,j,flip)
                best_distance = distance
                
            if best_distance < 3:
                sys.stderr.write("Terminating early!\n")
                pass
                return best_distance, best_pair

    return best_distance, best_pair

# create initial cluster_A, cluster_D matrices
# one for each trajectory

cluster_As = []
cluster_Ds = []

max_num_snps = snp_alignment.shape[1]
max_num_snps = 200

for snp_idx in xrange(0,max_num_snps):
    Ds = snp_alignment[:,snp_idx,:].sum(axis=1)
    As = snp_alignment[:,snp_idx,0]
    As = numpy.reshape(As, (1,len(As)))
    Ds = numpy.reshape(Ds, (1,len(Ds)))
    
    cluster_As.append(As)
    cluster_Ds.append(Ds)

clustering_operations = []

while True:

    sys.stderr.write("Iteration %d...\n" % (len(clustering_operations)+1))
    
    cluster_avg_Ds = []
    cluster_avg_fs = []
    cluster_sizes = []

    # first calculate cluster averages
    sys.stderr.write("Calculating cluster averages...\n")
    for cluster_idx in xrange(0,len(cluster_As)):
        
        cluster_avg_Ds.append( cluster_Ds[cluster_idx].mean(axis=0) )
    
        total_Ds = cluster_Ds[cluster_idx].sum(axis=0)
        cluster_avg_fs.append( cluster_As[cluster_idx].sum(axis=0)*1.0 / (total_Ds+(total_Ds==0)) )
        
        cluster_sizes.append( cluster_As[cluster_idx].shape[0] )
    sys.stderr.write("Done!\n")
    
    sys.stderr.write("Calculating cluster distances...\n")
    
    best_distance, best_pair = calculate_best_distance(cluster_avg_fs, cluster_avg_Ds, cluster_sizes)
    
    sys.stderr.write("Done!\n")
    print best_distance, best_pair            
    # stop if not a good distance
    if best_distance > 9 or (best_pair==None):
        break
        
    # otherwise, merge stuff!
    i,j,flip = best_pair
    
    clustering_operations.append((i,j,flip, best_distance))
    
    to_merge_As = cluster_As[j]
    to_merge_Ds = cluster_Ds[j]
    
    if flip:
        to_merge_As = to_merge_Ds-to_merge_As
    
    cluster_As[i] = numpy.vstack( [cluster_As[i],to_merge_As] )
    cluster_Ds[i] = numpy.vstack( [cluster_Ds[i],to_merge_Ds] )
    del cluster_As[j]
    del cluster_Ds[j]
    
# cluster avg fs and cluster sizes are set well! 

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
            pylab.plot(ts, fs,'-',color=color,alpha=0.5,linewidth=0.5)
            
    if cluster_sizes[cluster_idx]>1:
        
        pylab.plot(avg_ts, avg_fs,'o-',color=color,markeredgecolor='k',linewidth=2)

pylab.ylim([-0.05,1.05])
pylab.xlabel('Timepoint')
pylab.ylabel('Allele frequency')
pylab.savefig('clustering_example.png',bbox_inches='tight',dpi=300)
#pylab.show()

# 


