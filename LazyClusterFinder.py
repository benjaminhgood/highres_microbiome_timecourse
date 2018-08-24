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
    
    C = (df/2)**(0.5)
    
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
    
    
    #scaled_mse = (numpy.square(f2-f1)/(variances+(variances==0)))[good_idxs].sum()

    #scaled_mse_prime = (numpy.square(f2-(1-f1))/(variances_prime+(variances_prime==0)))[good_idxs].sum()

    #print "After mse"
    

    #print scaled_mse, scaled_mse_prime


    #distance = -1*log(chi2.sf(scaled_mse, df)+1e-300)-log(n3)
    #distance_prime = -1*log(chi2.sf(scaled_mse_prime, df)+1e-300)-log(n3)
    
    #print "After distances"
    
    #return distance, False
    
    distance = C*((numpy.square(f2-f1)/(variances+(variances==0)))[good_idxs].mean()-1)
    
    distance_prime = C*((numpy.square(f2-(1-f1))/(variances_prime+(variances_prime==0)))[good_idxs].mean() - 1 )
    
    if distance_prime < distance:
        distance = distance_prime
        flip = True
    else:
        flip = False
    
    #print "Returning"
    
        
    return distance, flip

# create initial cluster_A, cluster_D matrices
# one for each trajectory

cluster_As = []
cluster_Ds = []

#max_num_snps = snp_alignment.shape[1]
max_num_snps = 3000
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

for overall_iteration in xrange(0,2):

    finished_cluster_As = []
    finished_cluster_Ds = []

    iteration = 0
    while len(cluster_As)>0:
    
        iteration+=1
        if iteration%100==0:
            
            sys.stderr.write("Iteration %d...\n" % (iteration))
    
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
        #sys.stderr.write("Done!\n")
    
        #sys.stderr.write("Calculating cluster distances...\n")
    
        best_pair = None
        best_distance = 1e300
    
        target_cluster_As = cluster_As.pop(0)
        target_cluster_Ds = cluster_Ds.pop(0)
        target_cluster_avg_fs = cluster_avg_fs.pop(0)
        target_cluster_avg_Ds = cluster_avg_Ds.pop(0)
        target_cluster_size = cluster_sizes.pop(0)
    
        num_clusters = len(cluster_avg_fs)
            
        for i in reversed(xrange(0,num_clusters)):
        
            distance, flip = calculate_distance_between( cluster_avg_fs[i], cluster_avg_Ds[i], cluster_sizes[i], target_cluster_avg_fs, target_cluster_avg_Ds, target_cluster_size )
            
            if distance < best_distance:
                best_pair = (i,flip)
                best_distance = distance
                
            if best_distance < 3:
                break
                    
        #sys.stderr.write("Done!\n")
        #print best_distance, best_pair            
        
        if best_distance > 9 or (best_pair==None):
            # if bad distance add it to "finished" bin
            finished_cluster_As.append(target_cluster_As)
            finished_cluster_Ds.append(target_cluster_Ds)
        
        else:
            # otherwise, merge and move to end!
            
            # otherwise, merge stuff!
            i,flip = best_pair
            
            other_cluster_As = cluster_As.pop(i)
            other_cluster_Ds = cluster_Ds.pop(i)
    
            if flip:
                # re-polarize
                target_cluster_As = target_cluster_Ds-target_cluster_As
    
            cluster_As.append( numpy.vstack( [other_cluster_As, target_cluster_As] ) )  
            cluster_Ds.append( numpy.vstack( [other_cluster_Ds, target_cluster_Ds] ) )
            
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
            pylab.plot(ts, fs,'-',color=color,alpha=0.5,linewidth=0.5)
            
    if cluster_sizes[cluster_idx]>1:
        
        pylab.plot(avg_ts, avg_fs,'o-',color=color,markeredgecolor='k',linewidth=2)

pylab.ylim([-0.05,1.05])
pylab.xlabel('Timepoint')
pylab.ylabel('Allele frequency')
pylab.savefig('clustering_example.png',bbox_inches='tight',dpi=300)
#pylab.show()

# 


