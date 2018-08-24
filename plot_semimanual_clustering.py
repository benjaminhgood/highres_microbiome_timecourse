# Note to self, this is what I've been using to cluster B. vulgatus

import numpy
import sys
import cPickle
from scipy.stats import chi2
from math import log
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import squareform
from sklearn import metrics 

filename = sys.argv[1]
max_num_snps = long(sys.argv[2])

snp_alignment = cPickle.load( open( filename, "rb" ) )

# create initial cluster_A, cluster_D matrices
# one for each trajectory

#max_num_snps = snp_alignment.shape[1]
#max_num_snps = 100
dmax = 9
dmax = 5

max_num_snps = min([max_num_snps, snp_alignment.shape[1]])

sys.stderr.write("Clustering %d snps...\n" % max_num_snps)
cluster_As = []
cluster_Ds = []
cluster_sizes = []
for snp_idx in xrange(0,max_num_snps):
    Ds = snp_alignment[:,snp_idx,:].sum(axis=1)
    As = snp_alignment[:,snp_idx,0]
    As = numpy.reshape(As, (len(As),))
    Ds = numpy.reshape(Ds, (len(Ds),))
    
    cluster_As.append(As)
    cluster_Ds.append(Ds)
    cluster_sizes.append(1)

cluster_As = numpy.array(cluster_As)
cluster_Ds = numpy.array(cluster_Ds)

cluster_fs = cluster_As*1.0/(cluster_Ds+(cluster_Ds==0))


sys.stderr.write("Calculating distances...\n")

avg_fs = (cluster_fs*(cluster_Ds>0)).mean()

xs = cluster_fs-avg_fs

#good_idxs = (cluster_Ds>0)[:,None,:]*(cluster_Ds>0)[None,:,:]

#big_xs = xs[:,None,:]*good_idxs
#big_ys = xs[None,:,:]*good_idxs

#xdoty = (big_xs*big_ys).sum(axis=2)
#xdotx = (big_xs*big_xs).sum(axis=2)
#ydoty = (big_ys*big_ys).sum(axis=2)

#correlation_matrix = 1.0 - xdoty/numpy.sqrt(xdotx*ydoty)

distance_matrix_1 = []
distance_matrix_2 = []
for i in xrange(0,xs.shape[0]):
    
    target_fs = cluster_fs[i,:]
    target_D = cluster_Ds[i,:]
    
    good_idxs = (target_D>0)[None,:]*(cluster_Ds>0)
    
    dfs = good_idxs.sum(axis=1)
    
    mses = numpy.square(target_fs[None,:]-cluster_fs)*good_idxs
    mse_primes = numpy.square(1-target_fs[None,:]-cluster_fs)*good_idxs
    
    avg_D = (target_D[None,:] + cluster_Ds)/2.0
    avg_f = (target_fs[None,:] + cluster_fs)/2.0
    avg_f_prime = (1-target_fs[None,:] + cluster_fs)/2.0
    
    variances = ((avg_f*(1-avg_f)/(avg_D+(avg_D==0)))*good_idxs)
    variance_primes = ((avg_f_prime*(1-avg_f_prime)/(avg_D+(avg_D==0)))*good_idxs)
    
    
    distances = (mses/(variances+(variances==0))).sum(axis=1)/dfs
    distances_prime = (mse_primes/(variance_primes+(variance_primes==0))).sum(axis=1)/dfs
    
    #distances = numpy.fmin(distances, distances_prime)
    
    distance_matrix_1.append(distances)
    distance_matrix_2.append(distances_prime)
    
    
    
    
distance_matrix_1 = numpy.array(distance_matrix_1)
distance_matrix_2 = numpy.array(distance_matrix_2)
distance_matrix = numpy.fmin(distance_matrix_1,distance_matrix_2)

distance_matrix = (distance_matrix+distance_matrix.T)/2.0

Y = squareform(distance_matrix)

sys.stderr.write("Clustering...\n")
Z =  linkage(Y, method='average')
#Z = linkage(Y, method='complete')
sys.stderr.write("Done!\n")

#num_clusterss = numpy.arange(2,13)
#num_clusterss = numpy.array([3,4])
#silhouette_scores = []
#for num_clusters in num_clusterss:

#    nodes = fcluster(Z, num_clusters, criterion="maxclust")
#    S = metrics.silhouette_score(distance_matrix, nodes, metric = 'precomputed')
#    silhouette_scores.append(S)
    
#    print num_clusters, S
    
#silhouette_scores = numpy.array(silhouette_scores)
#num_clusters = num_clusterss[silhouette_scores.argmax()]

#nodes = fcluster(Z, num_clusters, criterion="maxclust")
nodes = fcluster(Z, 20, criterion="distance")

cluster_snp_map = {}
for snp_idx in xrange(0,len(nodes)):
    
    cluster_label = nodes[snp_idx]
    if cluster_label not in cluster_snp_map:
        cluster_snp_map[cluster_label] = []
        
    cluster_snp_map[cluster_label].append(snp_idx)

cluster_fs_map = {}

cluster_As_map = {}
cluster_Ds_map = {}
cluster_avg_fs_map = {}
cluster_total_Ds_map = {}

for cluster_label in cluster_snp_map.keys():
    
    anchor_idx = cluster_snp_map[cluster_label][0]
    
    cluster_As_map[cluster_label] = [cluster_As[anchor_idx,:]]
    cluster_Ds_map[cluster_label] = [cluster_Ds[anchor_idx,:]]
    
    if len(cluster_snp_map[cluster_label]) > 1:
    
        for snp_idx in cluster_snp_map[cluster_label][1:]:
            
            target_As = cluster_As[snp_idx,:]
            target_Ds = cluster_Ds[snp_idx,:]
            
            if distance_matrix_2[anchor_idx,snp_idx] < distance_matrix_1[anchor_idx,snp_idx]:
                # re-polarize
                target_As = target_Ds-target_As
                
                
            cluster_As_map[cluster_label].append(target_As)
            cluster_Ds_map[cluster_label].append(target_Ds)
            
    
    cluster_As_map[cluster_label] = numpy.array(cluster_As_map[cluster_label])
    cluster_Ds_map[cluster_label] = numpy.array(cluster_Ds_map[cluster_label])
    cluster_total_Ds_map[cluster_label] = cluster_Ds_map[cluster_label].sum(axis=0)
    cluster_avg_fs_map[cluster_label] = cluster_As_map[cluster_label].sum(axis=0)*1.0/(cluster_total_Ds_map[cluster_label]+(cluster_total_Ds_map[cluster_label]==0))
    
    if (cluster_avg_fs_map[cluster_label][0]+cluster_avg_fs_map[cluster_label][1])/2.0 > 0.5:
        cluster_avg_fs_map[cluster_label] = 1-cluster_avg_fs_map[cluster_label]
        cluster_As_map[cluster_label] = cluster_Ds_map[cluster_label] - cluster_As_map[cluster_label]
    
    
print "Num clusters:", len(cluster_avg_fs_map.keys())


# now plot!
import pylab

pylab.figure(1,figsize=(7,2))
pylab.figure(2) # dummy

ts = numpy.arange(0,cluster_fs.shape[1])

for cluster_label in cluster_avg_fs_map.keys():

    cluster_As = cluster_As_map[cluster_label]
    cluster_Ds = cluster_Ds_map[cluster_label]
    cluster_total_Ds = cluster_total_Ds_map[cluster_label]
    cluster_avg_fs = cluster_avg_fs_map[cluster_label]
    
    print cluster_As.shape[0]    

    if cluster_As.shape[0]<1000:
        continue
    
    pylab.figure(2)
    line, = pylab.plot([1],[1],'.')
    color = pylab.getp(line,'color')
    
    
    pylab.figure(1)
    for snp_idx in xrange(0,cluster_As.shape[0]):
        
        fs = cluster_As[snp_idx]*1.0/(cluster_Ds[snp_idx]+(cluster_Ds[snp_idx]==0))
        
        good_idxs = (cluster_Ds[snp_idx]>0)
        masked_fs = fs[good_idxs]
        masked_ts = ts[good_idxs]
        
        pylab.plot(masked_ts, masked_fs,'-',linewidth=0.5,alpha=0.5,color=color,zorder=0)
    
    good_idxs = (cluster_total_Ds>0)
    masked_fs = cluster_avg_fs[good_idxs]
    masked_ts = ts[good_idxs]
    pylab.plot(masked_ts, masked_fs, 'o-', linewidth=2, color=color, zorder=1,markeredgecolor='k')
pylab.figure(1)
pylab.xlabel('Timepoint')
pylab.ylabel('Frequency')
pylab.ylim([-0.05,1.05])
pylab.savefig('semimanual_clustering.png',bbox_inches='tight',dpi=300)
sys.exit(0)
import pylab as plt
pylab.figure(3)
plt.figure(figsize=(25, 10))
plt.title('Hierarchical Clustering Dendrogram')
plt.xlabel('sample index')
plt.ylabel('distance')
dendrogram(
    Z,
    leaf_rotation=90.,  # rotates the x axis labels
    leaf_font_size=8.,  # font size for the x axis labels
)
#pylab.savefig('dendrogram_example.png',bbox_inches='tight',dpi=300)