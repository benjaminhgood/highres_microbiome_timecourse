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


def calculate_covariance(f1,f2):
    
    f1_avg = f1.mean()
    f2_avg = f2.mean()
    
    return ((f1-f1_avg)*(f2-f2_avg)).sum()
    

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


correlation_matrix = []
for i in xrange(0,xs.shape[0]):
    
    x = xs[i,:]
    D = cluster_Ds[i,:]
    
    good_idxs = (D>0)[None,:]*(cluster_Ds>0)

    big_xs = x[None,:]*good_idxs
    big_ys = xs*good_idxs

    xdoty = (big_xs*big_ys).sum(axis=1)
    xdotx = (big_xs*big_xs).sum(axis=1)
    ydoty = (big_ys*big_ys).sum(axis=1)

    correlation_row = 1.0-numpy.fabs(xdoty/numpy.sqrt(xdotx*ydoty))
    correlation_matrix.append(correlation_row)
    
correlation_matrix = numpy.array(correlation_matrix)

Y = squareform(correlation_matrix)

#Y = pdist(cluster_fs,metric='correlation')
sys.stderr.write("Done!\n")

sys.stderr.write("Clustering...\n")
Z =  linkage(Y, method='average')
sys.stderr.write("Done!\n")

num_clusterss = numpy.arange(2,11)
silhouette_scores = []
for num_clusters in num_clusterss:

    nodes = fcluster(Z, num_clusters, criterion="maxclust")
    S = metrics.silhouette_score(correlation_matrix, nodes, metric = 'precomputed')
    silhouette_scores.append(S)
    
silhouette_scores = numpy.array(silhouette_scores)
num_clusters = num_clusterss[silhouette_scores.argmax()]

nodes = fcluster(Z, num_clusters, criterion="maxclust")
cluster_snp_map = {}
for snp_idx in xrange(0,len(nodes)):
    
    cluster_label = nodes[snp_idx]
    if cluster_label not in cluster_snp_map:
        cluster_snp_map[cluster_label] = []
        
    cluster_snp_map[cluster_label].append(snp_idx)

cluster_fs_map = {}
for cluster_label in cluster_snp_map.keys():
    
    anchor_fs = cluster_fs[cluster_snp_map[cluster_label][0]]
    if anchor_fs[0]>0.5:
        anchor_fs = 1-anchor_fs

    cluster_fs_map[cluster_label] = [anchor_fs]

    if len(cluster_snp_map[cluster_label]) > 1:
    
        for snp_idx in cluster_snp_map[cluster_label][1:]:
            
            target_fs = cluster_fs[snp_idx]
            
            cov = calculate_covariance(anchor_fs, target_fs)
            
            if cov<0:
                
                target_fs = 1-target_fs
                
            cluster_fs_map[cluster_label].append(target_fs)
            
    cluster_fs_map[cluster_label] = numpy.array(cluster_fs_map[cluster_label])
    
print len(cluster_fs_map.keys())


# now plot!
import pylab

pylab.figure(1,figsize=(7,2))
pylab.figure(2) # dummy

ts = numpy.arange(0,cluster_fs.shape[1])

for cluster_label in cluster_fs_map.keys():

    cluster_fs = cluster_fs_map[cluster_label]
    print cluster_fs.shape[0]    
    pylab.figure(2)
    line, = pylab.plot([1],[1],'.')
    color = pylab.getp(line,'color')
    
    avg_fs = numpy.median(cluster_fs, axis=0)
    
    if avg_fs[0]>0.5:
        avg_fs = 1-avg_fs
        cluster_fs = 1-cluster_fs
    
    
    pylab.figure(1)
    for snp_idx in xrange(0,cluster_fs.shape[0]):
        fs = cluster_fs[snp_idx]
        
        #if fs[8]>0.7 and fs[7]<0.3:
        #if numpy.fabs(fs-avg_fs).mean() > 0.2:
        if True:
            pylab.plot(ts,fs,'-',linewidth=0.5,alpha=0.5,color=color,zorder=0)
    
    pylab.plot(ts, avg_fs, 'o-', linewidth=2, color=color, zorder=1,markeredgecolor='k')
pylab.figure(1)
pylab.xlabel('Timepoint')
pylab.ylabel('Frequency')
pylab.ylim([-0.05,1.05])
pylab.savefig('cluster_example.png',bbox_inches='tight',dpi=300)