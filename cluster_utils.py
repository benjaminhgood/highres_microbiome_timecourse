import numpy
import sys
import cPickle
from scipy.stats import chi2
from math import log
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import squareform
from sklearn import metrics 

def calculate_distance_matrix(cluster_As, cluster_Ds):
    
    cluster_fs = cluster_As*1.0/(cluster_Ds+(cluster_Ds==0))

    distance_matrix_1 = [] # distance matrix for one polarization
    distance_matrix_2 = [] # distance matrix for the other
        
    for i in xrange(0,cluster_fs.shape[0]):
    
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
    
        distances = (mses/(variances+(variances==0))).sum(axis=1)/(dfs+(dfs==0))
        distances_prime = (mse_primes/(variance_primes+(variance_primes==0))).sum(axis=1)/(dfs+(dfs==0))

        distances[dfs==0] = 100
        distances_prime[dfs==0] = 100

        distance_matrix_1.append(distances)
        distance_matrix_2.append(distances_prime)
    
    distance_matrix_1 = numpy.array(distance_matrix_1)
    distance_matrix_2 = numpy.array(distance_matrix_2)
    distance_matrix = numpy.fmin(distance_matrix_1,distance_matrix_2)

    distance_matrix = (distance_matrix+distance_matrix.T)/2.0

    return distance_matrix, distance_matrix_1, distance_matrix_2


def cluster_snps(cluster_As, cluster_Ds, max_num_snps_to_cluster=2000):

    cluster_fs = cluster_As*1.0/(cluster_Ds+(cluster_Ds==0))

    sys.stderr.write("Calculating distances for %d snps...\n" % len(cluster_As))
     
    distance_matrix, distance_matrix_1, distance_matrix_2 = calculate_distance_matrix(cluster_As, cluster_Ds)   

    Y = squareform(distance_matrix)
    sys.stderr.write("Done!\n")
        

    if cluster_As.shape[0]>max_num_snps_to_cluster or cluster_As.shape[0]<2.5:
        # too many.. put them all in one big cluster
        nodes = numpy.ones(cluster_fs.shape[0])
    else:
        # do hierarchical clustering
        
        sys.stderr.write("SciPy hierarchical clustering...\n")
        #Z =  linkage(Y, method='average')
        Z = linkage(Y, method='complete')
        sys.stderr.write("Done!\n")

        max_num_clusters = min([4,cluster_As.shape[0]-1])
        num_clusterss = numpy.arange(2,max_num_clusters+1)
        silhouette_scores = []
        for num_clusters in num_clusterss:

            nodes = fcluster(Z, num_clusters, criterion="maxclust")
            num_realized_clusters = len(set(nodes))
            
            if num_realized_clusters==1:
                S = 0
            else:
                S = metrics.silhouette_score(distance_matrix, nodes, metric = 'precomputed')
            silhouette_scores.append(S)
    
            print num_clusters, num_realized_clusters, S
    
        silhouette_scores = numpy.array(silhouette_scores)
        num_clusters = num_clusterss[silhouette_scores.argmax()]
        Smax = silhouette_scores.max()
        print num_clusters, Smax
        if Smax < 0:
            nodes = numpy.ones(distance_matrix.shape[0])
        else:
            nodes = fcluster(Z, num_clusters, criterion="maxclust")
        
    # Now figure out polarizations and centroids
        
        
    cluster_snp_map = {}
    for snp_idx in xrange(0,len(nodes)):
    
        cluster_label = nodes[snp_idx]
        if cluster_label not in cluster_snp_map:
            cluster_snp_map[cluster_label] = []
        cluster_snp_map[cluster_label].append(snp_idx)

    snp_flip_map = {snp_idx: False for snp_idx in xrange(0,len(nodes))}
        
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
                    snp_flip_map[snp_idx] = not snp_flip_map[snp_idx]
                
                cluster_As_map[cluster_label].append(target_As)
                cluster_Ds_map[cluster_label].append(target_Ds)
            
    
        cluster_As_map[cluster_label] = numpy.array(cluster_As_map[cluster_label])
        cluster_Ds_map[cluster_label] = numpy.array(cluster_Ds_map[cluster_label])
        cluster_total_Ds_map[cluster_label] = cluster_Ds_map[cluster_label].sum(axis=0)
        cluster_avg_fs_map[cluster_label] = cluster_As_map[cluster_label].sum(axis=0)*1.0/(cluster_total_Ds_map[cluster_label]+(cluster_total_Ds_map[cluster_label]==0))
    
        # now polarize whole cluster if necessary
        if (cluster_avg_fs_map[cluster_label][0]+cluster_avg_fs_map[cluster_label][1])/2.0 > 0.5:
            cluster_avg_fs_map[cluster_label] = 1-cluster_avg_fs_map[cluster_label]
            cluster_As_map[cluster_label] = cluster_Ds_map[cluster_label] - cluster_As_map[cluster_label]
            for snp_idx in cluster_snp_map[cluster_label]:
                snp_flip_map[snp_idx] = not snp_flip_map[snp_idx]    
    
    # now write output
    cluster_map = {}
    for cluster_label in cluster_snp_map:
        cluster_map[cluster_label] = {}  
        cluster_map[cluster_label]['centroid'] = (cluster_avg_fs_map[cluster_label], cluster_total_Ds_map[cluster_label])
        cluster_map[cluster_label]['snps'] = []
        for snp_idx in cluster_snp_map[cluster_label]:
            cluster_map[cluster_label]['snps'].append((snp_idx, snp_flip_map[snp_idx]))
                
    return cluster_map
            

# Max d inferred from barcode linkage
def cluster_snps_by_distance(cluster_As, cluster_Ds, max_num_snps_to_cluster=2000, max_d=20):

    cluster_fs = cluster_As*1.0/(cluster_Ds+(cluster_Ds==0))

    sys.stderr.write("Calculating distances for %d snps...\n" % len(cluster_As))
     
    distance_matrix, distance_matrix_1, distance_matrix_2 = calculate_distance_matrix(cluster_As, cluster_Ds)   

    Y = squareform(distance_matrix)
    sys.stderr.write("Done!\n")
        

    if cluster_As.shape[0]<2.5: # cluster_As.shape[0]>max_num_snps_to_cluster or 
        # too few.. put them all in one big cluster
        nodes = numpy.ones(cluster_fs.shape[0])
    else:
        # do hierarchical clustering
        
        sys.stderr.write("SciPy hierarchical clustering...\n")
        Z =  linkage(Y, method='average')
        #Z = linkage(Y, method='complete')
        sys.stderr.write("Done!\n")

        nodes = fcluster(Z, 20, criterion="distance")

    # Now figure out polarizations and centroids    
        
    cluster_snp_map = {}
    for snp_idx in xrange(0,len(nodes)):
    
        cluster_label = nodes[snp_idx]
        if cluster_label not in cluster_snp_map:
            cluster_snp_map[cluster_label] = []
        cluster_snp_map[cluster_label].append(snp_idx)

    snp_flip_map = {snp_idx: False for snp_idx in xrange(0,len(nodes))}
        
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
                    snp_flip_map[snp_idx] = not snp_flip_map[snp_idx]
                
                cluster_As_map[cluster_label].append(target_As)
                cluster_Ds_map[cluster_label].append(target_Ds)
            
    
        cluster_As_map[cluster_label] = numpy.array(cluster_As_map[cluster_label])
        cluster_Ds_map[cluster_label] = numpy.array(cluster_Ds_map[cluster_label])
        cluster_total_Ds_map[cluster_label] = cluster_Ds_map[cluster_label].sum(axis=0)
        cluster_avg_fs_map[cluster_label] = cluster_As_map[cluster_label].sum(axis=0)*1.0/(cluster_total_Ds_map[cluster_label]+(cluster_total_Ds_map[cluster_label]==0))
    
        # now polarize whole cluster if necessary
        if (cluster_avg_fs_map[cluster_label][0]+cluster_avg_fs_map[cluster_label][1])/2.0 > 0.5:
            cluster_avg_fs_map[cluster_label] = 1-cluster_avg_fs_map[cluster_label]
            cluster_As_map[cluster_label] = cluster_Ds_map[cluster_label] - cluster_As_map[cluster_label]
            for snp_idx in cluster_snp_map[cluster_label]:
                snp_flip_map[snp_idx] = not snp_flip_map[snp_idx]    
    
    # now write output
    cluster_map = {}
    for cluster_label in cluster_snp_map:
        cluster_map[cluster_label] = {}  
        cluster_map[cluster_label]['centroid'] = (cluster_avg_fs_map[cluster_label], cluster_total_Ds_map[cluster_label])
        cluster_map[cluster_label]['snps'] = []
        for snp_idx in cluster_snp_map[cluster_label]:
            cluster_map[cluster_label]['snps'].append((snp_idx, snp_flip_map[snp_idx]))
                
    return cluster_map