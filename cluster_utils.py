import numpy
import sys
import cPickle
from scipy.stats import chi2
from math import log
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import squareform
from sklearn import metrics 
import config
from numpy.random import random

def calculate_secondary_distance_matrix(target_avg_fs, target_avg_Ds, cluster_As, cluster_Ds):
    safe_cluster_Ds = cluster_Ds+(cluster_Ds==0)
    cluster_fs = cluster_As*1.0/(safe_cluster_Ds)
    distance_matrix_1 = []
    distance_matrix_2 = []
    for i in xrange(0,target_avg_fs.shape[0]):
    
        #target_fs = target_avg_fs[i,:]
        #target_D = target_avg_Ds[i,:]
    
        target_fs = target_avg_fs[i,:]
        target_Ds = target_avg_Ds[i,:]
        target_As = target_fs*target_Ds
        safe_target_Ds = (target_Ds+(target_Ds==0))
        
    
        good_idxs = (target_Ds>0)[None,:]*(cluster_Ds>0)
    
        dfs = good_idxs.sum(axis=1)
    
        mses = numpy.square(target_fs[None,:]-cluster_fs)*good_idxs
        mse_primes = numpy.square(1-target_fs[None,:]-cluster_fs)*good_idxs
    
        # New version
        total_inverse_D = 1.0/safe_target_Ds[None,:]+1.0/safe_cluster_Ds
        avg_fs = (target_As[None,:]+cluster_As)*1.0/(safe_target_Ds[None,:]+safe_cluster_Ds)
        avg_f_primes = (target_Ds[None,:]-target_As[None,:]+cluster_As)*1.0/(safe_target_Ds[None,:]+safe_cluster_Ds)
        
        variances = (avg_fs*(1-avg_fs)*total_inverse_D )*good_idxs
        variance_primes = (avg_f_primes*(1-avg_f_primes)*total_inverse_D )*good_idxs
        
        # old version
        #avg_D = (target_D[None,:] + cluster_Ds)/2.0
        #avg_f = (target_fs[None,:] + cluster_fs)/2.0
        #avg_f_prime = (1-target_fs[None,:] + cluster_fs)/2.0
        #variances = ((avg_f*(1-avg_f)/(avg_D+(avg_D==0)))*good_idxs)
        #variance_primes = ((avg_f_prime*(1-avg_f_prime)/(avg_D+(avg_D==0)))*good_idxs)
    
        
    
        distances = (mses/(variances+(variances==0))).sum(axis=1)/(dfs+(dfs==0))
        distances_prime = (mse_primes/(variance_primes+(variance_primes==0))).sum(axis=1)/(dfs+(dfs==0))

        distances[dfs==0] = 100
        distances_prime[dfs==0] = 100

        distance_matrix_1.append(distances)
        distance_matrix_2.append(distances_prime)
    
    distance_matrix_1 = numpy.array(distance_matrix_1)
    distance_matrix_2 = numpy.array(distance_matrix_2)
    distance_matrix = numpy.fmin(distance_matrix_1,distance_matrix_2)

    #distance_matrix = (distance_matrix+distance_matrix.T)/2.0

    return distance_matrix, distance_matrix_1, distance_matrix_2
    

        

def calculate_distance_matrix(cluster_As, cluster_Ds):
    
    safe_cluster_Ds = cluster_Ds+(cluster_Ds==0)
    
    cluster_fs = cluster_As*1.0/(safe_cluster_Ds)

    distance_matrix_1 = [] # distance matrix for one polarization
    distance_matrix_2 = [] # distance matrix for the other
        
    for i in xrange(0,cluster_fs.shape[0]):
    
        target_fs = cluster_fs[i,:]
        target_As = cluster_As[i,:]
        target_Ds = cluster_Ds[i,:]
        safe_target_Ds = safe_cluster_Ds[i,:]
        
        good_idxs = (target_Ds>0)[None,:]*(cluster_Ds>0)
    
        dfs = good_idxs.sum(axis=1)
    
        mses = numpy.square(target_fs[None,:]-cluster_fs)*good_idxs
        mse_primes = numpy.square(1-target_fs[None,:]-cluster_fs)*good_idxs
    
    
        # New version
        total_inverse_D = 1.0/safe_target_Ds[None,:]+1.0/safe_cluster_Ds
        avg_fs = (target_As[None,:]+cluster_As)*1.0/(safe_target_Ds[None,:]+safe_cluster_Ds)
        avg_f_primes = (target_Ds[None,:]-target_As[None,:]+cluster_As)*1.0/(safe_target_Ds[None,:]+safe_cluster_Ds)
        
        variances = (avg_fs*(1-avg_fs)*total_inverse_D )*good_idxs
        variance_primes = (avg_f_primes*(1-avg_f_primes)*total_inverse_D )*good_idxs
        
        # Old version
        #avg_D = (target_D[None,:] + cluster_Ds)/2.0
        #avg_f = (target_fs[None,:] + cluster_fs)/2.0
        #avg_f_prime = (1-target_fs[None,:] + cluster_fs)/2.0
    
        #variances = ((avg_f*(1-avg_f)/(avg_D+(avg_D==0)))*good_idxs)
        #variance_primes = ((avg_f_prime*(1-avg_f_prime)/(avg_D+(avg_D==0)))*good_idxs)
    
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
def cluster_snps_by_distance(cluster_As, cluster_Ds, max_num_snps_to_cluster=2000, max_d=config.cluster_distance_threshold_reads,cluster_label_offset=0,min_coverage=10):

    good_idxs = (cluster_Ds>=min_coverage)
    cluster_As = cluster_As*good_idxs
    cluster_Ds = cluster_Ds*good_idxs

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

        sys.stderr.write("Forming flat clusters at distance=%g...\n" % max_d)
        nodes = fcluster(Z, max_d, criterion="distance")
        sys.stderr.write("Done!\n")
    # Now figure out polarizations and centroids    
        
    cluster_snp_map = {}
    for snp_idx in xrange(0,len(nodes)):
    
        cluster_label = nodes[snp_idx]+cluster_label_offset
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
def cluster_secondary_snps_by_distance(cluster_map, cluster_As, cluster_Ds, max_num_snps_to_cluster=2000, max_d=config.cluster_distance_threshold_reads,min_coverage=10):

    if len(cluster_As)==0:
        return {}

    good_idxs = (cluster_Ds>=min_coverage)
    cluster_As = cluster_As*good_idxs
    cluster_Ds = cluster_Ds*good_idxs

    cluster_labels = []
    cluster_avg_fs = []
    cluster_avg_Ds = []
    
    for cluster_label in sorted(cluster_map):
        cluster_labels.append(cluster_label)
        cluster_avg_fs.append(cluster_map[cluster_label]['centroid'][0])
        cluster_avg_Ds.append(cluster_map[cluster_label]['centroid'][1])
        
    cluster_labels = numpy.array(cluster_labels)
    cluster_avg_fs = numpy.array(cluster_avg_fs)
    cluster_avg_Ds = numpy.array(cluster_avg_Ds)
    
    cluster_fs = cluster_As*1.0/(cluster_Ds+(cluster_Ds==0))

    sys.stderr.write("Calculating distances for %d snps...\n" % len(cluster_As))
    
     
    distance_matrix, distance_matrix_1, distance_matrix_2 = calculate_secondary_distance_matrix(cluster_avg_fs, cluster_avg_Ds, cluster_As, cluster_Ds) 
    
    flipped_matrix = (distance_matrix_2<=distance_matrix_1)  

    
    # Loop through and find closest cluster (& polarization) for each SNV
    best_idxs = distance_matrix.argmin(axis=0)
    best_distances = distance_matrix.min(axis=0)
    new_cluster_map = {}
    
    leftover_cluster_As = []
    leftover_cluster_Ds = []
    idx_leftover_idx_map = {}
    leftover_idx_idx_map = {}
        
    for i in xrange(0,len(cluster_As)):
        best_distance = best_distances[i]
        best_idx = best_idxs[i]
        best_cluster_label = cluster_labels[best_idx]
        
        if best_distance > max_d:
            
            idx_leftover_idx_map[i] = len(leftover_cluster_As)   
            leftover_idx_idx_map[len(leftover_cluster_As)] = i
            leftover_cluster_As.append( cluster_As[i] )
            leftover_cluster_Ds.append( cluster_Ds[i] )
        
        else:    
        
            if best_cluster_label not in new_cluster_map:
                new_cluster_map[best_cluster_label] = {}
                new_cluster_map[best_cluster_label]['centroid'] = cluster_map[best_cluster_label]['centroid']
                new_cluster_map[best_cluster_label]['snps'] = []
       
          
            flip = flipped_matrix[best_idx,i]
            new_cluster_map[best_cluster_label]['snps'].append((i,flip))
    
    leftover_cluster_As = numpy.array(leftover_cluster_As)
    leftover_cluster_Ds = numpy.array(leftover_cluster_Ds)
    
    # Now do leftovers
    leftover_cluster_map = cluster_snps_by_distance(leftover_cluster_As, leftover_cluster_Ds,max_d=max_d,cluster_label_offset=len(cluster_labels),min_coverage=min_coverage)
    for cluster_label in sorted(leftover_cluster_map):
        
        new_cluster_map[cluster_label] = {}
        new_cluster_map[cluster_label]['centroid'] = leftover_cluster_map[cluster_label]['centroid']
        new_cluster_map[cluster_label]['snps'] = []
        
        for leftover_snp_idx, flip in leftover_cluster_map[cluster_label]['snps']:
            
            new_cluster_map[cluster_label]['snps'].append( (leftover_idx_idx_map[leftover_snp_idx] , flip) )
                
    return new_cluster_map
    
    
def fast_cluster_snps_by_distance(cluster_As, cluster_Ds, max_num_snps_to_cluster=1000, max_d=config.cluster_distance_threshold_reads,min_coverage=10):
    
    
    num_snvs = len(cluster_As)
    
    
    if num_snvs < max_num_snps_to_cluster:
        return cluster_snps_by_distance(cluster_As, cluster_Ds, max_d=max_d,min_coverage=min_coverage)
    
    good_idxs = (cluster_Ds>=min_coverage)
    cluster_As = cluster_As*good_idxs
    cluster_Ds = cluster_Ds*good_idxs

    seed_cluster_As = []
    seed_cluster_Ds = []
    seed_idx_idx_map = {}
    
    rest_cluster_As = []
    rest_cluster_Ds = []
    rest_idx_idx_map = {}
    
    p = max_num_snps_to_cluster*1.0/num_snvs
        
    for snp_idx in xrange(0,num_snvs):
            
        if random() < p:
            seed_idx_idx_map[len(seed_cluster_As)] = snp_idx
            seed_cluster_As.append(cluster_As[snp_idx])
            seed_cluster_Ds.append(cluster_Ds[snp_idx])
        else:
            rest_idx_idx_map[len(rest_cluster_As)] = snp_idx
            rest_cluster_As.append(cluster_As[snp_idx])
            rest_cluster_Ds.append(cluster_Ds[snp_idx])
                
    seed_cluster_As = numpy.array(seed_cluster_As)
    seed_cluster_Ds = numpy.array(seed_cluster_Ds)
    rest_cluster_As = numpy.array(rest_cluster_As)
    rest_cluster_Ds = numpy.array(rest_cluster_Ds)
    
    
    seed_cluster_map = cluster_snps_by_distance(seed_cluster_As, seed_cluster_Ds,max_d=max_d)
    rest_cluster_map = cluster_secondary_snps_by_distance(seed_cluster_map, rest_cluster_As, rest_cluster_Ds, max_d=(max_d))
    new_cluster_map = {}
    for cluster_label in sorted(seed_cluster_map):
        
        if cluster_label not in new_cluster_map:
            new_cluster_map[cluster_label] = {}
            new_cluster_map[cluster_label]['centroid'] = seed_cluster_map[cluster_label]['centroid']
            new_cluster_map[cluster_label]['snps'] = []
        
        for seed_snp_idx, flip in seed_cluster_map[cluster_label]['snps']:
            
            new_cluster_map[cluster_label]['snps'].append( (seed_idx_idx_map[seed_snp_idx] , flip) )
    
    for cluster_label in sorted(rest_cluster_map):
        
        if cluster_label not in new_cluster_map:
            new_cluster_map[cluster_label] = {}
            new_cluster_map[cluster_label]['centroid'] = rest_cluster_map[cluster_label]['centroid']
            new_cluster_map[cluster_label]['snps'] = []
        
        for rest_snp_idx, flip in rest_cluster_map[cluster_label]['snps']:
            
            new_cluster_map[cluster_label]['snps'].append( (rest_idx_idx_map[rest_snp_idx] , flip) )
                  
    return new_cluster_map
    
def cluster_clusters_by_distance(distance_matrix, max_d=0.04):

    Y = squareform(distance_matrix)
    Z =  linkage(Y, method='average')
    nodes = fcluster(Z, max_d, criterion="distance")
        
    cluster_snp_map = {}
    for snp_idx in xrange(0,len(nodes)):
    
        cluster_label = nodes[snp_idx]
        if cluster_label not in cluster_snp_map:
            cluster_snp_map[cluster_label] = []
        cluster_snp_map[cluster_label].append(snp_idx)
   
    return cluster_snp_map
  
    

                
        
    
    