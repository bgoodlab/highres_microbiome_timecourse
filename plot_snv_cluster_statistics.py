import numpy
import pylab
import sys
import cluster_utils
from numpy.random import shuffle

filename = sys.argv[1]
file = open(filename,"r")

self_fs = []
non_self_fs = []

cluster_idx_map = {}
idx_cluster_map = []

cluster_distance_map = {}

file.readline() # header
for line in file:
	items = line.split(",")
	
	cluster_1 = items[0].strip()
	cluster_2 = items[1].strip()
	
	if cluster_1 not in cluster_idx_map:
		cluster_idx_map[cluster_1] = len(idx_cluster_map)
		idx_cluster_map.append(cluster_1)
		
	if cluster_2 not in cluster_idx_map:
		cluster_idx_map[cluster_2] = len(idx_cluster_map)
		idx_cluster_map.append(cluster_2)
		
	
	Bs = [float(item) for item in items[2:6]]
	Btot = sum(Bs)
	third_gamete_fraction = float(items[-1])
	
	if Btot<100:
		continue
	
	if cluster_1 not in cluster_distance_map:
		cluster_distance_map[cluster_1] = {}
	
	if cluster_2 not in cluster_distance_map[cluster_1]:
		cluster_distance_map[cluster_1][cluster_2] = third_gamete_fraction
	
	if cluster_1 == cluster_2:
		
		self_fs.append(third_gamete_fraction)
		
	else:
		
		non_self_fs.append(third_gamete_fraction)
		
print sorted(self_fs)
print sorted(non_self_fs)

# Now form greedy clusters: 
used_clusters = set()

sorted_clusters = []

for cluster_1 in idx_cluster_map:
	if cluster_1 not in sorted_clusters:
		sorted_clusters.append(cluster_1)
	
	for cluster_2 in cluster_distance_map[cluster_1]:
		if cluster_distance_map[cluster_1][cluster_2] < 0.1:
			if cluster_2 not in sorted_clusters:
				sorted_clusters.append(cluster_2)
				
		
cluster_distance_matrix = numpy.zeros((len(sorted_clusters), len(sorted_clusters)))
#shuffle(sorted_clusters)

for idx_1 in xrange(0,len(sorted_clusters)):
	cluster_1 = sorted_clusters[idx_1]
	for idx_2 in xrange(0,len(sorted_clusters)):
		cluster_2 = sorted_clusters[idx_2]
		
		if cluster_1 in cluster_distance_map:
			if cluster_2 in cluster_distance_map[cluster_1]:
				cluster_distance_matrix[idx_1,idx_2] = cluster_distance_map[cluster_1][cluster_2]
				cluster_distance_matrix[idx_2,idx_1] = cluster_distance_map[cluster_1][cluster_2]
				
		if cluster_2 in cluster_distance_map:
			if cluster_1 in cluster_distance_map[cluster_2]:
				cluster_distance_matrix[idx_1,idx_2] = cluster_distance_map[cluster_2][cluster_1]
				cluster_distance_matrix[idx_2,idx_1] = cluster_distance_map[cluster_2][cluster_1]
		
	cluster_distance_matrix[idx_1,idx_1] = 0
	
supercluster_cluster_map = cluster_utils.cluster_clusters_by_distance(cluster_distance_matrix, max_d=0.04)

for supercluster in sorted(supercluster_cluster_map):
	print [sorted_clusters[idx] for idx in supercluster_cluster_map[supercluster]]
				
pylab.matshow(cluster_distance_matrix)
pylab.show()

	