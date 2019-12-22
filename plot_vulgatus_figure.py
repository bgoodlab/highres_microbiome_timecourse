###############################
#
# Rest of script begins here
#
################################
import matplotlib  
matplotlib.use('Agg') 
import pylab
import numpy
import sys
from math import log10
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from numpy.random import binomial, random_sample
import bz2
import parse_midas_data
import parse_timecourse_data
import matplotlib
import matplotlib.pyplot as plt
import timecourse_utils
import parse_patric
import calculate_preexisting_snps
import cluster_utils
import os
import cPickle
import config
from math import fabs
mpl.rcParams['font.size'] = 7
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

min_coverage = 10

species_name = species_name = "Bacteroides_vulgatus_57955"
output_directory = os.path.expanduser("~/strainfinder_output/")
filename_prefix = "%s/%s" % (output_directory, species_name)

sample_time_map = parse_timecourse_data.parse_sample_time_map()
ts = numpy.array([sample_time_map[sample] for sample in parse_timecourse_data.morteza_samples])

# Load location object
snp_locations = cPickle.load(open(filename_prefix+".strainfinder.locations.p",'rb'))
snp_alignment = cPickle.load(open(filename_prefix+".strainfinder.p",'rb'))

cluster_filename = filename_prefix+".clusters.p"

cluster_As = snp_alignment[:,:,0].T
cluster_Ds = (snp_alignment.sum(axis=2)).T

#cluster_As = cluster_As[0:1000,:]
#cluster_Ds = cluster_Ds[0:1000,:]

pylab.figure(figsize=(7,4))
fig = pylab.gcf()

outer_grid  = gridspec.GridSpec(3,2, height_ratios=[1,1,1], width_ratios=[2,1], hspace=0.15)

observed_axis = plt.Subplot(fig, outer_grid[0,0])
fig.add_subplot(observed_axis)
observed_axis.set_ylim([0,1])
observed_axis.set_xlim([-5,160])
observed_axis.set_ylabel('Allele\nfrequency')
observed_axis.spines['top'].set_visible(False)
observed_axis.spines['right'].set_visible(False)
observed_axis.get_xaxis().tick_bottom()
observed_axis.get_yaxis().tick_left()

for snp_idx in xrange(0,cluster_As.shape[0]):
    
    alts = cluster_As[snp_idx,:]
    depths = cluster_Ds[snp_idx,:]
    freqs = alts*1.0/(depths+(depths==0))
    
    good_idxs = (depths>min_coverage)
    
    observed_axis.plot(ts[good_idxs],freqs[good_idxs],'-',alpha=0.5,linewidth=0.25)

    
predicted_axis = plt.Subplot(fig, outer_grid[1,0])
fig.add_subplot(predicted_axis)
predicted_axis.set_ylim([0,1])
predicted_axis.set_xlim([-5,160])

predicted_axis.set_ylabel('Allele\nfrequency')
predicted_axis.spines['top'].set_visible(False)
predicted_axis.spines['right'].set_visible(False)
predicted_axis.get_xaxis().tick_bottom()
predicted_axis.get_yaxis().tick_left()

rest_axis = plt.Subplot(fig, outer_grid[2,0])
fig.add_subplot(rest_axis)
rest_axis.set_ylim([0,1])
rest_axis.set_xlim([-5,160])

rest_axis.set_ylabel('Allele\nfrequency')
rest_axis.spines['top'].set_visible(False)
rest_axis.spines['right'].set_visible(False)
rest_axis.get_xaxis().tick_bottom()
rest_axis.get_yaxis().tick_left()
rest_axis.set_xlabel('Time (days)')

distance_axis = plt.Subplot(fig, outer_grid[1,1])
fig.add_subplot(distance_axis)
sm = cmx.ScalarMappable(norm=colors.Normalize(vmin=0, vmax=0.4),cmap=cmx.get_cmap(name='Greys'))
sm.set_array([])
fig.colorbar(sm, ax=distance_axis)


def hex_to_rgb(value):
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))

def fraction_to_rgb(f,fmax=0.5):
	c = long((1-f/fmax)*255)
	return (c,c,c)
              
#cluster_map = cluster_utils.cluster_snps_by_distance(cluster_As, cluster_Ds,max_d=config.cluster_distance_threshold_barcodes)

#cluster_map = cluster_utils.cluster_snps_by_distance(cluster_As, cluster_Ds,max_d=config.cluster_distance_threshold_barcodes)
cluster_map = cluster_utils.fast_cluster_snps_by_distance(cluster_As, cluster_Ds,max_d=config.cluster_distance_threshold_barcodes)

total_clustered_snvs = 0
total_fractional_clustered_snvs = 0
        
cluster_data = []

cluster_color_map = {}  
good_clusters = []         
cluster_ids = cluster_map.keys()
cluster_sizes = [len(cluster_map[cluster_id]['snps']) for cluster_id in cluster_ids]

sorted_cluster_sizes, sorted_cluster_ids = zip(*sorted(zip(cluster_sizes, cluster_ids), reverse=True))

output_strs = []
output_str = ", ".join(['Cluster', 'Contig', 'Location'])
output_strs.append(output_str)

cluster_freq_map = {}
cluster_color_map = {}

for cluster_id in sorted_cluster_ids:
    color = None
    avg_fs, total_Ds = cluster_map[cluster_id]['centroid']
            
    good_idxs = (total_Ds>0)
            
    if good_idxs.sum()<1.5:
        continue
    
    masked_times = ts[good_idxs]
    masked_avg_fs = avg_fs[good_idxs]
            
    cluster_size = len(cluster_map[cluster_id]['snps'])
    fractional_size = cluster_size*1.0/len(cluster_As)     
    
    sys.stderr.write("Cluster %d (n=%d SNVs, %g)\n" % (cluster_id, cluster_size, fractional_size)) 
    
    if cluster_size>=100:
    #if fractional_size>0.01:
        keeper=True
    else:
        keeper=False
        for snp_idx, flip in sorted(cluster_map[cluster_id]['snps']):
            good_idxs = cluster_Ds[snp_idx,:]>0
            masked_times = ts[good_idxs]
            masked_As = cluster_As[snp_idx,good_idxs]
            masked_Ds = cluster_Ds[snp_idx,good_idxs]
            if flip:
                masked_As = masked_Ds-masked_As
     		masked_freqs = masked_As*1.0/masked_Ds
     		predicted_axis.plot(masked_times, masked_freqs, '-', color='0.7', alpha=0.5, markersize=2, markeredgecolor='none', linewidth=0.25,zorder=1)
 		
 	if not keeper: 
 		continue
 
    line, = predicted_axis.plot([1],[-1],'.')
    color = pylab.getp(line,"color")
    cluster_color_map[cluster_id] = color
    good_clusters.append(cluster_id)
    total_clustered_snvs += cluster_size
    total_fractional_clustered_snvs += fractional_size
                 
    cluster_freq_map[cluster_id] = (cluster_size,masked_times, masked_avg_fs)
    
    cluster_snp_locations = []  
    for snp_idx, flip in sorted(cluster_map[cluster_id]['snps']):
        good_idxs = cluster_Ds[snp_idx,:]>0
        masked_times = ts[:]
        masked_As = cluster_As[snp_idx,:]
        masked_Ds = cluster_Ds[snp_idx,:]
        
        contig,location,allele = snp_locations[snp_idx]        
                
        if flip:
            masked_As = masked_Ds-masked_As
            
            if allele=='A':
                allele='R'
            else:
                allele='A'
        
        masked_freqs = masked_As*1.0/masked_Ds
        
        cluster_data.append( ( contig, location, allele, cluster_id, masked_As, masked_Ds ) )
    
    	cluster_snp_locations.append((contig, location))
    	
     
    cluster_snp_locations.sort()
    for contig, location in cluster_snp_locations:
        output_str = ", ".join([str(color), str(contig), str(location)])
        output_strs.append(output_str)
    	
sys.stderr.write("Total clustered SNVs: %d\n" % total_clustered_snvs)
sys.stderr.write("Total fraction clustered: %g\n" % total_fractional_clustered_snvs)     

filename = parse_midas_data.analysis_directory+('%s_cluster_figure_clusters.txt' % species_name)

file = open(filename,"w")
for output_str in output_strs:
    file.write(output_str)
    file.write("\n")
file.close()


snp_cluster_map = {}
for contig, location, allele, cluster_id, As, Ds in cluster_data:
    snp_cluster_map[(contig,location)] = (allele,cluster_id, As, Ds)
cluster_S_map = {}
    
pickle_filename = ("%s/new_barcode_trajectories/%s_snp_barcodes.p" % (config.barcode_directory, species_name))
    
sys.stderr.write("Processing %s...\n" % species_name) 

sys.stderr.write("Loading pickle from disk...\n")
d = cPickle.load( open( pickle_filename,"rb") )
sys.stderr.write("Done!\n")

id_longsnp_map = d['id_longsnp_map']
allele_idx_map = d['allele_idx_map']
B_matrix = d['B']
S_matrix = d['S']
    
opposite_allele = {'A':'R', 'R':'A'}    
    
gamete_idx_map = {('A','A'): (allele_idx_map['A'], allele_idx_map['A']),  ('A','R'): (allele_idx_map['A'], allele_idx_map['R']), ('R','A'): (allele_idx_map['R'], allele_idx_map['A']), ('R','R'): (allele_idx_map['R'], allele_idx_map['R'])}
    
all_gametes = gamete_idx_map.values()
on_diagonal_gametes = [gamete_idx_map[('A','A')], gamete_idx_map[('R','R')]]
off_diagonal_gametes = [gamete_idx_map[('A','R')], gamete_idx_map[('R','A')]]
    
twogamete_haplotypes = [on_diagonal_gametes, off_diagonal_gametes]

for focal_longsnp_id in S_matrix:

    focal_species_name, focal_contig, focal_location = id_longsnp_map[focal_longsnp_id]
    
    if (focal_contig, focal_location) not in snp_cluster_map:
        continue
    
    focal_allele, focal_cluster, dummy, dummy = snp_cluster_map[(focal_contig, focal_location)]
    
        
    for target_longsnp_id in S_matrix[focal_longsnp_id]:
            
        if target_longsnp_id == focal_longsnp_id:
            continue
            
        target_species_name, target_contig, target_location = id_longsnp_map[target_longsnp_id] 
        
        if (target_contig, target_location) not in snp_cluster_map:
            continue
        
        target_allele, target_cluster, dummy, dummy = snp_cluster_map[(target_contig, target_location)]
            
        if target_contig!=focal_contig:
            continue
                
        if fabs(target_location-focal_location)>5e04:
            continue
                    
        Ss = S_matrix[focal_longsnp_id][target_longsnp_id][0]
        S0s = S_matrix[focal_longsnp_id][target_longsnp_id][1]
        
        cluster_pair = (focal_cluster, target_cluster)
        
        if cluster_pair not in cluster_S_map:
            cluster_S_map[cluster_pair] = numpy.zeros_like(Ss)

        cluster_S_map[cluster_pair][gamete_idx_map[('A','A')]] += Ss[gamete_idx_map[(focal_allele,target_allele)]]
        cluster_S_map[cluster_pair][gamete_idx_map[('A','R')]] += Ss[gamete_idx_map[(focal_allele,opposite_allele[target_allele])]]
        cluster_S_map[cluster_pair][gamete_idx_map[('R','A')]] += Ss[gamete_idx_map[(opposite_allele[focal_allele],target_allele)]]
        cluster_S_map[cluster_pair][gamete_idx_map[('R','R')]] += Ss[gamete_idx_map[(opposite_allele[focal_allele],opposite_allele[target_allele])]]

output_strs = []

# Header
output_str = ", ".join(['Clade1','Clade2','B_ab','B_AB','B_Ab','B_aB', 'third_gamete_fraction'])
output_strs.append(output_str)
        
cluster_distance_matrix = numpy.zeros((len(good_clusters), len(good_clusters)))

for i in xrange(0,len(good_clusters)):
    cluster_1 = good_clusters[i]
    color_1 = cluster_color_map[cluster_1]
    size_1, dummy1, dummy2 = cluster_freq_map[cluster_1]
    
    for j in xrange(i,len(good_clusters)):
        cluster_2 = good_clusters[j]
        color_2 = cluster_color_map[cluster_2]
        size_2, dummy1, dummy2 = cluster_freq_map[cluster_2]
    
        orientation_1 = cluster_S_map[(cluster_1,cluster_2)][gamete_idx_map[('R','R')]]+cluster_S_map[(cluster_1,cluster_2)][gamete_idx_map[('A','A')]]
        orientation_2 = cluster_S_map[(cluster_1,cluster_2)][gamete_idx_map[('A','R')]]+ cluster_S_map[(cluster_1,cluster_2)][gamete_idx_map[('R','A')]]
        
        max_orientation = max([orientation_1, orientation_2])
        min_orientation = min([orientation_1, orientation_2])
        
        third_gamete_fraction = min_orientation * 1.0 / (min_orientation+max_orientation)
        
        if True: #cluster_1!=cluster_2:
        	cluster_distance_matrix[i,j] = third_gamete_fraction
        	cluster_distance_matrix[j,i] = third_gamete_fraction
        
        output_str = ("%s (%g), %s (%g), %s, %s, %s, %s, %s" % ( color_1, size_1, color_2, size_2, cluster_S_map[(cluster_1,cluster_2)][gamete_idx_map[('R','R')]], cluster_S_map[(cluster_1,cluster_2)][gamete_idx_map[('A','A')]], cluster_S_map[(cluster_1,cluster_2)][gamete_idx_map[('A','R')]], cluster_S_map[(cluster_1,cluster_2)][gamete_idx_map[('R','A')]], third_gamete_fraction))
        
        output_strs.append(output_str)

filename = parse_midas_data.analysis_directory+('%s_cluster_figure_barcodes.txt' % species_name)

# Now find superclusters!
supercluster_cluster_map = cluster_utils.cluster_clusters_by_distance(cluster_distance_matrix-numpy.diag(numpy.diag(cluster_distance_matrix)), max_d=0.04)

cluster_supercluster_map = {}
supercluster_color_map = {}
supercluster_As_map = {}
supercluster_Ds_map = {}

cluster_idx_list = []
for supercluster in supercluster_cluster_map:
	line, = rest_axis.plot([1],[-1],'.')
	color = pylab.getp(line,'color')
	supercluster_color_map[supercluster] = color
	supercluster_As_map[supercluster] = numpy.zeros_like(ts)
	supercluster_Ds_map[supercluster] = numpy.zeros_like(ts)
	
	for cluster_idx in supercluster_cluster_map[supercluster]:
		cluster_supercluster_map[good_clusters[cluster_idx]] = supercluster	
		cluster_idx_list.append(cluster_idx)

plotted_cluster_distance_matrix = numpy.zeros((len(good_clusters)+1, len(good_clusters)+1, 3),dtype=numpy.int32)

plotted_cluster_distance_matrix[0,0] = fraction_to_rgb(0)

for i in xrange(0,len(cluster_idx_list)):
    cluster_1 = good_clusters[cluster_idx_list[i]]
    color_1 = cluster_color_map[cluster_1]
    plotted_cluster_distance_matrix[0,i+1] = hex_to_rgb(color_1)
    plotted_cluster_distance_matrix[i+1,0] = hex_to_rgb(color_1)
    
    for j in xrange(0,len(cluster_idx_list)):
        cluster_2 = good_clusters[cluster_idx_list[j]]
        color_2 = cluster_color_map[cluster_2]
    	
        plotted_cluster_distance_matrix[i+1,j+1] = fraction_to_rgb(cluster_distance_matrix[cluster_idx_list[i],cluster_idx_list[j]],fmax=0.4)
        plotted_cluster_distance_matrix[j+1,i+1] = plotted_cluster_distance_matrix[i+1,j+1]


print plotted_cluster_distance_matrix[0,1], plotted_cluster_distance_matrix[1,1] 
pos = distance_axis.imshow(plotted_cluster_distance_matrix)
distance_axis.set_xticks([])
distance_axis.set_yticks([])
distance_axis.plot([0.5,0.5],[-0.5,plotted_cluster_distance_matrix.shape[0]-0.5],'k-')
distance_axis.plot([-0.5,plotted_cluster_distance_matrix.shape[0]-0.5],[0.5,0.5],'k-')


for contig, location, allele, cluster_id, As, Ds in cluster_data:

	# First plot regular clusters
	good_idxs = Ds>0
	masked_ts = ts[good_idxs]
	masked_As = As[good_idxs]
	masked_Ds = Ds[good_idxs]
	masked_freqs = masked_As*1.0/masked_Ds
	
	color = cluster_color_map[cluster_id]
	
	cluster_size, dummy, dummy = cluster_freq_map[cluster_id]
	
	if cluster_size>=1000:
		colored_axis = predicted_axis
		grey_axis = rest_axis
	else:
		colored_axis = rest_axis
		grey_axis = predicted_axis
	
	colored_axis.plot(masked_ts, masked_freqs, '-', color=color, alpha=0.5, markersize=2, markeredgecolor='none', linewidth=0.25,zorder=2)
	grey_axis.plot(masked_ts, masked_freqs, '-', color='0.7', alpha=0.5, markersize=2, markeredgecolor='none', linewidth=0.25,zorder=1)
	
	# Now plot supercluster
	supercluster = cluster_supercluster_map[cluster_id]
	color = supercluster_color_map[supercluster]
	
	#rest_axis.plot(masked_ts, masked_freqs, '-', color=color, alpha=0.5, markersize=2, markeredgecolor='none', linewidth=0.25)

	# Add in to supercluster
	supercluster_As_map[supercluster] += As
	supercluster_Ds_map[supercluster] += Ds
	
# Now plot cluster averaged trajectories
for cluster_id in good_clusters:

	color = cluster_color_map[cluster_id]
	cluster_size, masked_times, masked_avg_fs = cluster_freq_map[cluster_id]
	markersize=log10(cluster_size)
	
	if cluster_size>=1000:
		axis = predicted_axis
	else:
		axis = rest_axis
	line, = axis.plot(masked_times, masked_avg_fs, 'k-', alpha=1, zorder=5, linewidth=0.5)
	#line, = predicted_axis.plot(masked_times, masked_avg_fs, 'o', alpha=1, markersize=markersize, zorder=6, linewidth=0.5,color=color)

# Now plot supercluster averaged trajectories
for supercluster in supercluster_As_map:

	color = supercluster_color_map[supercluster]
	As = supercluster_As_map[supercluster]
	Ds = supercluster_Ds_map[supercluster]
	good_idxs = Ds>0
	
	masked_times = ts[good_idxs]
	masked_As = As[good_idxs]
	masked_Ds = Ds[good_idxs]
	masked_avg_fs = masked_As*1.0/masked_Ds
	#line, = rest_axis.plot(masked_times, masked_avg_fs, 'k-', alpha=1, markersize=3, zorder=5, linewidth=0.5)
	#line, = rest_axis.plot(masked_times, masked_avg_fs, 'o', alpha=1, markersize=3, zorder=6, linewidth=0.5,color=color)
    

file = open(filename,"w")
for output_str in output_strs:
    file.write(output_str)
    file.write("\n")
file.close()

filename = parse_midas_data.analysis_directory+('%s_cluster_figure.png' % species_name)
fig.savefig(filename, bbox_inches='tight', dpi=300, transparent=True)    
sys.stderr.write("Done!\n")

