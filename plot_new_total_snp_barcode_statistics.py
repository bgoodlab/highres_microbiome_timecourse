import sys
import cPickle
import numpy
from math import fabs,log10
import parse_midas_data
import parse_timecourse_data
from numpy.random import random_sample
from scipy.stats import fisher_exact, chi2_contingency

#####
#
# Load trajectories 
#
#####
# create output filename
pickle_filename = sys.argv[1]
species_name_items = (pickle_filename.split("/")[-1]).split(".")[0].split("_")
species_name = "_".join(species_name_items[0:3])
print species_name
debug=False
chunk_size=100000000
location_trajectory_idx_map = {}

approx_genome_length = 5249241


snp_samples = parse_timecourse_data.morteza_samples
sample_size = len(snp_samples)

snp_map = {} # keys = (contig, location) 
             # values = (As,Ds) (filtered a little bit)

sys.stderr.write("Loading SNPs...\n") 

As = []
Ds = []
polarizations = []
location_trajectory_idx_map = {}

barcodes_instead_of_reads=True
if barcodes_instead_of_reads:
    snps_directory=parse_midas_data.barcode_directory+"/barcode_snps/"
else:
    snps_directory=parse_midas_data.default_snps_directory


final_line_number = 0
while final_line_number >= 0:
    
    sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
    dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_samples=snp_samples, chunk_size=chunk_size,initial_line_number=final_line_number,snps_directory=snps_directory)
    sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
    snp_samples = dummy_samples
        
    for gene_name in allele_counts_map.keys():    
        for var_type in allele_counts_map[gene_name].keys():
            chunk_polarizations = allele_counts_map[gene_name][var_type]['polarizations']
            locations = allele_counts_map[gene_name][var_type]['locations']
            allele_counts = allele_counts_map[gene_name][var_type]['alleles']
            if len(allele_counts)==0:
                continue
                    
            depths = allele_counts.sum(axis=2)
            depths = depths*(depths>=10)
            alts = allele_counts[:,:,0]
            
            freqs = alts*1.0/(depths+(depths==0))
            mafs = numpy.fmin(freqs,1-freqs)
                
            for snp_idx in xrange(0,len(locations)):
                
                if (depths[snp_idx,:]>10).sum()<10:
                    continue
                
                if (mafs[snp_idx,:]>0.2).sum() < 1.5:
                    continue
                
                
                trajectory_idx = len(As)
                location_trajectory_idx_map[locations[snp_idx]] = trajectory_idx
                As.append(alts[snp_idx,:])
                Ds.append(depths[snp_idx])
                polarizations.append(chunk_polarizations[snp_idx])
                        
sys.stderr.write("Done loading %d SNPs!\n" % len(As)) 

As = numpy.array(As)
Ds = numpy.array(Ds)

sys.stderr.write("Calculating distance matrix...\n")
import cluster_utils
#trajectory_distance_matrix, dummy_1, dummy_2 = cluster_utils.calculate_distance_matrix(As,Ds)
#del dummy_1
#del dummy_2

trajectory_distance_matrix = numpy.ones((len(As),len(As)))


sys.stderr.write("Done!\n")

#####
#
# Load precalculated fixation events
#
####

distance_bins = numpy.logspace(2,5,40)
distance_bins = numpy.hstack([[0],distance_bins,[1e06,1e07]])

#distance_bins = numpy.array([0,100]+[500*i for i in xrange(1,21)]+[2e04,1e05,1e07])
distances = distance_bins[1:]
distance_counts = numpy.zeros_like(distances)*1.0
distance_fractions = numpy.zeros_like(distance_counts)
distance_barcodes = numpy.zeros_like(distance_counts)
all_distance_counts = numpy.zeros_like(distances)*1.0


num_gamete_id_pair_map = {1: set(), 2: set(), 3: set(), 4: set()}

slightly_unbalanced_pairs = set()
extremely_unbalanced_pairs = set()

fourgamete_ids = {} # idxs that are involved in a fourgamete pair, # times observed

blacklisted_longsnp_ids = set()

sys.stderr.write("Loading pickle from disk...\n")
d = cPickle.load( open( pickle_filename,"rb") )
sys.stderr.write("Done!\n")

id_longsnp_map = d['id_longsnp_map']
allele_idx_map = d['allele_idx_map']
B_matrix = d['B']
S_matrix = d['S']

print S_matrix[0]

gamete_idx_map = {('A','A'): (allele_idx_map['A'], allele_idx_map['A']),  ('A','R'): (allele_idx_map['A'], allele_idx_map['R']), ('R','A'): (allele_idx_map['R'], allele_idx_map['A']), ('R','R'): (allele_idx_map['R'], allele_idx_map['R'])}

B_vector = numpy.array([B_matrix[longsnp_id].sum(axis=1) for longsnp_id in xrange(0,len(id_longsnp_map))])


# First loop through and find SNPs that themselves fail 2 gamete test
for focal_longsnp_id in xrange(0,len(id_longsnp_map)):
      
    is_bad_site = True
      
    # Can't look if there are no linked SNVs
    if (id_longsnp_map[focal_longsnp_id] in location_trajectory_idx_map)     and (len(S_matrix[focal_longsnp_id])>0):
    
        if focal_longsnp_id in S_matrix[focal_longsnp_id]:
        
            self_gametes = S_matrix[focal_longsnp_id][focal_longsnp_id]
    
            self_barcodes = self_gametes.sum()
    
            self_bad_barcodes = self_barcodes - self_gametes[gamete_idx_map[('A','A')]] - self_gametes[gamete_idx_map[('R','R')]]
    
            if self_bad_barcodes<0.5:
                # a good SNP!
                is_bad_site = False
                
    if is_bad_site:
        blacklisted_longsnp_ids.add(focal_longsnp_id)    
        
for focal_longsnp_id in xrange(0,len(id_longsnp_map)):
      
    if focal_longsnp_id in blacklisted_longsnp_ids:
         continue
         
    species_name, contig, location = id_longsnp_map[focal_longsnp_id]   
    
    focal_trajectory_idx = location_trajectory_idx_map[id_longsnp_map[focal_longsnp_id]]  
    
    alleles = B_vector[focal_longsnp_id]
    total_alleles = alleles.sum() 
    
    alt_alleles = alleles[allele_idx_map['A']]
    ref_alleles = alleles[allele_idx_map['R']]   
    
    for other_longsnp_id in S_matrix[focal_longsnp_id]:
    
        if other_longsnp_id==focal_longsnp_id:
            continue # don't look at mapping to itself
        
        if other_longsnp_id in blacklisted_longsnp_ids:
            continue # don't look at bad ones!
        
        other_species, other_contig, other_location = id_longsnp_map[other_longsnp_id]
        
        if location<other_location:
            L1 = location
            L2 = other_location
        else:
            L1 = other_location
            L2 = location
            
        
            
        distance = L2-L1
        distance = min([distance, approx_genome_length-distance])
        
        #distance = fabs(other_location-location) 
            
        gametes = S_matrix[focal_longsnp_id][other_longsnp_id][0]
    
        total_barcodes = gametes.sum()
        
        fraction = total_barcodes*1.0/total_alleles
        
        distance_idx = numpy.digitize([distance],distance_bins)[0]-1
        all_distance_counts[distance_idx] += 1
         
        if total_barcodes>2.5: # hopefully cut down on some errors...
        
            # get distance bin
            distance_counts[distance_idx]+=1
            distance_fractions[distance_idx]+=fraction
            distance_barcodes[distance_idx]+= total_barcodes
        
        if distance <0:
            continue
        
        gamete_freqs = gametes*1.0/total_barcodes
        # minimum number of barcodes
        if total_barcodes < 10:
            continue
        
        pair_alt_alleles = gametes[gamete_idx_map[('A','A')]] +  gametes[gamete_idx_map[('A','R')]]
        pair_ref_alleles = gametes[gamete_idx_map[('R','A')]] +  gametes[gamete_idx_map[('R','R')]]
        
        expected_pair_alt_alleles = total_barcodes*(alt_alleles*1.0/total_alleles)
        
        expected_pair_ref_alleles = total_barcodes*(ref_alleles*1.0/total_alleles)
        
        variance = expected_pair_alt_alleles*expected_pair_ref_alleles/total_barcodes
        
        scaled_deviation = ((pair_alt_alleles-expected_pair_alt_alleles)**2/variance)**0.5
        
        if scaled_deviation>2:
            slightly_unbalanced_pairs.add((focal_longsnp_id, other_longsnp_id))
            slightly_unbalanced_pairs.add((other_longsnp_id, focal_longsnp_id))
        
        if (expected_pair_alt_alleles > 20 and pair_alt_alleles==0) or (expected_pair_ref_alleles > 20 and pair_ref_alleles==0):
            # Extremely unbalanced alleles
            extremely_unbalanced_pairs.add((focal_longsnp_id, other_longsnp_id))
            extremely_unbalanced_pairs.add((other_longsnp_id, focal_longsnp_id))
            
        gamete_idxs = (gametes>0.5)
        num_gametes = gamete_idxs.sum()
        stringent_gamete_idxs = ((gametes>2.5)*(gamete_freqs>0.05))
        stringent_num_gametes = stringent_gamete_idxs.sum()
    
        if stringent_num_gametes==4:
            # fails the 4 gamete test!
            num_gamete_id_pair_map[4].add((focal_longsnp_id, other_longsnp_id))
            num_gamete_id_pair_map[4].add((other_longsnp_id, focal_longsnp_id))
            
            if focal_longsnp_id not in fourgamete_ids:
                fourgamete_ids[focal_longsnp_id] = 0
            if other_longsnp_id not in fourgamete_ids:
                fourgamete_ids[other_longsnp_id] = 0
                
            fourgamete_ids[focal_longsnp_id] += 1
            fourgamete_ids[other_longsnp_id] += 1
            
        elif num_gametes >= 2.5:
            # not perfectly linked!
            
            num_gamete_id_pair_map[3].add((focal_longsnp_id, other_longsnp_id))
            num_gamete_id_pair_map[3].add((other_longsnp_id, focal_longsnp_id))
            
                        
        elif num_gametes==2:
        
            # First check whether there are bad combinations of gametes
            if not ((gamete_idxs[gamete_idx_map[('A','A')]] and gamete_idxs[gamete_idx_map[('R','R')]]) or (gamete_idxs[gamete_idx_map[('A','R')]] and gamete_idxs[gamete_idx_map[('R','A')]])):
                # bad!
                num_gamete_id_pair_map[3].add((focal_longsnp_id, other_longsnp_id))
                num_gamete_id_pair_map[3].add((other_longsnp_id, focal_longsnp_id))
            
            else:
                
                if stringent_num_gametes==2:
                    num_gamete_id_pair_map[2].add((focal_longsnp_id, other_longsnp_id))
                    num_gamete_id_pair_map[2].add((other_longsnp_id, focal_longsnp_id))
            
                else:
                    num_gamete_id_pair_map[1].add((focal_longsnp_id, other_longsnp_id))
                    num_gamete_id_pair_map[1].add((other_longsnp_id, focal_longsnp_id))
            
        
        else: # (num_gametes==1)
            
            num_gamete_id_pair_map[1].add((focal_longsnp_id, other_longsnp_id))
            num_gamete_id_pair_map[1].add((other_longsnp_id, focal_longsnp_id))
  
# make this later
trajectory_distances = {1: [], 1.1: [], 2: [], 3: [], 4: [], 2.1: []}
for g in num_gamete_id_pair_map:
    
    for id1,id2 in num_gamete_id_pair_map[g]:
        
        if g < 3.5 and ((id1 in fourgamete_ids) or (id2 in fourgamete_ids)):
            continue
        
        # Get trajectory distance
        
        idx1 = location_trajectory_idx_map[id_longsnp_map[id1]]
        idx2 = location_trajectory_idx_map[id_longsnp_map[id2]]
        trajectory_distance = trajectory_distance_matrix[idx1, idx2]
        
        if g<2.5 and (id1,id2) in extremely_unbalanced_pairs:
            effective_g = g+0.1
        else:
            effective_g = g
        
        trajectory_distances[effective_g].append(trajectory_distance)
        
           
distance_fractions = distance_fractions*1.0/(
distance_counts+(distance_counts==0))

distance_barcodes = distance_barcodes*1.0/(
distance_counts+(distance_counts==0))

fourgamete_pairs = num_gamete_id_pair_map[4]
        
record_str_items = []
for id in sorted(fourgamete_ids):
    longsnp = id_longsnp_map[id]
    if fourgamete_ids[id]>2.5:
        record_str_items.append("%s|%d|%d" % (longsnp[0],longsnp[1],fourgamete_ids[id]))

print " ".join(record_str_items)
print ""
    
print "Distance stats:"
for idx in xrange(0,len(distances)):
    if distance_counts[idx]<0.5:
        continue
    
    print distances[idx], all_distance_counts[idx], distance_counts[idx], distance_fractions[idx], distance_barcodes[idx]

print "Gamete counts:"
min_distance = 0.3
max_distance = min_distance*2
for g in sorted(trajectory_distances):
    if len(trajectory_distances[g])==0:
        continue
        
    trajectory_distances[g] = numpy.array(trajectory_distances[g])
    
    max_distance = max([max_distance, trajectory_distances[g].max()])
    
    print g, len(trajectory_distances[g]), trajectory_distances[g].mean(), numpy.median(trajectory_distances[g])
    
# Time to plot distance survival function 
import pylab
import matplotlib as mpl
import stats_utils

mpl.rcParams['font.size'] = 7
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'


pylab.figure(1,figsize=(3.42,2))
pylab.semilogx([min_distance, max_distance],[0,0],'k:')

bins = numpy.logspace(log10(min_distance),log10(max_distance),50)

pylab.hist(trajectory_distances[1],bins=bins, histtype=u'step', label=('1 gamete (n=%d)'  % len(trajectory_distances[1])))
pylab.hist(trajectory_distances[2],bins=bins, histtype=u'step', label=('2 gametes (n=%d)' % len(trajectory_distances[2])))
pylab.hist(trajectory_distances[3],bins=bins, histtype=u'step', label=('3 gametes (n=%d)' % len(trajectory_distances[3])))
pylab.hist(trajectory_distances[4],bins=bins, histtype=u'step', label=('4 gametes (n=%d)' % len(trajectory_distances[4])))


pylab.xlabel('Trajectory distance, $d$')
pylab.ylabel('Number of SNV pairs')
pylab.yticks([])
pylab.legend(loc='upper left',frameon=False)
pylab.savefig('%s_num_gametes_vs_trajectory_distance.pdf' % (species_name),bbox_inches='tight')

pylab.figure(2,figsize=(3.42,2))
pylab.xlabel('Coordinate distance, $\ell$ (bp)')
pylab.ylabel('Fraction barcodes shared')
pylab.ylim([5e-04,1])
pylab.loglog(distances[distance_counts>0.5], distance_fractions[distance_counts>0.5],'k.-')
pylab.savefig('%s_shared_barcodes_vs_distance.pdf' % (species_name), bbox_inches='tight')

# When I make the plot 

# Distance vs average # of barcodes 