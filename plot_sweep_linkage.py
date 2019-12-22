import sys
import cPickle
import numpy
from math import fabs,log10
import parse_midas_data
import parse_timecourse_data
from numpy.random import random_sample
from scipy.stats import fisher_exact, chi2_contingency
import os.path

#####
#
# Load trajectories 
#
#####
# create output filename
pickle_directory = "/Users/bgood/highres_microbiome_timecourse_barcode_data/barcode_trajectories"
species_name = None
debug=False
chunk_size=100000000

features_filename = sys.argv[1]

print "Created by plot_sweep_linkage.py"
print "from", features_filename

features_file = open(features_filename,"r")
# each line contains a feature. starts with species name, then TSL of alleles
species_allele_map = {}
for line in features_file:
    
    # Ignore comments
    if line[0]=='#': 
        continue 
        
    items = line.split()
    species_name = items[0]
    
    if True:
        item = items[1]
        allele_name_prefix = item.strip()
        subitems = allele_name_prefix.split("|")
        allele = (subitems[0],long(subitems[1]))
        if species_name not in species_allele_map:
            species_allele_map[species_name] = set()
        species_allele_map[species_name].add(allele)
        
features_file.close()


species_names = sorted(species_allele_map.keys())

species_trajectory_distances = {}

for species_name in species_names:


    pickle_filename = "%s/%s.total.p" % (pickle_directory, species_name)
    if not os.path.isfile(pickle_filename):
        continue
    
    sys.stderr.write("Processing %s...\n" % species_name)    
    print species_name
    min_distance = 0

    distance_bins = numpy.array([0,100]+[500*i for i in xrange(1,21)]+[2e04,1e05,1e07])
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
    gamete_idx_map = d['gamete_idx_map']
    d = d['snp_barcode_timecourse']

    # First loop through and find SNPs that themselves fail 2 gamete test
    whitelisted_longsnp_ids = set()
    for focal_longsnp_id in d:
      
        if id_longsnp_map[focal_longsnp_id] in species_allele_map[species_name]:
            whitelisted_longsnp_ids.add(focal_longsnp_id)
      
    for focal_longsnp_id in sorted(d):
      
        if focal_longsnp_id not in whitelisted_longsnp_ids:
            continue
         
        location = id_longsnp_map[focal_longsnp_id][1]    
    
        alleles = d[focal_longsnp_id]['all']
        
        total_alleles = alleles.sum() 
    
        alt_alleles = alleles[allele_idx_map['A']]
        ref_alleles = alleles[allele_idx_map['R']]   
    
        for other_longsnp_id in sorted(d[focal_longsnp_id]['longsnps']):
        
            if other_longsnp_id==focal_longsnp_id:
                linkage_str = ("%s %s|%d <->self" % (species_name, id_longsnp_map[focal_longsnp_id][0], id_longsnp_map[focal_longsnp_id][1]) )
            
            else:
                linkage_str = ("%s %s|%d <-> %s|%d" % (species_name, id_longsnp_map[focal_longsnp_id][0], id_longsnp_map[focal_longsnp_id][1], id_longsnp_map[other_longsnp_id][0], id_longsnp_map[other_longsnp_id][1]) )
            
            #    continue # don't look at mapping to itself
        
            if other_longsnp_id not in whitelisted_longsnp_ids:
                continue # only look at linkage within sweep
        
            other_location = id_longsnp_map[other_longsnp_id][1]
            distance = fabs(other_location-location) 
            
            gametes = d[focal_longsnp_id]['longsnps'][other_longsnp_id]
            total_barcodes = gametes.sum()
        
            if total_barcodes<10:
                continue
        
            
            print linkage_str
            
            if gametes.min()>0.5:
                print "***"
        
            
            print "AA:%d, RR:%d, AR:%d, RA:%d" % (gametes[gamete_idx_map[('A','A')]], gametes[gamete_idx_map[('R','R')]], gametes[gamete_idx_map[('A','R')]], gametes[gamete_idx_map[('R','A')]])
            
            

    