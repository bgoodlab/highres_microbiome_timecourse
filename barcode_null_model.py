import sys
import numpy
import parse_midas_data
import parse_timecourse_data
import config
import cPickle as pickle
from math import ceil
import sys
import argparse
from scipy.special import gammaln as loggamma
from scipy.special import gammainc
from scipy.special import expm1
from math import log,exp
from scipy.stats import poisson

from numpy.random import shuffle
import barcode_utils

    
###########################################################################
#
# Standard header to read in argument information
#
###########################################################################


Dbins = numpy.logspace(1,4,16)-1e-06
Ds = numpy.array(list(Dbins[1:]))
Dbins[-1] = 1e09
Dbins[0] = 0
#desired_samples = ['6037.3']

#desired_samples = parse_timecourse_data.highcoverage_samples
#desired_samples = [parse_timecourse_data.highcoverage_end]
desired_samples = parse_timecourse_data.morteza_samples
#desired_samples = parse_timecourse_data.morteza_samples[1:3]
#desired_samples = [parse_timecourse_data.highcoverage_antibiotic, parse_timecourse_data.highcoverage_end]
#desired_samples = [parse_timecourse_data.highcoverage_start_2, parse_timecourse_data.highcoverage_antibiotic, parse_timecourse_data.highcoverage_postantibiotic, parse_timecourse_data.highcoverage_end]

#pylab.figure(1)
#pylab.xlabel('D')
#pylab.ylabel('G(D)')
        
species_coverage_matrix, samples, species = parse_midas_data.parse_global_marker_gene_coverages()

desired_species = numpy.array(species)

new_species_coverage_matrix = []
for sample_name in desired_samples:
    idx = samples.index(sample_name)
    new_species_coverage_matrix.append(species_coverage_matrix[:,idx])
    
new_species_coverage_matrix = numpy.array(new_species_coverage_matrix)
    

#for i in xrange(0,len(species)):
#        
#        species_coverages = species_coverage_matrix[i,:]
#        if (species_coverages>=min_marker_coverage)

total_barcode_map = barcode_utils.parse_total_barcode_map()
        
outputG = {}
 
G_numeratorss = [] 
G_denominatorss = []

for sample_idx in xrange(0,len(desired_samples)):
    
    # Process each sample separately
    
    sample_name = desired_samples[sample_idx]

    sys.stderr.write("Processing sample %s...\n" % sample_name)

    good_species = parse_midas_data.parse_good_species_list()

    #good_species = desired_species[new_species_coverage_matrix[sample_idx,:]>=config.good_species_min_coverage]
    
    # map from species_name to numeric species id
    # (used to save memory)
    species_id_map = {good_species[i]:i for i in xrange(0,len(good_species))}
    # reverse map from species id to species_name
    id_species_map = good_species
            
    
    sys.stderr.write("Loading depth map...\n")
    barcode_depth_map = barcode_utils.parse_barcode_depth_map(sample_name)    
    sys.stderr.write("Done!\n")
        
    # map from barcode id to total number of SNVs covered in each species
    # organized as a map from species -> n_s
    barcode_species_n_map = {}
    species_Q_map = {}
    
    Btot = total_barcode_map[sample_name][0]
    
    # create barcode->species map
    sys.stderr.write("Collating species barcodes...\n")
    for species_name in good_species:
        
        species_id = species_id_map[species_name]
        
        species_Q_map[species_id] = 0.0
        
        sys.stderr.write("%s...\n" % species_name)
                
        # Make sure barcodes exist for this species at this timepoint.
        if not barcode_utils.barcodes_exist(species_name, sample_name):
            continue
         
        # Load barcodes      
        allele_barcode_map = barcode_utils.parse_allele_barcode_tuples(species_name, sample_name)

        for allele in allele_barcode_map:
        
            if not (allele[-1]=='A' or allele[-1]=='R'):
                # Not a snp, ignore
                continue
                        
            # Add in the barcodes
            for barcode_id, barcode_weight in allele_barcode_map[allele]:
                
                
                if barcode_id not in barcode_species_n_map:
                    barcode_species_n_map[barcode_id] = {}
                
                if species_id not in barcode_species_n_map[barcode_id]:
                    barcode_species_n_map[barcode_id][species_id] = 0
                    
                barcode_species_n_map[barcode_id][species_id] += 1
                
                species_Q_map[species_id] += 1.0/Btot
                
    sys.stderr.write("Done!\n")
    
    Qtot = sum(species_Q_map.values())
    species_complement_Q_map = {species_id: Qtot-species_Q_map[species_id] for species_id in species_Q_map}
    
    G_numerators = numpy.zeros_like(Ds)*1.0
    G_denominators = numpy.zeros_like(Ds)*1.0
    
    for barcode_id in barcode_species_n_map: 
    
        ntot = sum(barcode_species_n_map[barcode_id].values())
        
        numerator_term = 0
        for species_id in barcode_species_n_map[barcode_id]:
            ns = barcode_species_n_map[barcode_id][species_id]
            Qs = species_Q_map[species_id]
            
            numerator_term += ns*(ntot-ns)*1.0/(Qtot-Qs)
            
        depth = barcode_depth_map[barcode_id]
        
        bin_idx = numpy.digitize([depth],bins=Dbins)[0]-1
        
        G_numerators[bin_idx] += numerator_term        
        G_denominators[bin_idx] += ntot
    
    G_numeratorss.append( G_numerators)
    G_denominatorss.append( G_denominators )
        
    
sys.stderr.write("Saving...\n")
print ", ".join([str(d) for d in Dbins])
for sample_idx in xrange(0,len(desired_samples)):
    sample_name = desired_samples[sample_idx]
    G_numerators = G_numeratorss[sample_idx]
    G_denominators = G_denominatorss[sample_idx]
    
    print ", ".join([sample_name]+["%g|%g" % (G_numerators[idx], G_denominators[idx]) for idx in xrange(0,len(G_numerators))])
sys.stderr.write("Done!\n")
#pylab.savefig('test.pdf',bbox_inches='tight')
