import sample_utils 
import config
import parse_midas_data
import parse_timecourse_data
import os.path
import pylab
import sys
import numpy

import diversity_utils
import gene_diversity_utils

import stats_utils
from math import log10,ceil
from numpy.random import randint

import core_gene_utils
import gzip
import os

import calculate_snp_prevalences

min_coverage = config.min_median_coverage
min_sample_size = 10

#species_name = 'Bacteroides_vulgatus_57955'
#species_name = 'Bacteroides_uniformis_57318'
#target_samples = parse_timecourse_data.highcoverage_samples
species_name = 'Parabacteroides_merdae_56972'
target_samples = ['4023', '6041']

#species_name = 
debug = False

snp_samples = parse_timecourse_data.highcoverage_samples

sys.stderr.write("Loading core genes...\n")
core_genes = core_gene_utils.parse_core_genes(species_name)
non_shared_genes = core_gene_utils.parse_non_shared_reference_genes(species_name)
shared_pangenome_genes = core_gene_utils.parse_shared_genes(species_name)
sys.stderr.write("Done! Core genome consists of %d genes\n" % len(core_genes))
sys.stderr.write("%d shared genes and %d non-shared genes\n" % (len(shared_pangenome_genes), len(non_shared_genes)))

snv_freq_map = calculate_snp_prevalences.parse_population_freqs(species_name,polarize_by_consensus=True)
    
# Now calculate gene differences
# Load gene coverage information for species_name
sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name,allowed_samples=target_samples, disallowed_genes=shared_pangenome_genes)
sys.stderr.write("Done! Loaded %d genes\n" % len(gene_names))

gene_copynums = gene_depth_matrix*1.0/marker_coverages[None,:]

gene_copynum_map = {}
for gene_idx in xrange(0,len(gene_names)):
    gene_copynum_map[gene_names[gene_idx]] = (gene_copynums[gene_idx].min(),numpy.median(gene_copynums), gene_copynums[gene_idx].max())

snp_samples = target_samples

sys.stderr.write("Loading %d samples...\n" % len(snp_samples))
     
snp_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_samples=snp_samples,allowed_genes=core_genes)
        
target_genes = []
target_gene_numbers = []
for gene_name in allele_counts_map:
    items = gene_name.split(".")
    number = long(items[-1])
    if gene_copynum_map[gene_name][0]>0.7:
        target_genes.append(gene_name)
        target_gene_numbers.append(number)
    
target_gene_numbers, target_genes = zip(*sorted(zip(target_gene_numbers, target_genes)))

within_locations = []
within_genes = []
within_freqs = []
within_prevalences = []
within_targets = []

gene_sites_map = {}

for gene_name in target_genes:
    for var_type in passed_sites_map[gene_name]:
        if gene_name not in gene_sites_map:
            gene_sites_map[gene_name] = 0
            
        gene_sites_map[gene_name] += numpy.median(passed_sites_map[gene_name][var_type]['sites'])
        
for gene_name in target_genes:
    gene_number = long(gene_name.split(".")[-1])
    for var_type in allele_counts_map[gene_name]:
        allele_counts_matrix = allele_counts_map[gene_name][var_type]['alleles']
        if len(allele_counts_matrix) == 0:
            continue
            
        
        locations = allele_counts_map[gene_name][var_type]['locations']
        alts = allele_counts_matrix[:,:,0]
        depths = allele_counts_matrix[:,:,:].sum(axis=2)
        
        freqs = alts*1.0/(depths+(depths==0))
        
        for snp_idx in xrange(0,len(locations)):
            location_tuple = locations[snp_idx]
            if location_tuple in snv_freq_map:
                prevalence = snv_freq_map[location_tuple]
            else:
                prevalence = 0
            location = locations[snp_idx][1]
            
            good_idxs = (depths[snp_idx,:]>0)
            
            if len(good_idxs)==0:
                continue
            
            freq_vector = freqs[snp_idx,good_idxs]
            polarized_freq_vector = numpy.fmin(freq_vector,1-freq_vector)
            freq = freq_vector[polarized_freq_vector.argmax()]
            
            if numpy.median(depths[snp_idx,:])>0: 
        
                
                within_locations.append((gene_number,location))
                within_freqs.append(freq)
                within_genes.append(gene_name)
                within_targets.append(gene_sites_map[gene_name])
                within_prevalences.append(prevalence)

within_locations, within_genes, within_freqs, within_prevalences, within_targets = zip(*sorted(zip(within_locations, within_genes, within_freqs, within_prevalences, within_targets)))

for snp_idx in xrange(0,len(within_locations)):
    print within_locations[snp_idx][1], within_genes[snp_idx], within_freqs[snp_idx], within_prevalences[snp_idx], within_targets[snp_idx]