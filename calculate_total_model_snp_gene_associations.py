###########################################################################
#
# The purpose of this script is to link within-species SNV changes to core  genes of different species.
# prints data to stdout (small enough to fit in dropbox)
# conventially fixation_barcode_output.txt
#
###########################################################################


import sys
import pylab
import numpy
import parse_midas_data
import parse_timecourse_data
import stats_utils
import barcode_utils
import gzip
import collections

import os.path
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

pangenome_species = parse_midas_data.parse_pangenome_species()
#pangenome_species = ['Phascolarctobacterium_sp_59817']
    
###########################################################################
#
# Load fixations to track
#
###########################################################################

        

corrected = True
core_genes_only = True
external_core_genes=True
min_S = 1 # or 1
Pstar = 0.1


# load candidate SNPs from supplied file

# long snp map
# long snp = (species, "contig|location") tuple (allele is stored separately)
# id = int
longsnp_id_map = {}
id_longsnp_map = []

longsnp_self_longgenes = []

# long gene map
# long gene = (species, gene_name) tuple 
# id = int
longgene_id_map = {}
id_longgene_map = []

# core longgenes
# Load core genes across those species
sys.stderr.write("Loading core genes...\n")
import core_gene_utils
core_longgenes = set()
for species_name in pangenome_species:
    core_genes = core_gene_utils.parse_core_genes(species_name, external_filtering=external_core_genes)
    for gene_name in core_genes:
        core_longgenes.add((species_name, gene_name))
sys.stderr.write("Done!\n")

sys.stderr.write("Loading fixations...\n")
import calculate_barcode_within_species_fixations as calculate_within_species_fixations
for species_name in pangenome_species:
    
    snp_changes, gene_changes = calculate_within_species_fixations.load_within_species_fixations(species_name,allowed_epochs=parse_timecourse_data.initial_plus_previous_epoch_intervals)
    snp_changes = list(snp_changes)
    shuffle(snp_changes)    
    if len(snp_changes)>1000:
        snp_changes = snp_changes[0:1000] 
        
    for snp_change in sorted(snp_changes):
        contig = snp_change[0]
        location = snp_change[1]
        gene_name = snp_change[2]
        
        longsnp = (species_name, contig, location)
    
        if longsnp not in longsnp_id_map:
            id = len(id_longsnp_map)
            longsnp_id_map[longsnp]=id
            id_longsnp_map.append(longsnp)
            longsnp_self_longgenes.append((species_name, gene_name))
            
sys.stderr.write("Done!\n")


#desired_samples = parse_timecourse_data.highcoverage_samples
#desired_samples = [parse_timecourse_data.highcoverage_end]
desired_samples = parse_timecourse_data.morteza_samples
#desired_samples = parse_timecourse_data.morteza_samples[1:3]
#desired_samples = [parse_timecourse_data.highcoverage_antibiotic, parse_timecourse_data.highcoverage_end]
#desired_samples = [parse_timecourse_data.highcoverage_start_2, parse_timecourse_data.highcoverage_antibiotic, parse_timecourse_data.highcoverage_postantibiotic, parse_timecourse_data.highcoverage_end]

from barcode_utils import calculate_M_from_D


gamete_idx_map = {('R','R'): 0, ('R','A'):1, ('A','R'):2, ('A','A'):3  }
allele_idx_map = {'R':0, 'A':1}
# Memory requirements are big, so we have to break things up into smaller batches
# of unassigned genes. The output is then grouped into a single file.
batch_size = 100000
num_batches = long(ceil(len(id_longsnp_map)*1.0/batch_size))

if num_batches>1.5 and (len(id_longsnp_map)%batch_size)*1.0/batch_size < 0.1:
    # just expand the batch size a little bit
    batch_size = long(ceil(len(id_longsnp_map)*1.0/(num_batches-1)))
    num_batches = long(ceil(len(id_longsnp_map)*1.0/batch_size))

sys.stderr.write("Divided %d snps into %d batches of size %d.\n" % (len(id_longsnp_map), num_batches, batch_size))

for batch in xrange(0,num_batches):
    
    sys.stderr.write("Processing batch %d...\n" % batch)
    
    snp_shared_barcodes = {}
    snp_total_barcodes = {}
    snp_total_M = {}
    for i in xrange(0,len(id_longsnp_map)):
        if (i/batch_size) == batch:
            snp_shared_barcodes[i] = {'A': collections.Counter(), 'R': collections.Counter()}
            snp_total_barcodes[i] = {'A': 0, 'R': 0}
            snp_total_M[i] = {'A': 0, 'R': 0}
    
    sys.stderr.write("(%d SNVs)\n" % len(snp_total_barcodes))
    
    total_total_barcodes = 0
    target_longgene_barcode_fraction_map = {}
    
    for sample_idx in xrange(0,len(desired_samples)):
        # Process each sample separately
        # (combine results at end)
    
        sample_name = desired_samples[sample_idx]

        num_significant_hits = 0 
        
        sys.stderr.write("Processing sample %s...\n" % sample_name)
        sys.stderr.write("Loading depth map...\n")
        barcode_depth_map = barcode_utils.parse_barcode_depth_map(sample_name,corrected=corrected,min_depth=2.5)
        sys.stderr.write("Done!\n")
        
        if len(barcode_depth_map)==0:
            continue
        
        sys.stderr.write("Postprocessing depth map...\n")
        
        total_barcodes = len(barcode_depth_map)
        total_total_barcodes += total_barcodes
        avg_M = 0.0
        for D in barcode_depth_map.itervalues():
            avg_M += calculate_M_from_D(D)
        avg_M /= total_barcodes
        sys.stderr.write("Done!\n")
         
        # create barcode->species map
        sys.stderr.write("Collating species barcodes...\n")
        
        barcode_longgene_ids_map = {} # map from barcode id to longgenes
        
        focal_longsnp_barcode_map = {} # Map from longsnp to barcode ids that that feature has
        
        for species_name in pangenome_species:
        
            # 'new_species' is set of de-novo assembled genes
            # do not use these for now
            if species_name=='new_species':
                continue
            
            #if species_name not in desired_speciess:
            #    continue
        
            # Make sure barcodes exist for this species at this timepoint.
            if not barcode_utils.barcodes_exist(species_name, sample_name):
                continue
         
            # Load barcodes      
            allele_barcode_map = barcode_utils.parse_allele_barcode_tuples(species_name, sample_name,corrected=corrected)

            for allele in allele_barcode_map:
        
                if (allele[-1]=='A' or allele[-1]=='R'):
                    # A SNV!
                    
                    snp_str = allele[:-2]
                    polarization = allele[-1]
                    items = snp_str.split("|")
                    contig = items[0]
                    location = long(items[1])
                    longsnp = (species_name, contig, location)
                
                    #if location==68991:
                    #    print longsnp, polarization, (longsnp in longsnp_id_map)
                
                    #print longsnp, allele, longsnp_id_map.keys()
                
                    # Not one of the SNVs we are tracking
                    if longsnp not in longsnp_id_map: 
                        #print longsnp
                        continue
                    
                    longsnp_id = longsnp_id_map[longsnp]
                    
                    # Not one of the ones in this batch
                    if longsnp_id not in snp_shared_barcodes: 
                        continue
                        
                    if len(allele_barcode_map[allele])==0:
                        continue
                        
                    if longsnp_id not in focal_longsnp_barcode_map:
                        focal_longsnp_barcode_map[longsnp_id] = {}
                            
                    if polarization not in focal_longsnp_barcode_map[longsnp_id]:
                        focal_longsnp_barcode_map[longsnp_id][polarization] = set()
                            
                    
                    # Add in the barcodes to track.
                    for barcode_id, barcode_weight in allele_barcode_map[allele]:
                        #print barcode_id
                        # don't include barcode if too little coverage
                        # Poor man's error correction
                        if barcode_id not in barcode_depth_map: 
                            continue
                           
                        # Add barcode to list
                        focal_longsnp_barcode_map[longsnp_id][polarization].add(barcode_id)
                    
                
                else:
                    # A gene!
                
                    gene_name = allele
                    longgene = (species_name, gene_name)
                
                    # If specified, only look at genes in the core genome
                    if (core_genes_only) and (longgene not in core_longgenes):
                        continue
                
                    if longgene not in longgene_id_map:
                        # add a record
                        longgene_id = len(id_longgene_map)
                        id_longgene_map.append(longgene)
                        longgene_id_map[longgene] = longgene_id
                    
                    longgene_id = longgene_id_map[longgene]
                    
                
                    # Only look at snps that have barcodes
                    if len(allele_barcode_map[allele])==0:
                        continue
            
                    # estimated # of fragments per barcode
                    total_B = 0

                    for barcode_id, barcode_weight in allele_barcode_map[allele]:
                
                        # don't include barcode if too little coverage
                        # Poor man's error correction
                        if barcode_id not in barcode_depth_map: 
                            continue
                    
                        total_B += 1
                        
                        if barcode_id not in barcode_longgene_ids_map:
                            barcode_longgene_ids_map[barcode_id] = set()
                            #barcode_depth_map[barcode_id] = 0
                    
                        barcode_longgene_ids_map[barcode_id].add(longgene_id)
                    
                    if longgene_id not in target_longgene_barcode_fraction_map:
                        target_longgene_barcode_fraction_map[longgene_id]=0
                        
                    target_longgene_barcode_fraction_map[longgene_id] += total_B
        
        
    
        # Make sure we loaded some barcodes...
        if len(barcode_longgene_ids_map)==0:
            continue
                
        total_longsnps = len(target_longgene_barcode_fraction_map)
    
        log_corrected_Pstar = log(Pstar/total_longsnps/total_longsnps)
                
        sys.stderr.write("Done! Loaded %d barcodes\n" % len(barcode_longgene_ids_map))
        sys.stderr.write("(other count: %d)\n" % total_barcodes)
        sys.stderr.write("Checking linkage with %d target snps\n" % total_longsnps)
        
        sys.stderr.write("Looping through target features...\n")
        
        for focal_longsnp_id in focal_longsnp_barcode_map:
        
            for focal_polarization in focal_longsnp_barcode_map[focal_longsnp_id]:
        
                # B = total # barcodes for this feature
                B = len(focal_longsnp_barcode_map[focal_longsnp_id][focal_polarization])
            
                # Get total M
                Mtot = 0.0
                for barcode_id in focal_longsnp_barcode_map[focal_longsnp_id][focal_polarization]:
                    Mtot += calculate_M_from_D(barcode_depth_map[barcode_id])
            
                snp_total_barcodes[focal_longsnp_id][focal_polarization] += B
                snp_total_M[focal_longsnp_id][focal_polarization] += Mtot
            
                # Get barcode ids that map to this gene    
                barcode_ids = []
                for barcode_id in focal_longsnp_barcode_map[focal_longsnp_id][focal_polarization]:
                
                    # If it doesn't map to anywhere we know of, 
                    # skip it (saves some time later)
                    if barcode_id not in barcode_longgene_ids_map:
                        continue
                     
                    barcode_ids.append(barcode_id)
            
                num_barcodes = len(barcode_ids)
            
                # only create entries if there are at least two barcodes
                if num_barcodes<1.5:
                    continue
                
                #sys.stderr.write("Processing %d barcodes for %s...\n" % (num_barcodes,allele))
               
                # Counts number of barcodes that map to different genes
                # within species
                for barcode_id in barcode_ids:
                    snp_shared_barcodes[focal_longsnp_id][focal_polarization].update( barcode_longgene_ids_map[barcode_id] )
                    
            # Done looping through focal polarizations
        # Done looping through focal snps              
    # Done looping through samples!          
    
    # Now output stuff!
    for focal_longsnp_id in sorted(snp_total_barcodes):
        species, contig, location = id_longsnp_map[focal_longsnp_id]
        self_longgene = longsnp_self_longgenes[focal_longsnp_id]
        gene_name = self_longgene[1]
        print "%s|%s|%s|%d" % (species, contig, gene_name, location)
        for polarization in snp_total_barcodes[focal_longsnp_id]:
            
            output_strs = [polarization, str(snp_total_barcodes[focal_longsnp_id][polarization]), str(snp_total_M[focal_longsnp_id][polarization])]
            
            # Tuple of key,counts sorted in descending order of counts    
            for longgene_id, longgene_count in snp_shared_barcodes[focal_longsnp_id][polarization].most_common(10):
                longgene = id_longgene_map[longgene_id]
                
                # Don't look at linkage to its own gene!
                if longgene==self_longgene:
                    continue
                
                if longgene_count < 5:
                    continue
                    
                if longgene_count < 0.05*snp_total_barcodes[focal_longsnp_id][polarization]:
                    continue
                            
                q = target_longgene_barcode_fraction_map[longgene_id]*1.0/total_total_barcodes
                
                output_strs.append("%s|%s|%d|%.2f" % (longgene[0], longgene[1], longgene_count, q))
            
            print ", ".join(output_strs)    
            
                    
            
sys.stderr.write("Done!\n")

