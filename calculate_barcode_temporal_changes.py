import sample_utils
import config
import parse_midas_data
import os.path
import pylab
import sys
import numpy
import sfs_utils
import stats_utils
  
from math import log,sqrt
from scipy.special import erfc,bdtr      

import diversity_utils
import gene_diversity_utils
import core_gene_utils
import gzip
import os
import parse_timecourse_data

temporal_change_directory = '%stemporal_changes/' % (config.barcode_directory)
intermediate_filename_template = '%s%s.txt.gz'  

snps_directory=config.barcode_directory+"/barcode_snps/"

min_sample_size = 2

import stats_utils
from math import log10,ceil
from numpy.random import randint


def load_temporal_change_map(species_name):
    
    intermediate_filename = intermediate_filename_template % (temporal_change_directory, species_name)

    temporal_change_map = {}


    if not os.path.isfile(intermediate_filename):
        return temporal_change_map
    
    file = gzip.open(intermediate_filename,"r")
    file.readline() # header
    for line in file:
        items = line.split(",")
        if items[0].strip()!=species_name:
            continue
            
        sample_1 = items[1].strip()
        sample_2 = items[2].strip()
        type = items[3].strip()
        n_tested = float(items[4])
        n_obs = float(items[5])
        n_err = float(items[6])
        sample_pair = (sample_1, sample_2)
        if sample_pair not in temporal_change_map:
            temporal_change_map[sample_pair] = {}
        
        changes = []
        if len(items)<8:
            pass
        else:
            change_strs = items[7:]
            for change_str in change_strs:
            
                subitems = change_str.split(";")
                
                # switch on type of change
                if type=='snps':    
                    gene_name = subitems[0].strip()
                    contig = subitems[1].strip()
                    position = long(subitems[2])
                    variant_type = subitems[3].strip()
                    A1 = float(subitems[4])
                    D1 = float(subitems[5])
                    A2 = float(subitems[6])
                    D2 = float(subitems[7])
                    f_initial = float(subitems[8])
                    f_final = float(subitems[9])
                    changes.append( (gene_name, contig, position, variant_type, A1, D1, A2, D2, f_initial, f_final) )
                    
        temporal_change_map[sample_pair][type] = n_tested, n_obs, n_err, changes
    
    return temporal_change_map

def calculate_mutations_from_temporal_change_map(temporal_change_map,sample_1,sample_2,type='snps'):
    sample_pair = (sample_1,sample_2)
    if sample_pair not in temporal_change_map:
        return 0,0,0,[]
        
    n_tested,n_obs, n_err, snp_changes = temporal_change_map[sample_pair][type]
    return n_tested,n_obs,n_err,snp_changes



if __name__=='__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
    parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
    parser.add_argument("species", help="Name of specific species to run code on")

    args = parser.parse_args()

    debug = args.debug
    chunk_size = args.chunk_size
    species_name=args.species
    good_species_list = [species_name]

    os.system('mkdir -p %s' % temporal_change_directory)
 
    # Load subject and sample metadata
    sys.stderr.write("Loading sample metadata...\n")
    subject_sample_map = sample_utils.parse_subject_sample_map()
    sample_order_map = sample_utils.parse_sample_order_map()
    sys.stderr.write("Done!\n")

    intermediate_filename = intermediate_filename_template % (temporal_change_directory, species_name)
    
    output_file = gzip.open(intermediate_filename,"w")
    # header!
    output_file.write(", ".join(['Species', 'Sample1', 'Sample2', 'Type', 'n_tested','n_obs','n_err', 'Change1', '...']))
    output_file.write("\n")
    
    for species_name in good_species_list:

        #super_good_species_list = parse_midas_data.parse_super_good_species_list()

        #if species_name not in super_good_species_list:
        #    sys.stderr.write("Not in super good species list!\n")
        #    continue
            
        barcode_sample_coverage_map = parse_midas_data.parse_barcode_sample_coverage_map(species_name)
        
        # sys.stderr.write("Loading temporal samples...\n")
        # Only plot samples above a certain depth threshold that are involved in timecourse
        desired_samples = []
        for sample in parse_timecourse_data.all_samples:
            
            if sample not in barcode_sample_coverage_map:
                continue
                
            if barcode_sample_coverage_map[sample] >= config.barcode_min_median_coverage:
                desired_samples.append(sample)
        # sys.stderr.write("Done!\n")
          
        # Analyze SNPs, looping over chunk sizes. 
        # Clunky, but necessary to limit memory usage on cluster

        sys.stderr.write("Loading whitelisted genes...\n")
        personal_core_genes = core_gene_utils.parse_personal_core_genes(species_name)
        sys.stderr.write("Done! %d personal core genes\n" % (len(personal_core_genes)))        
        
        sys.stderr.write("Loading SNVs...\n")
        snp_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_samples=desired_samples, chunk_size=chunk_size,allowed_genes=personal_core_genes, min_present_fraction=config.snv_sweep_min_present_fraction,snps_directory=snps_directory,min_Dbar=config.barcode_min_median_coverage)
        sys.stderr.write("Done! Loaded SNVs from %d genes\n" % len(allele_counts_map))
        
        
        initial_timepoints = []
        initial_final_timepoints = []
        disease_timepoints = []
        disease_final_timepoints = []
        for sample in snp_samples:
            if parse_timecourse_data.sample_epoch_map[sample]=='initial':
                initial_timepoints.append(True)
                initial_final_timepoints.append(False)
                disease_timepoints.append(False)
                disease_final_timepoints.append(False)
                
            elif parse_timecourse_data.sample_epoch_map[sample]=='previous':
                initial_timepoints.append(False)
                initial_final_timepoints.append(False)
                disease_timepoints.append(False)
                disease_final_timepoints.append(False)
            elif parse_timecourse_data.sample_epoch_map[sample]=='disease':
                initial_timepoints.append(False)
                initial_final_timepoints.append(True)
                disease_timepoints.append(True)
                disease_final_timepoints.append(False)
            else:
                initial_timepoints.append(False)
                initial_final_timepoints.append(True)
                disease_timepoints.append(False)
                disease_final_timepoints.append(True)
                
        initial_timepoints = numpy.array(initial_timepoints)
        initial_final_timepoints = numpy.array(initial_final_timepoints)
        disease_timepoints = numpy.array(disease_timepoints)
        disease_final_timepoints = numpy.array(disease_final_timepoints)
        
        # Using disease timepoints as both initial and final   
        #temporal_idxs = numpy.logical_or(initial_timepoints[:,None]*initial_final_timepoints[None,:],disease_timepoints[:,None]*disease_final_timepoints[None,:])
        
        # Without using disease timepoints as initial ones
        temporal_idxs = initial_timepoints[:,None]*initial_final_timepoints[None,:]
        
        snp_difference_map = {}

        snp_difference_matrix = numpy.zeros((len(snp_samples),len(snp_samples)))*1.0
        snp_expected_difference_matrix = numpy.zeros_like(snp_difference_matrix)
        snp_tested_matrix = numpy.zeros_like(snp_difference_matrix)
        
        
        any_differences = 0
        any_expected_differences = 0
        any_tested = 0
        
        Dbars = numpy.array([barcode_sample_coverage_map[sample] for sample in snp_samples])
                
        for gene_name in allele_counts_map.keys():
            if gene_name not in personal_core_genes:
                sys.stderr.write("Shouldn't get here: %s\n" % gene_name)
                continue
                
            for variant_type in allele_counts_map[gene_name].keys():
            
                allele_counts_matrix = allele_counts_map[gene_name][variant_type]['alleles']
                
                if len(allele_counts_matrix)==0:
                    continue
                    
                depths = allele_counts_matrix.sum(axis=2)
                safe_depths = depths+(depths==0)
                alts = allele_counts_matrix[:,:,0]
                freqs = alts*1.0/(safe_depths)
                logcopynums = numpy.log(safe_depths/Dbars) # I feel like it's ok to use safe_depths here, because Dbars is at least 10, and depths/Dbars should be at most 3 for depths>0
                
                passed_depths = depths*(depths>=config.barcode_min_fixation_coverage)
                passed_alts = alts*(depths>=config.barcode_min_fixation_coverage)
                passed_refs = (depths-alts)*(depths>=config.barcode_min_fixation_coverage)
                
                weak_polymorphic_sites = numpy.logical_not(numpy.logical_or((passed_alts<=0.2*passed_depths).all(axis=1),(passed_refs<=0.2*passed_depths).all(axis=1))) 
                
                # Things with at least 5 reads at both sites, and not too bad of a copynum difference
                passed_sites = passed_depths[:,:,None]*passed_depths[:,None,:]*weak_polymorphic_sites[:,None,None]
                
                # Make sure copynum doesn't change too much
                passed_sites = numpy.logical_and(passed_sites,  (numpy.fabs(logcopynums[:,:,None]-logcopynums[:,None,:])<=log(2)))
                
                # Calculate apparent avg fbar between pairs of sites
                fbars = (alts[:,:,None]+alts[:,None,:])*1.0/(safe_depths[:,:,None]+safe_depths[:,None,:])
                safe_fbars = numpy.clip(fbars,1e-02,1-1e-02)
                
                total_alts = alts.sum(axis=1)
                total_depths = depths.sum(axis=1)
                fbarbar = (total_alts)*1.0/(total_depths+(total_depths==0))
                
                harmonic_depths = 2.0/(1.0/safe_depths[:,:,None]+1.0/safe_depths[:,None,:])
                
                f0s = freqs[:,:,None]*numpy.ones_like(fbars)
                ffs = freqs[:,None,:]*numpy.ones_like(fbars)
                
                d0s = numpy.array(depths[:,:,None]*numpy.ones_like(fbars),dtype=numpy.int32)
                dfs = numpy.array(depths[:,None,:]*numpy.ones_like(fbars),dtype=numpy.int32)
                
                # The number of alts you'd need at the initial timepoint
                # to have initial frequency <=0.2
                initial_lower_ks = numpy.array(numpy.floor(0.2*depths)[:,:,None]*numpy.ones_like(fbars),dtype=numpy.int32)
                # ... initial frequency >=0.8
                initial_upper_ks = numpy.array(numpy.ceil(0.8*depths)[:,:,None]*numpy.ones_like(fbars),dtype=numpy.int32)
                 
                # The number of alts you'd need at the final timepoint
                # to have final frequency <=0.3
                final_lower_ks = numpy.array(numpy.floor(0.3*depths)[:,None,:]*numpy.ones_like(fbars),dtype=numpy.int32)
                # ... >=0.7
                final_upper_ks = numpy.array(numpy.ceil(0.7*depths)[:,None,:]*numpy.ones_like(fbars),dtype=numpy.int32)
                
                # Probability of having freq <=0.2 at initial timepoint
                # conditioned on true frequency being fbar
                p_initial_lower = bdtr(initial_lower_ks, d0s, fbars)
                # ... >= 0.8
                p_initial_upper = bdtr(d0s-initial_upper_ks, d0s, 1-fbars)
                
                # ... final freq <= 0.3
                p_final_lower = bdtr(final_lower_ks, dfs, fbars)
                # ... final freq >= 0.7
                p_final_upper = bdtr(dfs-final_upper_ks, dfs, 1-fbars)
                
                # Probability of an apparent change under the null hypothesis
                null_probabilities = (p_initial_lower*p_final_upper)+(p_initial_upper*p_final_lower)
                
                # Based on harmonic depths (i.e., Gaussian approx)
                #passed_sites = numpy.logical_and(passed_sites, harmonic_depths>=10)
                
                # Based on binomial p-value calculation
                passed_sites = numpy.logical_and(passed_sites, null_probabilities<config.barcode_max_null_fixation_pvalue)
                
                #sigmas = numpy.sqrt(2*safe_fbars*(1-safe_fbars)/harmonic_depths)
                  
                #dfs = numpy.fabs(freqs[:,:,None]-freqs[:,None,:])
                #diffs = numpy.logical_or((freqs[:,:,None]>=0.7)*(freqs[:,None,:]<=0.3),(freqs[:,:,None]<=0.3)*(freqs[:,None,:]>=0.7))*(dfs>=0.6)*passed_sites
                
                alt_diffs = (f0s<=0.2)*(ffs>=0.7)*passed_sites
                ref_diffs = (f0s>=0.8)*(ffs<=0.3)*passed_sites  
                diffs = numpy.logical_or(alt_diffs, ref_diffs)      
                
                num_diffs = diffs.sum(axis=0)
                
                # of sites where you expect to see a change,
                # even under null hypothesis
                num_expected_diffs = (null_probabilities*passed_sites).sum(axis=0)
                
                snp_difference_matrix += num_diffs
                snp_expected_difference_matrix += num_expected_diffs
                snp_tested_matrix += passed_sites.sum(axis=0)
                
                any_differences += ((diffs.sum(axis=2).sum(axis=1))>0).sum(axis=0)
                any_expected_differences += num_expected_diffs.sum()
                any_tested += ((passed_sites.sum(axis=2).sum(axis=1))>0).sum(axis=0)
                
        
                for idx in xrange(0,diffs.shape[0]):
                    # Looping over sites...
                    diff_matrix = diffs[idx,:,:]
                    
                    if (diff_matrix[temporal_idxs].any()):
                        for i in xrange(0,len(snp_samples)):
                            for j in xrange(i+1,len(snp_samples)):
                                if temporal_idxs[i,j] and diffs[idx,i,j]>0:
                                
                                    
                                
                                    if (i,j) not in snp_difference_map:
                                        snp_difference_map[(i,j)] = []
                                    
                                    # record final and initial frequency
                                    final_freq = freqs[idx,passed_depths[idx,:]>0][-1]
                                    initial_freq = freqs[idx,passed_depths[idx,:]>0][0]
                                
                                    if ref_diffs[idx,i,j]>0:
                                        # need to re-polarize
                                        final_freq = 1-final_freq
                                        initial_freq = 1-initial_freq
                                
                                    snp_change = (gene_name, allele_counts_map[gene_name][variant_type]['locations'][idx], variant_type, (alts[idx,i],depths[idx,i]), (alts[idx,j],depths[idx,j]), (initial_freq, final_freq) )
                        
                                    snp_difference_map[(i,j)].append( snp_change)
                                    
        record_str_items = [species_name, 'any', 'any', 'snps', "%g" % any_tested, "%g" % any_differences, "%g" % any_expected_differences]
        record_str = ", ".join(record_str_items)
        output_file.write(record_str)
        output_file.write("\n")
            
        for i in xrange(0,len(snp_samples)):
            sample_i = snp_samples[i]
            for j in xrange(i+1,len(snp_samples)):
                sample_j = snp_samples[j]
                
                if not temporal_idxs[i,j]:
                    continue
                
                #print (i,j), sample_i, sample_j
                    
                snp_strs = []
                if (i,j) in snp_difference_map:
                    
                    #print snp_difference_matrix[i,j], len(snp_difference_map[(i,j)])
                
                    for snp_change in snp_difference_map[(i,j)]:
        
                        gene_name, location, variant_type, allele_counts_1, allele_counts_2, bookend_freqs = snp_change
                        contig = location[0]
                        position = location[1]
            
                        A1,D1 = allele_counts_1
                        A2,D2 = allele_counts_2
            
                        f_initial, f_final = bookend_freqs
            
                        snp_str = ('%s;%s;%d;%s;%d;%d;%d;%d;%g;%g' % (gene_name, contig, position, variant_type, A1, D1, A2, D2, f_initial, f_final))
            
                        snp_strs.append(snp_str)
            
                record_str_items = [species_name, sample_i, sample_j, 'snps', "%g" % snp_tested_matrix[i,j], "%g" % snp_difference_matrix[i,j], "%g" % snp_expected_difference_matrix[i,j]] + snp_strs
                record_str = ", ".join(record_str_items)
                output_file.write(record_str)
                output_file.write("\n")
            
            
    sys.stderr.write("Done looping over species!\n")
    output_file.close()
    sys.stderr.write("Done!\n")
    
    # testing loading of intermediate file
    temporal_change_map = load_temporal_change_map(good_species_list[0])
 
