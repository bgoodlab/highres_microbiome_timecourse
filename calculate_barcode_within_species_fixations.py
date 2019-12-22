import numpy
import config
import diversity_utils
import parse_timecourse_data
import parse_midas_data
import sys
import gzip
import core_gene_utils
from scipy.special import bdtr

import calculate_barcode_temporal_changes 
default_filename = config.data_directory+"within_species_fixations_barcodes.txt.gz"

def load_within_species_fixations(desired_species, allowed_epochs=[]):
    
    allowed_epochs = set(allowed_epochs)
    snp_changes = []
    epochs = []
    
    file = gzip.open(default_filename,"r")
    for line in file:
        items = line.split(",")
        if len(items)<2:
            continue
            
        species = items[0].strip()
        type = items[1].strip()
        epoch_strs = items[2].strip()
        epochs = [subitem.strip() for subitem in epoch_strs.split(";")]
        change_strs = items[3:]
        
        if species==desired_species:
                        
            if type=='snps':
                
                for snp_change_str in change_strs:
                    
                    subitems = snp_change_str.split(";")
                    contig = subitems[0].strip()
                    location = long(subitems[1])
                    gene_name = subitems[2].strip()
                    variant_type = subitems[3].strip()
                    #initial_h = float(subitems[4])
                    initial_freq = float(subitems[4])
                    final_freq = float(subitems[5])
                    epochs = set()
                    for epoch_item in subitems[6:]:
                        epochs.add(epoch_item.strip())
                    
                    if (len(allowed_epochs)==0) or (len(epochs & allowed_epochs) > 0):
                        snp_changes.append((contig, location, gene_name, variant_type, initial_freq, final_freq))
                    
            if type=='genes':
                 

                 for gene_change_str in change_strs:
                    
                    subitems = gene_change_str.strip().split(";")
                    gene_name = subitems[0]
                    epochs = set(subitems[1:])
                    
                    if len(epochs & allowed_epochs) > 0:
                        gene_changes.append(gene_name)
    
    snp_changes = set(snp_changes)
    
    return snp_changes, epochs
    
if __name__=='__main__':
    
    sample_time_map = parse_timecourse_data.parse_sample_time_map()
    ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, parse_timecourse_data.all_samples, min_time=-10000)
    all_samples = numpy.array(parse_timecourse_data.all_samples)[sample_idxs]

    desired_samples = parse_timecourse_data.all_samples
    
    file = gzip.open(default_filename,"w")
        
    good_species_list = parse_midas_data.parse_good_species_list()
    for species_name in good_species_list:
    
        sys.stderr.write("Processing %s...\n" % species_name)
    
        #sys.stderr.write("Loading pre-computed temporal changes for %s...\n" % species_name)
        temporal_change_map = calculate_barcode_temporal_changes.load_temporal_change_map(species_name)
        #sys.stderr.write("Done!\n")
    
    
        snp_change_map = {}
        passed_epochs = set()
        
        for i in xrange(0,len(desired_samples)):
    
            sample_i = desired_samples[i]
            epoch_i = parse_timecourse_data.sample_epoch_map[sample_i]
        
            for j in xrange(i+1,len(desired_samples)):
            
                sample_j = desired_samples[j]
                epoch_j = parse_timecourse_data.sample_epoch_map[sample_j]
    
                current_epoch= ('%s_%s' % (epoch_i, epoch_j))
    
                num_tested, num_snp_changes, nerr, snp_changes = calculate_barcode_temporal_changes.calculate_mutations_from_temporal_change_map(temporal_change_map, sample_i, sample_j)
                
                if num_tested<0.5:
                    # No data, skip!
                    continue
                
                failed_fpr = (nerr > max([0.1, 0.1*num_snp_changes]))
                passed_fpr = (not failed_fpr)
                
                if passed_fpr:
                    passed_epochs.add(current_epoch)
                
                if (num_snp_changes<0.5):
                    # Only want to make a list of things that changed.
                    continue
                
                passed_fpr = passed_fpr and (bdtr(num_tested-num_snp_changes, num_tested, 1.0-nerr*1.0/num_tested) < 5e-2/60)
                
                if not passed_fpr:
                    continue
                
                for gene_name, contig, position, variant_type, A1, D1, A2, D2, f0,ff in snp_changes:
                
                    record = (contig, position, gene_name, variant_type)
            
                    #Dmin = min([D1,D2])
                    
                    f1 = A1*1.0/D1
                    f2 = A2*1.0/D2
                    h1 = 2*f1*(1-f1) 
                
                    if f1<f2:
                        pair = (sample_i, sample_j)
                    else:
                        pair = (sample_j, sample_i)
                
                    
                    if record not in snp_change_map:
                        snp_change_map[record] = {'epochs': set(), 'initial_hs': [], 'pairs': [set(),set()], 'initial_freqs': [], 'final_freqs':[]}
                    
                    snp_change_map[record]['epochs'].add(current_epoch)
                    snp_change_map[record]['initial_hs'].append(h1)
                    snp_change_map[record]['pairs'][0].add( pair[0])
                    snp_change_map[record]['pairs'][1].add( pair[1])
                    snp_change_map[record]['initial_freqs'].append(f0)
                    snp_change_map[record]['final_freqs'].append(ff)
                    
                    
                                        
        # Done looping over samples
        # print record for species
        
        all_records = set(snp_change_map.keys())
        
        output_str = "%s, snps" % species_name
        
        output_str += (", "+ "; ".join([e for e in sorted(passed_epochs)]))
        
        for record in sorted(snp_change_map):
        
            contig, position, gene_name, variant_type = record
           
            epochs = list(sorted(snp_change_map[record]['epochs']))
        
            snp_change_map[record]['initial_hs'] = numpy.array(snp_change_map[record]['initial_hs'])
            median_h = numpy.median(snp_change_map[record]['initial_hs'])
        
            # Transform to freq
            initial_freq = (1-(1-2*median_h)**0.5)/2
            
            #initial_freq = snp_change_map[record]['initial_freqs'][0]
            final_freq = snp_change_map[record]['final_freqs'][0]
        
            record_str = ";".join([contig, str(position), gene_name, variant_type, str(initial_freq), str(final_freq)]+epochs)
        
            output_str += (", "+record_str)
        
        file.write(output_str)
        file.write("\n")    
    
    file.close()         
    