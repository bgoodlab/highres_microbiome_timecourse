import numpy
import config
import diversity_utils
import parse_timecourse_data
import parse_midas_data
import sys
import gzip
import core_gene_utils

import calculate_temporal_changes as calculate_temporal_changes
default_filename = config.data_directory+"within_species_fixations.txt.gz"

def load_within_species_fixations(desired_species, allowed_epochs=[]):
    
    allowed_epochs = set(allowed_epochs)
    snp_changes = []
    gene_changes = [] 
    
    file = gzip.open(default_filename,"r")
    for line in file:
        items = line.split(",")
        if len(items)<2:
            continue
            
        species = items[0].strip()
        type = items[1].strip()
        change_strs = items[2:]
        
        if species==desired_species:
            
            
            if type=='snps':
                
                for snp_change_str in change_strs:
                    
                    subitems = snp_change_str.split(";")
                    contig = subitems[0].strip()
                    location = long(subitems[1])
                    gene_name = subitems[2].strip()
                    variant_type = subitems[3].strip()
                    initial_h = float(subitems[4])
                    epochs = set()
                    for epoch_item in subitems[5:]:
                        epochs.add(epoch_item.strip())
                    
                    if (len(allowed_epochs)==0) or (len(epochs & allowed_epochs) > 0):
                        snp_changes.append((contig, location, gene_name, variant_type, initial_h))
                    
            if type=='genes':
                 

                 for gene_change_str in change_strs:
                    
                    subitems = gene_change_str.strip().split(";")
                    gene_name = subitems[0]
                    epochs = set(subitems[1:])
                    
                    if len(epochs & allowed_epochs) > 0:
                        gene_changes.append(gene_name)
    
    snp_changes = set(snp_changes)
    gene_changes = set(gene_changes)
    
    return snp_changes, gene_changes
    
if __name__=='__main__':
    
    sample_time_map = parse_timecourse_data.parse_sample_time_map()
    ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, parse_timecourse_data.all_samples, min_time=-10000)
    all_samples = numpy.array(parse_timecourse_data.all_samples)[sample_idxs]

    desired_samples = parse_timecourse_data.morteza_samples+parse_timecourse_data.kuleshov_samples
    
    file = gzip.open(default_filename,"w")
        
    good_species_list = parse_midas_data.parse_good_species_list()
    #good_species_list = ['Phascolarctobacterium_sp_59817']
    for species_name in good_species_list:
    
        sys.stderr.write("Processing %s...\n" % species_name)
        
        sys.stderr.write("Loading core genes...\t")
        personal_core_genes = core_gene_utils.parse_personal_core_genes(species_name)
        shared_genes = core_gene_utils.parse_shared_genes(species_name)
        sys.stderr.write("Done!\n")
        
    
        #sys.stderr.write("Loading pre-computed temporal changes for %s...\n" % species_name)
        temporal_change_map = calculate_temporal_changes.load_temporal_change_map(species_name)
        #sys.stderr.write("Done!\n")
    
        snp_change_map = {}
        gene_change_map = {}
    
        for i in xrange(0,len(desired_samples)):
    
            sample_i = desired_samples[i]
        
            if sample_i not in parse_timecourse_data.epoch_samples['initial']:
                continue # only look at changes between initial and final
    
            for j in xrange(i+1,len(desired_samples)):
            
                sample_j = desired_samples[j]
    
                if sample_j in parse_timecourse_data.kuleshov_samples:
                    target_sample_i = sample_j
                    target_sample_j = sample_i
                else:
                    target_sample_i = sample_i
                    target_sample_j = sample_j
                    
                L, perr, mutations, reversions = calculate_temporal_changes.calculate_mutations_reversions_from_temporal_change_map(temporal_change_map, target_sample_i, target_sample_j)
        
                nerr = L*perr
        
                snp_changes = mutations+reversions
                num_snp_changes = len(snp_changes)
        
                gene_L, gene_perr, gains, losses = calculate_temporal_changes.calculate_gains_losses_from_temporal_change_map(temporal_change_map, target_sample_i, target_sample_j)
        
                gene_nerr = gene_L*gene_perr
        
                gene_changes = gains+losses
                num_gene_changes = len(gene_changes)
        
                #print sample_i, sample_j, perr, gene_perr
        
                if (perr<-0.5) or (gene_perr < -0.5):
                    continue
        
                if (nerr > max([0.5, 0.1*num_snp_changes])) or (gene_nerr > max([0.5, 0.1*num_gene_changes])):
                    continue # Only take things with low-ish FPR
            
                #print "Made it this far!"
            
        
                if (num_snp_changes<0.5) and (num_gene_changes<0.5):
                    # Only want to make a list of things that changed.
                    continue
            
                #print "Made it!"
                current_epoch_id = -1
                # Assign an epoch based on final timepoint
                for epoch in parse_timecourse_data.epoch_samples:
                    if sample_j in parse_timecourse_data.epoch_samples[epoch]:
                        current_epoch = parse_timecourse_data.epoch_order_map[epoch]
                        break
            
                for gene_name, contig, position, variant_type, A1, D1, A2, D2 in snp_changes:
                
                    record = (contig, position, gene_name, variant_type)
                
                    f1 = A1*1.0/D1
                    h1 = 2*f1*(1-f1) 
                
                    if gene_name not in personal_core_genes:
                        continue
                 
                    if record not in snp_change_map:
                        snp_change_map[record] = {'epochs': set(), 'initial_hs': []}
                    
                    snp_change_map[record]['epochs'].add(current_epoch)
                    snp_change_map[record]['initial_hs'].append(h1)
            
                for gene_name, D1, Dm1, D2, Dm2 in gene_changes:
                
                    record = gene_name
                 
                    if gene_name in shared_genes:
                        continue
                 
                    if record not in gene_change_map:
                        gene_change_map[record] = set()
                    
                    gene_change_map[record].add(current_epoch)    
            
        # Done looping over samples
        # print record for species
        
        output_str = "%s, snps" % species_name
        for record in sorted(snp_change_map):
        
            contig, position, gene_name, variant_type = record
        
            epochs = [parse_timecourse_data.epochs[i] for i in sorted(snp_change_map[record]['epochs'])]
        
            snp_change_map[record]['initial_hs'] = numpy.array(snp_change_map[record]['initial_hs'])
            median_h = numpy.median(snp_change_map[record]['initial_hs'])
        
            record_str = ";".join([contig, str(position), gene_name, variant_type, str(median_h)]+epochs)
        
            output_str += (", "+record_str)
        
        file.write(output_str)
        file.write("\n")    
    
        output_str = "%s, genes" % species_name
        for record in sorted(gene_change_map):
        
            gene_name = record
        
            epochs = [parse_timecourse_data.epochs[i] for i in sorted(gene_change_map[record])]
        
        
            record_str = ";".join([gene_name]+epochs)
        
            output_str += (", "+record_str)
        
        file.write(output_str)
        file.write("\n")    
        
    file.close()         
    