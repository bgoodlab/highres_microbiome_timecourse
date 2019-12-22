import os
import sys
import parse_midas_data
import parse_timecourse_data
import config
import numpy
import barcode_utils
import bz2
import core_gene_utils

base_directory = config.barcode_directory+"barcode_snps/"     
    
if __name__=='__main__':
    
    
    if len(sys.argv)==2:
        species_names = [sys.argv[1]]
    else:
        species_names = parse_midas_data.parse_good_species_list()

    # Make base directory
    os.system('mkdir -p %s' % base_directory)
    
    # Make output files for all samples (including Kuleshov one)
    all_samples = parse_timecourse_data.all_samples
 
    # Set of just the samples with barcode sequencing
    barcode_samples = set(parse_timecourse_data.morteza_samples)
 
    for species_name in species_names:
        
        sys.stderr.write("Loading core genes...\t")
        personal_core_genes = core_gene_utils.parse_personal_core_genes(species_name)
        sys.stderr.write("Done!\n")
        
        sys.stderr.write("Processing %s...\n" % species_name)
        desired_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name,allowed_samples=all_samples,allowed_genes=personal_core_genes)
        sys.stderr.write("Loaded SNVs!\n")
        
        # make list of SNPs and their starting keys
        sys.stderr.write("Creating SNV dictionary...\t")
        
        snp_info_strs = {}
        snp_As = {}
        snp_Ds = {}
        
        for gene_name in allele_counts_map.keys():
            for variant_type in allele_counts_map[gene_name].keys():
            
                allele_counts_matrix = allele_counts_map[gene_name][variant_type]['alleles']
                
                if len(allele_counts_matrix)==0:
                    continue
                    
                depths = allele_counts_matrix.sum(axis=2)
                alts = allele_counts_matrix[:,:,0]
            
                for site_idx in xrange(0,len(allele_counts_map[gene_name][variant_type]['locations'])):
                
                    chromosome, location = allele_counts_map[gene_name][variant_type]['locations'][site_idx]
                    alternate_allele = allele_counts_map[gene_name][variant_type]['polarizations'][site_idx]
                
                    snp = (chromosome, location)
                    info_str = "%s|%d|%s|%s|R|0" % (chromosome, location, gene_name, variant_type)
                
                    snp_info_strs[snp] = info_str
                    snp_As[snp] = numpy.array(alts[site_idx,:],dtype=numpy.int32)
                    snp_Ds[snp] = numpy.array(depths[site_idx,:],dtype=numpy.int32)
                    
                    if alternate_allele=='R':
                        # unpolarize 
                        # (since we'll be comparing to unpolarized barcode alleles
                        snp_As[snp] = snp_Ds[snp]-snp_As[snp]
        
        sys.stderr.write("Done!\n")
                
        for sample_idx in xrange(0,len(desired_samples)):
        
            sample_name = desired_samples[sample_idx]
            
            if sample_name in barcode_samples:
            
                sys.stderr.write("Loading barcodes from sample %s...\n" % sample_name)
                
                # First zero out entries for this sample
                for snp in snp_As:
                    snp_As[snp][sample_idx] = 0
                    snp_Ds[snp][sample_idx] = 0
            
                allele_barcode_map, allele_error_map = barcode_utils.parse_allele_barcode_tuples(species_name, sample_name, corrected=True)
        
                for allele_str in allele_barcode_map:
            
                    if not (allele_str[-1]=='A' or allele_str[-1]=='R'):
                        # Not a snp, ignore
                        continue
                
                    snp_items = allele_str.split("|")
                    chromosome = snp_items[0]
                    location = long(snp_items[1])
                    allele = snp_items[2]
                    snp = (chromosome,location)
            
                    if snp not in snp_info_strs:
                        continue
            
                    num_barcodes = len(allele_barcode_map[allele_str])
            
            
                    if allele=='A':
                        snp_As[snp][sample_idx] += num_barcodes
            
                    snp_Ds[snp][sample_idx] += num_barcodes
            
        Dss = []    
        for snp in sorted(snp_info_strs):
            if snp_Ds[snp].sum() < 0.5:
                continue # Don't print something with no reads!
            Dss.append(snp_Ds[snp])
        Dss = numpy.array(Dss)
        
        # Now we need to write the SNPs to disk
        species_directory = base_directory+species_name
        os.system('mkdir -p %s' % species_directory)
        filename = "%s/annotated_snps.txt.bz2" % species_directory
        output_file = bz2.BZ2File(filename,"w")
        # Write header
        output_strs = ["site_id"] + list(desired_samples)
        output_file.write("\t".join(output_strs))
        output_file.write("\n")
        
        if len(Dss)==0:
            output_file.close()
            continue
        
        Dbars = []
        for t_idx in xrange(0,Dss.shape[1]):
            
            Ds = Dss[:,t_idx]*1.0
            if (Ds>0).sum()==0:
                Dbars.append(0)
            else:
                Dbars.append(numpy.median(Ds[Ds>0]*1.0))
                
        Dbars = numpy.array(Dbars)   
        
        Dbars = Dbars*(Dbars>=config.barcode_min_median_coverage)
        
        # Need to calculate SFS as well
        allowed_variant_type_list = ['1D','2D','3D','4D']
        allowed_variant_types = set(allowed_variant_type_list)  


        # We shouldn't be doing this for raw data 
        #samples = parse_midas_data.parse_merged_sample_names(items)
    
        site_map = [{} for sample in desired_samples]
        for sample_idx in xrange(0,len(desired_samples)):
            site_map[sample_idx] = {variant_type:{} for variant_type in allowed_variant_types}

        for snp in sorted(snp_info_strs):
            As = snp_As[snp]
            Ds = snp_Ds[snp]
            
            good_idxs = (Ds>=0.3*Dbars)*(Ds<=3*Dbars)*(Dbars>0)
           
            As = As*good_idxs
            Ds = Ds*good_idxs
           
            highcoverage_idxs = (Ds>=config.barcode_min_good_coverage)
            
            if highcoverage_idxs.any():
                good_idxs = good_idxs*(As[highcoverage_idxs]>=0.1*Ds[highcoverage_idxs]).any()
            else:
                good_idxs = good_idxs*0
           
            As = As*good_idxs
            Ds = Ds*good_idxs
            
            snp_info_str = snp_info_strs[snp]   
            variant_type = snp_info_str.split("|")[3] 
            for i in xrange(0,len(As)):
                site = (As[i],Ds[i])
                if site not in site_map[i][variant_type]:
                    site_map[i][variant_type][site] = [0,0.0]
                site_map[i][variant_type][site][0] += 1
                site_map[i][variant_type][site][1] += 0 
            
            if Ds.sum()<0.5:
                continue # Don't print something with no reads!
            
            output_strs = [snp_info_strs[snp]]+["%d,%d" % item for item in zip(As,Ds)]
            output_file.write("\t".join(output_strs))
            output_file.write("\n")
        
        output_file.close()
        
        # Median coverage file
        filename = "%s/median_coverage.txt.bz2" % species_directory
        output_file = bz2.BZ2File(filename,"w")
        # Write header
        output_strs = list(desired_samples)
        output_file.write("\t".join(output_strs))
        output_file.write("\n")
        output_file.write("\t".join(["%g" % Dbar for Dbar in Dbars]))
        output_file.close()
        
        # Now write SFS
        output_file = bz2.BZ2File("%s/within_sample_sfs.txt.bz2" % (species_directory),"w")
        output_file.write("\t".join(["SampleID", "variant_type", "D,A,count,reverse_count", "..."]))
        for sample_idx in xrange(0,len(desired_samples)):
            sample = desired_samples[sample_idx]
            for variant_type in allowed_variant_type_list:
                output_file.write("\n")
                output_file.write("\t".join([sample, variant_type]+["%d,%d,%d,%g" % (site[0],site[1],site_map[sample_idx][variant_type][site][0],site_map[sample_idx][variant_type][site][1]) for site in sorted(site_map[sample_idx][variant_type].keys())]))
        output_file.close()

sys.stderr.write("Done!\n")
            