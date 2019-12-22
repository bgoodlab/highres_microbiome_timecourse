import parse_midas_data
import config
import parse_timecourse_data
import pylab
import sfs_utils
import stats_utils
import numpy
import sys
import bz2

good_species_list = parse_midas_data.parse_good_species_list()
desired_samples = parse_timecourse_data.all_samples

for species_name in good_species_list:
        
    Dbar_map = parse_midas_data.parse_sample_coverage_map(species_name)
    
    samples = []
    Dbars = []
    for sample in desired_samples:
        if sample not in Dbar_map:
            Dbar=0
        else:
            Dbar = Dbar_map[sample]
    
        samples.append(sample)
        Dbars.append(Dbar)      
    
    filename = config.data_directory + ("/snps/%s/median_coverage.txt.bz2" % (species_name))
    # Median coverage file
    output_file = bz2.BZ2File(filename,"w")
    # Write header
    output_strs = list(samples)
    output_file.write("\t".join(output_strs))
    output_file.write("\n")
    output_file.write("\t".join(["%g" % Dbar for Dbar in Dbars]))
    output_file.close()
        
    
    
    