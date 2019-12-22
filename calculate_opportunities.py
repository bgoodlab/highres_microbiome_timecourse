import parse_midas_data
import config
import parse_timecourse_data
import pylab
import sfs_utils
import stats_utils
import numpy
import sys
import bz2

species_name = sys.argv[1]

desired_samples = parse_timecourse_data.all_samples


    
    filename = config.data_directory + ("/snps/%s/median_coverage.txt.bz2" % (species_name))
    # Median coverage file
    output_file = bz2.BZ2File(filename,"w")
    # Write header
    output_strs = list(samples)
    output_file.write("\t".join(output_strs))
    output_file.write("\n")
    output_file.write("\t".join(["%g" % Dbar for Dbar in Dbars]))
    output_file.close()
        
    
    
    