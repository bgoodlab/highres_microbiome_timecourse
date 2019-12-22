import numpy
import pylab
import diversity_utils
import parse_timecourse_data
import parse_midas_data
import sys
import calculate_temporal_changes

sample_time_map = parse_timecourse_data.parse_sample_time_map()
ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, parse_timecourse_data.morteza_samples)
all_samples = numpy.array(parse_timecourse_data.morteza_samples)[sample_idxs]


desired_samples = parse_timecourse_data.morteza_samples

num_all = 0
num_five = 0
num_ten = 0
num_good = 0
        
good_species_list = parse_midas_data.parse_good_species_list()
#good_species_list = ['Phascolarctobacterium_sp_59817']
for species_name in good_species_list:
    
    # Only plot samples above a certain depth threshold that are confidently phaseable.
    highcoverage_samples = diversity_utils.calculate_highcoverage_samples(species_name)
    
    haploid_samples = diversity_utils.calculate_haploid_samples(species_name)
    
    if len(highcoverage_samples)>=5:
        num_five+=1
    
    if len(highcoverage_samples)>=10:
        num_ten+=1
        
    num_all+=1
    
    has_all_epochs = parse_timecourse_data.has_all_epochs(highcoverage_samples)
    
    if has_all_epochs:
        num_good+=1
                
        print species_name, len(haploid_samples), "/", len(highcoverage_samples), "/", len(desired_samples)
        
print num_all, num_five, num_ten, num_good