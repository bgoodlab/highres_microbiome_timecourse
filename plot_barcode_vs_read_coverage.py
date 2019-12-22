import parse_midas_data
import config
import parse_timecourse_data
import pylab
import sfs_utils
import stats_utils
import numpy
import sys

pylab.figure(1)
pylab.xlabel('Read coverage')
pylab.ylabel('Barcode coverage')
pylab.loglog([1e-01,1e04],[1e-01,1e04],'k:')
pylab.loglog([1e-01,1e04],[1e-01/2,1e04/2],'k:')
pylab.xlim([1e-01,1e04])
pylab.ylim([1e-01,1e04])
good_species_list = parse_midas_data.parse_good_species_list()

#good_species_list = ['Eubacterium_eligens_61678']
#good_species_list = ['Oscillibacter_sp_60799']
#good_species_list = ['Eubacterium_ventriosum_61474']
good_species_list = ['Bacteroides_stercoris_56735']
desired_samples = parse_timecourse_data.morteza_samples

snps_directory=config.barcode_directory+"/barcode_snps/"

ratios = []

for species_name in good_species_list:
        
    read_Dbars = parse_midas_data.parse_sample_coverage_map(species_name)
    barcode_Dbars = parse_midas_data.parse_barcode_sample_coverage_map(species_name)
     
    #print barcode_Dbars 
        
    #barcode_samples, barcode_sfs_map = parse_midas_data.parse_within_sample_sfs(species_name, allowed_variant_types=set(['1D','2D','3D','4D']),snps_directory=snps_directory) 
        
    #barcode_Dbars = {}
    #for sample in desired_samples:
    #    if sample in barcode_sfs_map: 
    #        xs,ns = sfs_utils.calculate_depth_distribution_from_sfs_map(barcode_sfs_map[sample])
    #        if len(xs)==0:
    #            Dbar=0
    #        else:
    #            Dbar = stats_utils.calculate_median_from_distribution(xs,ns)
    #    else:
    #        Dbar = 0
    #    barcode_Dbars[sample] = Dbar
        
    read_Dbar_vector = []
    barcode_Dbar_vector = []
        
    for sample in desired_samples:
            
        if sample not in read_Dbars:
            continue
                
        if sample not in barcode_Dbars:
            sys.stderr.write("No barcodes for %s in %s (Dbar=%g)\n" % (species_name,sample,read_Dbars[sample]))
            continue
        
        if (read_Dbars[sample]<config.min_median_coverage) and (barcode_Dbars[sample]<config.barcode_min_median_coverage):
            continue
            
        if barcode_Dbars[sample]<config.barcode_min_median_coverage:
            print species_name, sample, read_Dbars[sample], barcode_Dbars[sample]
                
        read_Dbar_vector.append(read_Dbars[sample])
        barcode_Dbar_vector.append(barcode_Dbars[sample])
            
    read_Dbar_vector = numpy.array(read_Dbar_vector)
    barcode_Dbar_vector = numpy.array(barcode_Dbar_vector)
    
    ratios.extend( (barcode_Dbar_vector*1.0/read_Dbar_vector)[read_Dbar_vector>=config.min_median_coverage] )
    
    read_Dbar_vector = numpy.clip(read_Dbar_vector,1e-01,1e04)
    barcode_Dbar_vector = numpy.clip(barcode_Dbar_vector,1e-01,1e04)
    
    pylab.plot(read_Dbar_vector, barcode_Dbar_vector,'bo',alpha=0.5)

pylab.savefig('read_vs_barcode_coverages.pdf',bbox_inches='tight')
pylab.figure(2)
pylab.hist(ratios,bins=20)
pylab.savefig('barcode_coverage_ratio_histogram.pdf',bbox_inches='tight')