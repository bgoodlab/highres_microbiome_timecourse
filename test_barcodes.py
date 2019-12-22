import matplotlib.gridspec as gridspec
import matplotlib as mpl
from numpy.random import binomial, random_sample, shuffle, multinomial, poisson
import bz2
import parse_midas_data
import parse_timecourse_data
import matplotlib
import matplotlib.pyplot as plt
import timecourse_utils
import parse_patric
#import calculate_preexisting_snps
import calculate_snp_prevalences
import cluster_utils
import core_gene_utils
import calculate_within_species_fixations
import numpy
import pylab
from math import exp
       
desired_samples = parse_timecourse_data.morteza_samples
 
species_coverage_matrix, samples, species = parse_midas_data.parse_global_marker_gene_coverages()
sample_time_map = parse_timecourse_data.parse_sample_time_map()
species_ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, samples)
samples = numpy.array(samples)[sample_idxs]
species_coverage_matrix = species_coverage_matrix[:,sample_idxs]
total_coverage = species_coverage_matrix.sum(axis=0)
species_freq_matrix = numpy.clip(species_coverage_matrix*1.0/total_coverage,0, 2)    


ps = species_freq_matrix[:,0]

theory_Ms = numpy.arange(1,40)
# Delta version
#conditional_Sbars = numpy.array([(1-numpy.power(1-ps, M)).sum() for M in theory_Ms])
# Poisson version
conditional_Sbars = numpy.array([(1-numpy.exp(-ps*M)).sum() for M in theory_Ms])

# make fake dataset
Mbar = 8
dbar = 50

num_barcodes = 1000000
Ms = []
Ss = []
Ds = []
for barcode_idx in xrange(0,num_barcodes):
    M = poisson(Mbar)
    ns = multinomial(M,ps)
    ds = poisson(dbar*ns)
    ds = ds[ds>0.5]
    barcode_ps = ds*1.0/ds.sum()
    entropy_S = exp(-(barcode_ps*numpy.log(barcode_ps)).sum())
    threshold_S = (ds>2.5).sum()
    raw_S = (ns>0.5).sum()
    #print raw_S, entropy_S
    S = threshold_S
    D = poisson(M*dbar)
    Ms.append(M)
    Ss.append(S)
    Ds.append(D)
    
Ms = numpy.array(Ms)
Ss = numpy.array(Ss)
Ds = numpy.array(Ds)

observed_Sbar = Ss.mean()

# Estimate Mhat with crude point estimate
Mhat = theory_Ms[numpy.square(conditional_Sbars-observed_Sbar).argmin()]
print "Mbar =", Mbar
print "Overall Mhat =", Mhat

# Do binning thing (even though we don't need to in this case)
Dbins = numpy.logspace(0,3,20)
binned_Ds = []
binned_Sbars = []
binned_Mbars = []
binned_Mhats = []
for bin_idx in xrange(1,len(Dbins)):
    
    good_barcodes = (Ds<Dbins[bin_idx])*(Ds>=Dbins[bin_idx-1])
    
    if good_barcodes.sum() > 0.5:
        
        binned_Sbar = Ss[good_barcodes].mean()
        binned_Mbar = Ms[good_barcodes].mean()    
        binned_Mhat = theory_Ms[numpy.square(conditional_Sbars-binned_Sbar).argmin()]

        binned_Ds.append(Dbins[bin_idx])
        binned_Sbars.append(binned_Sbar)
        binned_Mbars.append(binned_Mbar)
        binned_Mhats.append(binned_Mhat)

# Plot conditional avg S vs M
pylab.figure(1)
pylab.xlabel('M')
pylab.ylabel('S')
pylab.plot(theory_Ms, conditional_Sbars)
pylab.plot(theory_Ms, numpy.ones_like(theory_Ms)*Ss.mean(),label='All D')
pylab.legend(loc='lower right',frameon=False)



pylab.figure(2)
pylab.semilogx(binned_Ds,binned_Mbars,label='True')
pylab.semilogx(binned_Ds,binned_Mhats,label='Inferred')
pylab.legend(loc='lower right',frameon=False)
pylab.xlabel('binned D')
pylab.ylabel('Mavg|D')
pylab.savefig('barcode_test.pdf',bbox_inches='tight')


