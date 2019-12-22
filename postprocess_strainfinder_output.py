from StrainFinder import *
import cPickle
import sys
import numpy
import parse_midas_data
import parse_timecourse_data
from math import fabs

species_name = sys.argv[1]
Nmax = 20
output_directory = os.path.expanduser("~/strainfinder_output/")

filename_prefix = "%s/%s" % (output_directory, species_name)

# Get best (min) AIC version
# Get filenames of EM objects
fns = [filename_prefix+(".%d.em.cpickle" % N) for N in range(2,Nmax+1)]    

species_name+".strainfinder.locations.p"

# Load location object
snp_locations = cPickle.load(open(filename_prefix+".strainfinder.locations.p",'rb'))

# Load EM objects
ems = [cPickle.load(open(fn, 'rb')) for fn in fns]

# Get the best AIC in each EM object
aics = numpy.array([em.select_best_estimates(1)[0].aic for em in ems])

print "AIC scores:", aics.min()-aics

# Select EM with the minimum AIC
em = ems[numpy.argmin(aics)]
#em = ems[1]

# M = timepoint, L = along genome, 4 = alleles
input_alignment = em.data.x # alignment data, dim = (M x L x 4)

best_estimate = em.select_best_estimates(1)[0]

strain_genotypes = best_estimate.p # genotypes of first estimate, dim = (N x L x 4)
strain_freqs = best_estimate.z # (M x N)

#output_alignment = numpy.dot(strain_freqs, strain_genotypes)

strain_genotypes_derived = strain_genotypes[:,:,0]
strain_genotypes_ancestral = strain_genotypes[:,:,1]

strain_counts = strain_genotypes_derived.sum(axis=0)
strain_counts = numpy.fmin(strain_counts,strain_genotypes.shape[0]-strain_counts)
print "SFS strains"

singletons = (strain_counts==1)

for k in xrange(0,strain_genotypes.shape[0]):
    print k, (strain_counts==k).sum()

output_alignment = numpy.einsum('ij,jkl', strain_freqs, strain_genotypes)

print input_alignment.shape
print output_alignment.shape

print "Best # of strains:", em.select_best_estimates(1)[0].z.shape[1]
print em.select_best_estimates(1)[0].z # frequencies of first estimate, dim = (M x N)

sample_time_map = parse_timecourse_data.parse_sample_time_map()
ts = numpy.array([sample_time_map[sample] for sample in parse_timecourse_data.morteza_samples])

import pylab
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint
import matplotlib.colors as mcolors

mpl.rcParams['font.size'] = 6
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

pylab.figure(figsize=(7,6))
fig = pylab.gcf()

outer_grid  = gridspec.GridSpec(3,1, height_ratios=[1,1,1], hspace=0.15)

strain_axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(strain_axis)
strain_axis.set_ylabel('Strain frequencies')

strain_colors = []
for i in xrange(0,strain_freqs.shape[1]):
    line, = strain_axis.plot(ts,strain_freqs[:,i])
    strain_colors.append( pylab.getp(line,'color'))
strain_axis.set_ylim([0,1])

print strain_colors

predicted_axis = plt.Subplot(fig, outer_grid[1])
fig.add_subplot(predicted_axis)
predicted_axis.set_ylim([0,1])
predicted_axis.set_ylabel('Predicted\n SNV frequencies')


observed_axis = plt.Subplot(fig, outer_grid[2])
fig.add_subplot(observed_axis)
observed_axis.set_ylim([0,1])
observed_axis.set_ylabel('Observed\n SNV frequencies')
observed_axis.set_xlabel('Time (days)')

print output_alignment.shape
print input_alignment.shape

num_total = 0
num_flipped = 0

strain_singletons = numpy.array([0 for color in strain_colors])

for i in xrange(0,output_alignment.shape[1]):
    
    allele_counts = output_alignment[:,i,:]
    depths = allele_counts.sum(axis=1)
    predicted_freqs = (allele_counts[:,0]/(depths+(depths==0)))
    
    allele_counts = input_alignment[:,i,:]
    depths = allele_counts.sum(axis=1)
    observed_freqs = (allele_counts[:,0]/(depths+(depths==0)))
    
    good_idxs = (depths>0)
    
    if good_idxs.sum() < 2:
        continue
    
    # Polarize by predicted freq
    num_total+=1
    if predicted_freqs[0]>0.5:
        num_flipped+=1
        predicted_freqs = 1-predicted_freqs
        observed_freqs = 1-observed_freqs
    
    if singletons[i]:
        alt_allele = 1-long(numpy.median(strain_genotypes_derived[:,i]))
        #print alt_allele
        #print strain_genotypes_derived[:,i]
        strain_idx = (strain_genotypes_derived[:,i]==alt_allele).argmax()
        #print strain_idx
        color = strain_colors[strain_idx]
        strain_singletons[strain_idx]+=1
        zorder=2
    else:
        color = '0.7'
        zorder = 1
    predicted_axis.plot(ts, predicted_freqs,'-',color=color,alpha=0.5,zorder=zorder)
    observed_axis.plot(ts[good_idxs], observed_freqs[good_idxs],'-',color=color,alpha=0.5,zorder=zorder)

print num_flipped, num_total  
print strain_singletons
pylab.savefig('strainfinder_output.pdf',bbox_inches='tight')
pylab.savefig('strainfinder_output.png',bbox_inches='tight',dpi=300)


print "Calculating four gamete test"
total_tested = 0
num_haplotypes = numpy.array([0 for i in xrange(0,5)])
for i in xrange(0,strain_genotypes_derived.shape[1]):
    
    contig_i,location_i,polarization_i = snp_locations[i]
    #print contig_i,location_i
    #if singletons[i]:
    #    continue
    for j in xrange(i+1,strain_genotypes_derived.shape[1]):
    
        contig_j,location_j,polarization_j = snp_locations[j]
        
        if contig_i != contig_j:
            continue
            
        if fabs(location_i-location_j)>50000:
            continue
        
        #if singletons[j]:
        #    continue
        
        total_tested+=1
        
        haplotypes = 10*strain_genotypes_derived[:,i]+strain_genotypes_derived[:,j]
        
        unique_haplotypes = numpy.unique(haplotypes)
        
        
        num_haplotypes[len(unique_haplotypes)] += 1
        
print num_haplotypes
print num_haplotypes.sum()
print num_haplotypes*1.0/num_haplotypes.sum()  
        