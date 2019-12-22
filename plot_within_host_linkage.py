import pylab
import numpy
from numpy.random import normal, shuffle
import stats_utils
import sys
#filename = "linkage_output.txt"
#filename="uniformis_within_host_linkage.txt"
filename=sys.argv[1]
file = open(filename,"r")
locations = []
gene_names = []
freqs = []
prevalences = []
targets = []
for line in file:
    items = line.split()
    locations.append( long(items[0]) )
    gene_names.append( items[1] )
    freqs.append( float(items[2]) )
    prevalences.append( float(items[3]) )
    targets.append( float(items[4]) )
    
idxs = numpy.arange(0,len(locations))
locations = numpy.array(locations)
freqs = numpy.array(freqs)
prevalences = numpy.array(prevalences)
targets = numpy.array(targets)

gene_freqs_map = {}
gene_targets_map = {}
unique_gene_names = []
for idx in idxs:
    gene_name = gene_names[idx]
    if gene_name not in gene_freqs_map:
        gene_freqs_map[gene_name] = []
        gene_targets_map[gene_name] = targets[idx]
        unique_gene_names.append(gene_name)
    gene_freqs_map[gene_name].append(freqs[idx])
    

gene_polymorphisms = []
gene_divergences = []
gene_num_sites = []
for gene_name in unique_gene_names:
    
    gene_freqs_map[gene_name] = numpy.array(gene_freqs_map[gene_name])
    if gene_targets_map[gene_name]<300:
        continue
    gene_num_sites.append( gene_targets_map[gene_name]  )
    gene_polymorphisms.append( ((gene_freqs_map[gene_name]>=0.2)*(gene_freqs_map[gene_name]<=0.8)).sum()*1.0/gene_targets_map[gene_name] )
    gene_divergences.append( ((gene_freqs_map[gene_name]>=0.8)).sum()*1.0/gene_targets_map[gene_name] )

good_idxs = (freqs>0.2)
private_idxs = (prevalences<=0.2)*(freqs>0.2)

pylab.figure(1,figsize=(12,1))
pylab.plot(idxs[good_idxs], freqs[good_idxs],'.',markersize=2,color='k')
pylab.plot(idxs[private_idxs], freqs[private_idxs],'.',markersize=3,color='r')

pylab.ylim([0,1.05])
pylab.savefig('site_polymorphism.png',bbox_inches='tight',dpi=300)
pylab.figure(2,figsize=(36,1))
pylab.plot(locations, freqs,'.',markersize=2)
pylab.plot(locations[private_idxs], freqs[private_idxs],'.',markersize=3,color='r')

pylab.ylim([0,1.05])
pylab.savefig('location_polymorphism.png',bbox_inches='tight',dpi=300)
pylab.figure(3,figsize=(12,1))
pylab.plot(numpy.arange(0,len(gene_polymorphisms)), gene_polymorphisms,'.',markersize=2)
pylab.ylim([-0.05,0.1])
pylab.savefig('gene_polymorphism.png',bbox_inches='tight',dpi=300)

pylab.figure(4,figsize=(3.42,2))
clipped_gene_polymorphisms = numpy.clip(gene_polymorphisms,1e-03,0.99)
pylab.hist(numpy.log10(clipped_gene_polymorphisms))
pylab.savefig('gene_polymorphism_hist.png',bbox_inches='tight',dpi=300)

pylab.figure(5,figsize=(3.42,2))

gene_polymorphisms_x = clipped_gene_polymorphisms
permuted_gene_polymorphisms_x = numpy.array(clipped_gene_polymorphisms,copy=True)
shuffle(permuted_gene_polymorphisms_x)
gene_polymorphisms_y = numpy.roll(clipped_gene_polymorphisms,1)
gene_polymorphisms_z = numpy.roll(clipped_gene_polymorphisms,2)

double_zeros = []
permuted_double_zeros = []
for k in xrange(1,20):
    gene_polymorphisms_k = numpy.roll(gene_polymorphisms,k)
    permuted_gene_polymorphisms_k = numpy.roll(permuted_gene_polymorphisms_x,k)
    
    double_zeros.append(((gene_polymorphisms_x<=1e-03)*(gene_polymorphisms_k<=1e-03)).sum())
    permuted_double_zeros.append(((gene_polymorphisms_x<=1e-03)*(permuted_gene_polymorphisms_k<=1e-03)).sum())

print double_zeros  
print permuted_double_zeros  
 
pylab.loglog(gene_polymorphisms_x*(1+normal(0,0.05,size=len(gene_polymorphisms_x))), gene_polymorphisms_y*(1+normal(0,0.05,size=len(gene_polymorphisms_y))),'k.',markersize=2)
pylab.savefig('gene_diversity_correlation.png',bbox_inches='tight')
print (gene_polymorphisms_x<=1e-03).sum(), ((gene_polymorphisms_x<=1e-03)*(gene_polymorphisms_z<=1e-03)).sum(), (gene_polymorphisms_x<=1e-02).sum()**2/(len(gene_polymorphisms))

for gene_idx in xrange(0,len(gene_polymorphisms)):
    if gene_polymorphisms[gene_idx]<=1e-04:
        print "*",
    else:
        print " ",
    print unique_gene_names[gene_idx], gene_polymorphisms[gene_idx], gene_divergences[gene_idx], gene_num_sites[gene_idx]  
    
    
# Calculate RUNS!
runs = []
current_run = 0
gimmie=False
for p in gene_polymorphisms:
    if p<=1e-04:
        current_run+=1
    else:
        if current_run>1.5 and p<3e-03 and (gimmie==False):
            # keep it alive
            gimmie=True
        else:
            if current_run>0.5:
                runs.append(current_run)
                current_run=0
            gimmie=False
            
            

print runs

permuted_runs = []
current_run = 0
gimmie=False

permuted_polymorphisms = []
gene_polymorphisms_copy = numpy.array(gene_polymorphisms,copy=True)
num_bootstraps = 100
for n in xrange(0,num_bootstraps):
    shuffle(gene_polymorphisms_copy)
    permuted_polymorphisms.extend(gene_polymorphisms_copy)

for p in permuted_polymorphisms:

    if p<=1e-04:
        current_run+=1
    else:
        if current_run>1.5 and p<3e-03 and (gimmie==False):
            # keep it alive
            gimmie=True
        else:
            if current_run>0.5:
                permuted_runs.append(current_run)
                current_run=0
            gimmie=False
                    
pylab.figure(6,figsize=(3.42,2))

null_xs, null_ns = stats_utils.calculate_unnormalized_survival_from_vector(permuted_runs)

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(runs)
pylab.step(null_xs,null_ns/null_ns[0]*ns[0],color='0.7',where='pre')
pylab.step(xs,ns,where='pre')


pylab.semilogy([1],[1])
pylab.ylim([0.8,300])
pylab.xlim([0,max(runs)+1])
pylab.xlabel("Consecutive 'homozygous' genes, g")
pylab.ylabel('Number of runs >=g')
pylab.savefig('runs.pdf',bbox_inches='tight')
         