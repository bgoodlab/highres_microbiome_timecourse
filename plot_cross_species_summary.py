import matplotlib  
matplotlib.use('Agg') 
import config
import parse_midas_data
import parse_timecourse_data
import bacterial_phylogeny_utils
import diversity_utils
import sfs_utils
import pylab
import sys
import numpy
import hmp_utils
from math import log10, fabs, log
from scipy.stats import gaussian_kde
import figure_utils

import calculate_barcode_within_species_fixations as calculate_within_species_fixations
import calculate_barcode_temporal_changes
from numpy.random import normal

import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint

mpl.rcParams['font.size'] = 6
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

min_abundance = 3e-04
max_abundance = 1

min_species_abundance_threshold = 1e-03

epoch_colors = {'initial': '#8856a7', 'antibiotic': '#d7191c', 'final' : '#2b83ba', 'disease': '#fdae61'}

initial_color = '#8856a7' 
antibiotic_color = '#d7191c'
final_color = '#2b83ba'
disease_color = '#fdae61'

PLOT_DISEASE_DIFFS = True

NUM_SPECIES_GENETIC_CHANGES_CHECKED = 0
NUM_SPECIES_GENETIC_CHANGES_ABX = 0

abcd_x_loc=-0.12

good_species_list = parse_midas_data.parse_super_good_species_list()

good_species_pretty_list = figure_utils.get_pretty_species_names(good_species_list)

pretty_species_name_map = {species_name: pretty_species_name for species_name, pretty_species_name in zip(good_species_list, good_species_pretty_list)}

species_abundance_distribution_map = hmp_utils.parse_species_abundance_distributions()

species_coverage_matrix, samples, species = parse_midas_data.parse_global_marker_gene_coverages()
sample_time_map = parse_timecourse_data.parse_sample_time_map()
ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, samples)
samples = numpy.array(samples)[sample_idxs]
species_coverage_matrix = species_coverage_matrix[:,sample_idxs]
total_coverage = species_coverage_matrix.sum(axis=0)
species_freq_matrix = numpy.clip(species_coverage_matrix*1.0/total_coverage,0, 2)    

initial_idx, antibiotic_idx, final_idx = parse_timecourse_data.get_initial_antibiotic_final_idxs(samples)

sample_list = list(samples)


initial_idx_1 = sample_list.index(parse_timecourse_data.highcoverage_start_1)
initial_idx_2 = sample_list.index(parse_timecourse_data.highcoverage_start_2)
antibiotic_idx = sample_list.index(parse_timecourse_data.highcoverage_antibiotic)
final_idx = sample_list.index(parse_timecourse_data.highcoverage_end)
disease_idx = sample_list.index(parse_timecourse_data.highcoverage_lyme)

initial_idxs = []
for sample_name in parse_timecourse_data.epoch_samples['initial']:
	initial_idxs.append(sample_list.index(sample_name))
initial_idxs = numpy.array(initial_idxs)

abx_idxs = []
for sample_name in parse_timecourse_data.epoch_samples['antibiotic']:
	abx_idxs.append(sample_list.index(sample_name))
	break
abx_idxs = numpy.array(abx_idxs)

print parse_timecourse_data.epoch_samples['antibiotic'], abx_idxs

initial_species_freqs = numpy.median(species_freq_matrix[:,initial_idxs],axis=1)
abx_species_freqs = species_freq_matrix[:,abx_idxs].min(axis=1)

#initial_species_freqs = numpy.fmax(species_freq_matrix[:,initial_idx_1], species_freq_matrix[:,initial_idx_2])

antibiotic_species_freqs = species_freq_matrix[:,antibiotic_idx]
final_species_freqs = species_freq_matrix[:,final_idx]
disease_species_freqs = species_freq_matrix[:,disease_idx]

antibiotic_fold_changes = numpy.clip(antibiotic_species_freqs/(initial_species_freqs+(initial_species_freqs<1e-07)),1.3e-02,1e02/1.3)
final_fold_changes = numpy.clip(final_species_freqs/(initial_species_freqs+(initial_species_freqs<1e-07)),1.3e-02,1e02/1.3)

print "Calculating species level abundance changes"
print (initial_species_freqs>=1e-03).sum(), "species with median freq over initial timepoints >= .1%"
print ((abx_species_freqs<=(0.1*initial_species_freqs))*(initial_species_freqs>=1e-03)).sum(), "species with 10-fold reduction at ABX"
print ((abx_species_freqs<=(0.5*initial_species_freqs))*(initial_species_freqs>=1e-03)).sum(), "species with 2-fold reduction at ABX"




good_idxs = []
for species_idx in xrange(0,len(species)):
    if species[species_idx] in good_species_list:
        good_idxs.append(True)
    else:
        good_idxs.append(False)
        

#good_idxs = initial_species_freqs>=1e-04
#good_idxs = species_freq_matrix.mean(axis=1)>=min_species_abundance_threshold
species = numpy.array(species)[good_idxs]
initial_species_freqs = initial_species_freqs[good_idxs]
antibiotic_species_freqs = antibiotic_species_freqs[good_idxs]
disease_species_freqs = disease_species_freqs[good_idxs]
final_species_freqs = final_species_freqs[good_idxs]
antibiotic_fold_changes = antibiotic_fold_changes[good_idxs]
final_fold_changes = final_fold_changes[good_idxs]
avg_freqs = numpy.median(species_freq_matrix[good_idxs],axis=1)

print "Big changes: ", len(antibiotic_fold_changes), (antibiotic_fold_changes<0.1).sum()


species_freq_trajectories = numpy.clip(species_freq_matrix[good_idxs,:],4e-04,max_abundance)
#species_freq_trajectories = [species_freq_trajectories[idx,:] for idx in xrange(0,species_freq_trajectories.shape[0])]

species_idxs = range(0,len(species))

######
#
# Do computations
#
######
abundance_data = [{} for species_idx in species_idxs]
polymorphism_data = [{} for species_idx in species_idxs]
fixation_data = [{'disease': -1, 'antibiotic':-1, 'final':-1, 'all':-1} for species_idx in species_idxs]
initial_freq_data = [{'disease':[], 'antibiotic': [], 'final': []} for species_idx in species_idxs]
final_freq_data = [{'disease':[], 'antibiotic': [], 'final': []} for species_idx in species_idxs]


retention_data = [{'disease':-1, 'antibiotic':-1, 'final': -1} for species_idx in species_idxs]

for species_idx in species_idxs:
    species_name = species[species_idx]

    if species_name not in good_species_list:
        continue


    if species_freq_trajectories[species_idx, initial_idx_2]>species_freq_trajectories[species_idx, initial_idx_1]:
        initial_idx = initial_idx_2
    else:
        initial_idx = initial_idx_1

    for sample_idx in xrange(0,len(sample_list)):
    
        sample = sample_list[sample_idx]
        
        if sample not in parse_timecourse_data.morteza_samples:
            continue
        
        abundance = species_freq_trajectories[species_idx,sample_idx]
            
        epoch = parse_timecourse_data.sample_epoch_map[sample]
        
        if epoch == 'postantibiotic':
            epoch = 'antibiotic'
        
        if epoch == 'prefinal':
            epoch = 'final'
            
        if epoch not in abundance_data[species_idx]:
            abundance_data[species_idx][epoch] = []
            
        abundance_data[species_idx][epoch].append(abundance)
            
    for epoch in abundance_data[species_idx]:    
        abundance_data[species_idx][epoch].sort()
    
    sys.stderr.write("Loading SFSs for %s...\t" % species_name)
    sfs_samples, sfs_map = parse_midas_data.parse_within_sample_sfs(species_name) 
    sample_coverage_map = parse_midas_data.parse_sample_coverage_map(species_name)
    sys.stderr.write("Done!\n")
    
    polymorphism_data[species_idx]['max'] = 0
    for sample in parse_timecourse_data.morteza_samples:
        
        if sample not in sfs_map:
            continue
            
        if sample not in sample_coverage_map:
            continue
            
        Dbar = sample_coverage_map[sample]
        
        if Dbar < 10:
            continue
            
        ### Add part about coverage
            
        within_sites, between_sites, total_sites = sfs_utils.calculate_polymorphism_rates_from_sfs_map(sfs_map[sample])
        
        within_polymorphism = within_sites*1.0/total_sites
        if within_polymorphism<1e-05:
            within_polymorphism=1e-05
        
        epoch = parse_timecourse_data.sample_epoch_map[sample]
        
        if epoch == 'postantibiotic':
            epoch = 'antibiotic'
        
        if epoch == 'prefinal':
            epoch = 'final'
            
        if epoch not in polymorphism_data[species_idx]:
            polymorphism_data[species_idx][epoch] = []
            
        polymorphism_data[species_idx][epoch].append(within_polymorphism*1e06)
        polymorphism_data[species_idx]['max'] = max([within_polymorphism*1e06, polymorphism_data[species_idx]['max']])
            
    for epoch in polymorphism_data[species_idx]:
        if epoch!='max':
            polymorphism_data[species_idx][epoch].sort()
    print "Max polymorphism", polymorphism_data[species_idx]['max']
            
    # Put it here!
    # Now load and plot Fst data
    temporal_change_map = calculate_barcode_temporal_changes.load_temporal_change_map(species_name)

    n_tested, n_observed, n_expected, snp_changes = calculate_barcode_temporal_changes.calculate_mutations_from_temporal_change_map(temporal_change_map, 'any', 'any',type='snps')

    fixation_data[species_idx]['all'] = n_tested

    disease_snps, disease_genes = calculate_within_species_fixations.load_within_species_fixations(species_name,allowed_epochs=['initial_disease'])
    antibiotic_snps, antibiotic_genes = calculate_within_species_fixations.load_within_species_fixations(species_name,allowed_epochs=['initial_antibiotic','initial_postantibiotic'])
    final_snps, final_genes = calculate_within_species_fixations.load_within_species_fixations(species_name,allowed_epochs=['initial_final'])
    
    retained_antibiotic_snps = (antibiotic_snps & final_snps)
        
    plot_disease = (n_tested>0.5)            
    plot_antibiotic = (n_tested>0.5) 
    plot_final = (n_tested>0.5) 
    
    if not (plot_antibiotic or plot_final):
        fixation_data[species_idx]['all'] = 0    
    
    if plot_disease:
        
        fixation_data[species_idx]['disease'] = len(disease_snps)
        
        if len(disease_snps)>0:
            initial_freqs = []
            final_freqs = []
            for antibiotic_snp in disease_snps:
                contig, location, gene_name, variant_type, initial_freq, final_freq = antibiotic_snp
                initial_freqs.append(initial_freq)
                final_freqs.append(final_freq)
                
            initial_freqs.sort()
            initial_freqs = numpy.array(initial_freqs)
        
            # transform to frequencies
            #initial_hs = (1-numpy.sqrt(1-2*initial_hs))/2
        
            initial_freq_data[species_idx]['disease'] = initial_freqs
    
            final_freqs.sort()
            final_freqs = numpy.array(final_freqs)
        
            # transform to frequencies
            #initial_hs = (1-numpy.sqrt(1-2*initial_hs))/2
        
            final_freq_data[species_idx]['disease'] = final_freqs
            retention_data[species_idx]['disease'] = (final_freqs>=0.7).sum()*1.0/len(final_freqs)
            
    if plot_antibiotic:
        
        fixation_data[species_idx]['antibiotic'] = len(antibiotic_snps)
        
        if len(antibiotic_snps)>0:
            initial_freqs = []
            final_freqs = []
            for antibiotic_snp in antibiotic_snps:
                contig, location, gene_name, variant_type, initial_freq, final_freq = antibiotic_snp
                initial_freqs.append(initial_freq)
                final_freqs.append(final_freq)
                
            initial_freqs.sort()
            initial_freqs = numpy.array(initial_freqs)
        
            # transform to frequencies
            #initial_hs = (1-numpy.sqrt(1-2*initial_hs))/2
        
            initial_freq_data[species_idx]['antibiotic'] = initial_freqs
    
            final_freqs.sort()
            final_freqs = numpy.array(final_freqs)
        
            # transform to frequencies
            #initial_hs = (1-numpy.sqrt(1-2*initial_hs))/2
        
            final_freq_data[species_idx]['antibiotic'] = final_freqs
            retention_data[species_idx]['antibiotic'] = (final_freqs>=0.7).sum()*1.0/len(final_freqs)
            
    if plot_final:
        
        fixation_data[species_idx]['final'] = len(final_snps)
        
        if len(final_snps)>0:
            initial_freqs = []
            final_freqs = []
            for antibiotic_snp in final_snps:
                contig, location, gene_name, variant_type, initial_freq, final_freq = antibiotic_snp
                initial_freqs.append(initial_freq)
                final_freqs.append(final_freq)
                
            initial_freqs.sort()
            initial_freqs = numpy.array(initial_freqs)
        
            # transform to frequencies
            #initial_hs = (1-numpy.sqrt(1-2*initial_hs))/2
        
            initial_freq_data[species_idx]['final'] = initial_freqs
    
            final_freqs.sort()
            final_freqs = numpy.array(final_freqs)
        
            # transform to frequencies
            #initial_hs = (1-numpy.sqrt(1-2*initial_hs))/2
        
            final_freq_data[species_idx]['final'] = final_freqs
            retention_data[species_idx]['final'] = (final_freqs>=0.7).sum()*1.0/len(final_freqs)
            

def polymorphism_category(num_changes,allow_negatives=True):
    if num_changes<7000:
        fixation_category = 0
    else:
        fixation_category = 1
    
    return fixation_category

    
def possible_fixation_category(num_changes,allow_negatives=True):
    if num_changes<0.5:
        fixation_category = 0
    elif num_changes<7e02:
        fixation_category = 1
    else:
        fixation_category = 3
        
    return fixation_category

def fixation_category(num_changes, allow_negatives=True):
    
    if num_changes<0.5:
        return 0
    elif num_changes>1e04:
        return 1e04
    else:
        return num_changes
    
    if num_changes<0 and allow_negatives:
        fixation_category = 0
    elif num_changes<0.5:
        fixation_category = 1
    elif num_changes<50:
        fixation_category = 2
    elif num_changes<1000:
        fixation_category = 3
    else:
        fixation_category = 4
        
    return fixation_category

polymorphism_categories = [polymorphism_category(polymorphism_data[species_idx]['max']) for species_idx in species_idxs]
        
possible_fixation_categories = [possible_fixation_category(fixation_data[species_idx]['all']) for species_idx in species_idxs]

antibiotic_fixation_categories = [fixation_category(fixation_data[species_idx]['antibiotic'],allow_negatives=False) for species_idx in species_idxs]
final_fixation_categories = [fixation_category(fixation_data[species_idx]['final']) for species_idx in species_idxs]


### Sort and plot!

species_idxs = numpy.array(sorted(species_idxs,key=lambda x: (polymorphism_categories[x],antibiotic_fixation_categories[x], final_fixation_categories[x]), reverse=True))

# sort everything by descending order of XXX
#initial_species_freqs, species, antibiotic_fold_changes, final_fold_changes, species_freq_trajectories, antibiotic_species_freqs, final_species_freqs = (numpy.array(x) for x in zip(*sorted(zip(initial_species_freqs, species, antibiotic_fold_changes, final_fold_changes, species_freq_trajectories, antibiotic_species_freqs, final_species_freqs), key=lambda pair: (-pair[2],-pair[3]), reverse=True)))

#initial_species_freqs, species, antibiotic_fold_changes, final_fold_changes = (numpy.array(x) for x in zip(*sorted(zip(initial_species_freqs, species, antibiotic_fold_changes, final_fold_changes), key=lambda pair: (-pair[0],-pair[2]), reverse=True)))

pylab.figure(1,figsize=(7,7))
fig = pylab.gcf()

outer_grid  = gridspec.GridSpec(7,1,height_ratios=[0.05, 0.7,0.9,0.3, 1,0.5, 0.5],hspace=0.1)

legend_axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(legend_axis)

legend_axis.set_ylim([0,1])
legend_axis.set_xlim([0,1])

legend_axis.spines['top'].set_visible(False)
legend_axis.spines['right'].set_visible(False)
legend_axis.spines['left'].set_visible(False)
legend_axis.spines['bottom'].set_visible(False)

legend_axis.set_xticks([])
legend_axis.set_yticks([])

legend_axis.plot([-2],[-1],'.-',markersize=3,zorder=1,markeredgewidth=0,label='Baseline',color=initial_color)
legend_axis.plot([-2],[-1],'.-',markersize=3,zorder=3,markeredgewidth=0,color=disease_color,label='Disease') 
legend_axis.plot([-2],[-1],'.-',markersize=3,zorder=3,markeredgewidth=0,label='Antibiotics (ABX)',color=antibiotic_color)
legend_axis.plot([-2],[-1],'.-',markersize=3,zorder=2,markeredgewidth=0,label='Final',color=final_color)

#legend_axis.legend(loc=(0,-0.1),frameon=False,fontsize=7,numpoints=1,ncol=4,handlelength=1)

legend_axis.legend(loc='lower right',frameon=False,fontsize=7,numpoints=1,ncol=5,handlelength=1)   


abundance_axis = plt.Subplot(fig, outer_grid[1])
fig.add_subplot(abundance_axis)
abundance_axis.set_ylabel('Relative\nabundance',fontsize=7)

abundance_axis.text(abcd_x_loc, 1, 'a', horizontalalignment='center', verticalalignment='center', transform=abundance_axis.transAxes,fontweight='bold',fontsize=7)

polymorphism_axis = plt.Subplot(fig, outer_grid[2])
fig.add_subplot(polymorphism_axis)
polymorphism_axis.set_ylabel('Within-species\ndiversity (SNVs/Mbp)',fontsize=7)

polymorphism_axis.text(abcd_x_loc, 1, 'b', horizontalalignment='center', verticalalignment='center', transform=polymorphism_axis.transAxes,fontweight='bold',fontsize=7)

#fst_axis = plt.Subplot(fig, outer_grid[3])
#fig.add_subplot(fst_axis)
#fst_axis.set_ylabel('Fst',fontsize=7)
#fst_axis.set_xlim([1e-01,1e05])
#fst_axis.text(abcd_x_loc, 1, 'b', horizontalalignment='center', verticalalignment='center', transform=fst_axis.transAxes,fontweight='bold',fontsize=7)

legend2_axis = plt.Subplot(fig, outer_grid[3])
fig.add_subplot(legend2_axis)

legend2_axis.set_ylim([0,1])
legend2_axis.set_xlim([0,1])

legend2_axis.spines['top'].set_visible(False)
legend2_axis.spines['right'].set_visible(False)
legend2_axis.spines['left'].set_visible(False)
legend2_axis.spines['bottom'].set_visible(False)

legend2_axis.set_xticks([])
legend2_axis.set_yticks([])

legend2_axis.plot([-2],[-1],'s',markersize=3,zorder=2,markeredgewidth=0,label='SNV differences between baseline & disease',color=disease_color)
legend2_axis.plot([-2],[-1],'s',markersize=3,zorder=3,label='baseline & ABX',color=antibiotic_color,markeredgewidth=0)
legend2_axis.plot([-2],[-1],'s',markersize=3,zorder=2,markeredgewidth=0,label='baseline & final',color=final_color)

legend2_axis.legend(loc=(0,-0.2),frameon=False,fontsize=7,numpoints=1,ncol=3,handlelength=1)   


fixation_axis = plt.Subplot(fig, outer_grid[4])
fig.add_subplot(fixation_axis)
fixation_axis.set_ylabel('Total number of\n SNV differences',fontsize=7)

fixation_axis.text(abcd_x_loc, 1, 'c', horizontalalignment='center', verticalalignment='center', transform=fixation_axis.transAxes,fontweight='bold',fontsize=7)


initial_axis = plt.Subplot(fig, outer_grid[5])
fig.add_subplot(initial_axis)
initial_axis.set_ylabel("Initial\nfreqs", fontsize=7)

initial_axis.text(abcd_x_loc, 1, 'd', horizontalalignment='center', verticalalignment='center', transform=initial_axis.transAxes,fontweight='bold',fontsize=7)


retention_axis = plt.Subplot(fig, outer_grid[6])
fig.add_subplot(retention_axis)
retention_axis.set_ylabel("Fraction \n retained", fontsize=7)

retention_axis.text(abcd_x_loc, 1, '(e)', horizontalalignment='center', verticalalignment='center', transform=retention_axis.transAxes,fontweight='bold',fontsize=7)




xticks = []
xticklabels = []

for x in xrange(0,len(species_idxs)):
    species_idx = species_idxs[x]
    species_name = species[species_idx]

    if species_freq_trajectories[species_idx, initial_idx_2]>species_freq_trajectories[species_idx, initial_idx_1]:
        initial_idx = initial_idx_2
    else:
        initial_idx = initial_idx_1
    
    if x in [15]:
        lower_linewidth=0.75
        lower_color = 'k'
    else:
        lower_linewidth=0.25
        lower_color = '0.9'
        
    if x in []: #[len(species_idxs)-1]:
        upper_linewidth=0.75
        upper_color = 'k'
    else:
        upper_linewidth=0.25
        upper_color = '0.9'
     
    # Plot grid
    abundance_axis.plot([x-0.5,x-0.5], [1e-04,1],'-',linewidth=lower_linewidth, color=lower_color)
    abundance_axis.plot([x+0.5,x+0.5], [1e-04,1],'-',linewidth=upper_linewidth,color=upper_color)
    
    polymorphism_axis.plot([x-0.5,x-0.5], [1e-01,1e06],'-',linewidth=lower_linewidth,color=lower_color)
    polymorphism_axis.plot([x+0.5,x+0.5], [1e-01,1e06],'-',linewidth=upper_linewidth,color=upper_color)
    
    #fst_axis.plot([x-0.5,x-0.5], [1e-01,1e05],'-',linewidth=lower_linewidth,color=lower_color)
    #fst_axis.plot([x+0.5,x+0.5], [1e-01,1e05],'-',linewidth=upper_linewidth,color=upper_color)
    
    
    fixation_axis.plot([x-0.5, x-0.5], [3e-01,1e05],'-',linewidth=lower_linewidth,color=lower_color)
    fixation_axis.plot([x+0.5, x+0.5], [3e-01,1e05],'-',linewidth=upper_linewidth,color=upper_color)
    
    initial_axis.plot([x-0.5, x-0.5], [-0.5,1.1],'-',linewidth=lower_linewidth,color=lower_color)
    initial_axis.plot([x+0.5, x+0.5], [-0.5,1.1],'-',linewidth=upper_linewidth,color=upper_color)
    
    
    retention_axis.plot([x-0.5,x-0.5], [-0.5,1.1],'-',linewidth=lower_linewidth,color=lower_color)
    retention_axis.plot([x+0.5,x+0.5], [-0.5,1.1],'-',linewidth=upper_linewidth,color=upper_color)
    
    
    # First plot abundances
    f0 = initial_species_freqs[species_idx]
    fa = max([min_abundance, antibiotic_species_freqs[species_idx]])
    ff = max([min_abundance, final_species_freqs[species_idx]])
    fd = max([min_abundance, disease_species_freqs[species_idx]])
    
    species_items = species_name.split("_")
    #label = "%s %s (%s)" % (species_items[0], species_items[1], species_items[2])
    label = pretty_species_name_map[species_name]
    
    change_str = ''
    
    if antibiotic_species_freqs[species_idx]<0.1*initial_species_freqs[species_idx]:
        change_str = "*"+change_str
        if final_species_freqs[species_idx]<0.1*initial_species_freqs[species_idx]:
            change_str = "*"+change_str
    xticks.append(x)
    xticklabels.append(label)
    
    
    #print species_name, f0, fa, ff, species_idx
    f_trajectory = species_freq_trajectories[species_idx,:]
    
    hmp_fs = species_abundance_distribution_map[species_name]
        
    abundance_axis.text(x,1.05,change_str,ha='center',fontsize=5)
    
    good_idxs = (hmp_fs>1e-04)
    if good_idxs.sum()>2:
    
        hmp_fs = hmp_fs[good_idxs]
        log_hmp_fs = numpy.log(hmp_fs)
    
        
        kernel = gaussian_kde(log_hmp_fs)
    
        theory_fs = numpy.logspace(-4,0,100)
        theory_log_fs = numpy.log(theory_fs)
        theory_pdf = kernel(theory_log_fs)
        theory_pdf = theory_pdf / theory_pdf.max() * 0.35
    
        abundance_axis.fill_betweenx(theory_fs,x-theory_pdf, x+theory_pdf,linewidth=0,facecolor='0.9',zorder=0)    
    
    if 'initial' in abundance_data[species_idx]:
        
        fs = abundance_data[species_idx]['initial']
        abundance_axis.semilogy([x-0.3]*len(fs), fs,'.-',color=epoch_colors['initial'],markersize=2,linewidth=0.25)
    
    if 'disease' in abundance_data[species_idx]:
        
        fs = abundance_data[species_idx]['disease']
        abundance_axis.semilogy([x-0.1]*len(fs), fs,'.-',color=epoch_colors['disease'],markersize=2,linewidth=0.25)
    
    if 'antibiotic' in abundance_data[species_idx]:
        
        fs = abundance_data[species_idx]['antibiotic']
        abundance_axis.semilogy([x+0.1]*len(fs), fs,'.-',color=epoch_colors['antibiotic'],markersize=2,linewidth=0.5)
    
    if 'final' in abundance_data[species_idx]:
        
        fs = abundance_data[species_idx]['final']
        abundance_axis.semilogy([x+0.3]*len(fs), fs,'.-',color=epoch_colors['final'],markersize=2,linewidth=0.25)
        
    # Now plot polymorphisms 
    
    if 'initial' in polymorphism_data[species_idx]:
        
        pis = polymorphism_data[species_idx]['initial']
        polymorphism_axis.semilogy([x-0.3]*len(pis), pis,'.-',color=epoch_colors['initial'],markersize=2,linewidth=0.25)
    
    if 'disease' in polymorphism_data[species_idx]:
        
        pis = polymorphism_data[species_idx]['disease']
        polymorphism_axis.semilogy([x-0.1]*len(pis), pis,'.-',color=epoch_colors['disease'],markersize=2,linewidth=0.25)
    
    if 'antibiotic' in polymorphism_data[species_idx]:
        
        pis = polymorphism_data[species_idx]['antibiotic']
        polymorphism_axis.semilogy([x+0.1]*len(pis), pis,'.-',color=epoch_colors['antibiotic'],markersize=2,linewidth=0.5)
    
    if 'final' in polymorphism_data[species_idx]:
        
        pis = polymorphism_data[species_idx]['final']
        polymorphism_axis.semilogy([x+0.3]*len(pis), pis,'.-',color=epoch_colors['final'],markersize=2,linewidth=0.25)
        
      
    if fixation_data[species_idx]['all']<0.5:
        fixation_axis.fill_between([x-0.5,x+0.5], [3e-01,3e-01],[1e05,1e05],color='0.9')
        initial_axis.fill_between([x-0.5,x+0.5], [-0.1,-0.1],[1.1,1.1],color='0.9')  
        retention_axis.fill_between([x-0.5,x+0.5], [-0.1,-0.1],[1.1,1.1],color='0.9')  
        continue
    else:
        fixation_axis.fill_between([x-0.5,x+0.5], [fixation_data[species_idx]['all'], fixation_data[species_idx]['all']],[1e05,1e05],color='0.9')      
            
    # Now plot fixation data
    if (fixation_data[species_idx]['antibiotic']<-0.1) and (fixation_data[species_idx]['final']<-0.1):
        fixation_axis.fill_between([x-0.5,x+0.5], [3e-01,3e-01],[3e04,3e04],color='0.9')
        initial_axis.fill_between([x-0.5,x+0.5], [-0.1,-0.1],[1.1,1.1],color='0.9')  
        retention_axis.fill_between([x-0.5,x+0.5], [-0.1,-0.1],[1.1,1.1],color='0.9')  
        continue
       
    #if (fixation_data[species_idx]['all']>-0.1):
    #    fixation_axis.semilogy([x], numpy.clip([fixation_data[species_idx]['all']],5e-01, 1e06),'k_',markersize=2,zorder=3) 
    
    if PLOT_DISEASE_DIFFS:
        if (fixation_data[species_idx]['disease']>0.5):
        #fixation_axis.semilogy([x-0.3,x-0.3], numpy.clip([0,fixation_data[species_idx]['disease']],5e-02, 1e06),'o-',markersize=2,zorder=3,color=disease_color,linewidth=0.25) 
            fixation_axis.bar([x-0.3], [fixation_data[species_idx]['disease']],bottom=[1e-02],width=0.2,facecolor=disease_color,linewidth=0)
        
        NUM_SPECIES_GENETIC_CHANGES_CHECKED += 1
        if (fixation_data[species_idx]['antibiotic']>0.5):
            #fixation_axis.semilogy([x,x], numpy.clip([0,fixation_data[species_idx]['antibiotic']],5e-02, 1e06),'o-',markersize=2,zorder=3,color=antibiotic_color,linewidth=0.25) 
    
            fixation_axis.bar([x-0.1], [fixation_data[species_idx]['antibiotic']],bottom=[1e-02],width=0.2,facecolor=antibiotic_color,linewidth=0)
            
            NUM_SPECIES_GENETIC_CHANGES_ABX += 1

    
        if (fixation_data[species_idx]['final']>0.5):
            #fixation_axis.semilogy([x+0.3,x+0.3], numpy.clip([0,fixation_data[species_idx]['final']],5e-02, 1e06),'o-',markersize=2,zorder=2,color=final_color,markeredgewidth=0,linewidth=0.25)  
            #print 'plotting', fixation_data[species_idx]['final']
            fixation_axis.bar([x+0.1], [fixation_data[species_idx]['final']],bottom=[1e-02], width=0.2,facecolor=final_color,linewidth=0)
     
    
    
     
        
    # Now plot initial freq data   
    if (len(initial_freq_data[species_idx]['antibiotic'])<0.5) and (len(initial_freq_data[species_idx]['final'])<0.5):
        initial_axis.fill_between([x-0.5,x+0.5], [-0.1,-0.1],[1.1,1.1],color='0.9')  
    
    if (len(initial_freq_data[species_idx]['disease'])>0.5):
        
        initial_hs = initial_freq_data[species_idx]['disease']
        
        if len(initial_hs)>10:
            # plot mean and stddev
            lower_h = initial_hs[long(0.25*len(initial_hs))]
            upper_h = initial_hs[long(0.75*len(initial_hs))]
            median_h = initial_hs[long(0.5*len(initial_hs))]
            
            initial_axis.plot([x-0.3,x-0.3], [lower_h, upper_h],'-',linewidth=0.5,zorder=1,color=disease_color)    
            initial_axis.plot([x-0.3], [median_h],'s',markersize=2,zorder=1,markeredgewidth=0,color=disease_color)    
        
        else:
            # just plot all the points
            for initial_h in initial_hs:
                initial_axis.plot([x-0.3], [initial_h],'o',markersize=1,zorder=1,markeredgewidth=0,color=disease_color)    
        
    if (len(initial_freq_data[species_idx]['antibiotic'])>0.5):
        
        initial_hs = initial_freq_data[species_idx]['antibiotic']
        
        if len(initial_hs)>10:
            # plot mean and stddev
            lower_h = initial_hs[long(0.25*len(initial_hs))]
            upper_h = initial_hs[long(0.75*len(initial_hs))]
            median_h = initial_hs[long(0.5*len(initial_hs))]
            
            initial_axis.plot([x,x], [lower_h, upper_h],'-',linewidth=0.5,zorder=1,color=antibiotic_color)    
            initial_axis.plot([x], [median_h],'s',markersize=2,zorder=1,markeredgewidth=0,color=antibiotic_color)    
        
        else:
            # just plot all the points
            for initial_h in initial_hs:
                initial_axis.plot([x], [initial_h],'o',markersize=1,zorder=1,markeredgewidth=0,color=antibiotic_color)    
    
    
    
    if (len(initial_freq_data[species_idx]['final'])>0.5):
        
        initial_hs = initial_freq_data[species_idx]['final']
        
        if len(initial_hs)>10:
            # plot mean and stddev
            lower_h = initial_hs[long(0.25*len(initial_hs))]
            upper_h = initial_hs[long(0.75*len(initial_hs))]
            median_h = initial_hs[long(0.5*len(initial_hs))]
            
            initial_axis.plot([x+0.3,x+0.3], [lower_h, upper_h],'-',linewidth=0.5,zorder=1, color=final_color)    
            initial_axis.plot([x+0.3], [median_h],'s',markersize=2,zorder=1,markeredgewidth=0, color=final_color)    
        
        else:
            # just plot all the points
            for initial_h in initial_hs:
                initial_axis.plot([x+0.3], [initial_h],'o',markersize=1,zorder=1,markeredgewidth=0, color=final_color)    
    
    # Now plot retention data
    
    if retention_data[species_idx]['disease']>-0.5:
        
        retention_axis.plot([x-0.3,x-0.3], [-1,retention_data[species_idx]['disease']],'o-',markersize=2,zorder=2,color=disease_color,markeredgewidth=0,linewidth=0.25)   
    
    if retention_data[species_idx]['antibiotic']>-0.5:
        
        retention_axis.plot([x,x], [-1,retention_data[species_idx]['antibiotic']],'o-',markersize=3,zorder=2,color=antibiotic_color,markeredgewidth=0,linewidth=0.5)   
    
    if False and retention_data[species_idx]['final']>-0.5:
        
        retention_axis.plot([x+0.3,x+0.3], [-1,retention_data[species_idx]['final']],'o-',markersize=2,zorder=2,color=final_color,markeredgewidth=0,linewidth=0.25)   
    
    
        
    
abundance_axis.set_ylim([min_abundance, max_abundance])
abundance_axis.fill_between([len(species_idxs)-0.5,len(species_idxs)+1],[min_abundance, min_abundance],[max_abundance, max_abundance],color='0.9')
abundance_axis.fill_between([-2,-0.5],[min_abundance, min_abundance],[max_abundance, max_abundance],color='0.9')

xmin = -1
xmax = len(species_idxs)

abundance_axis.set_xlim([xmin,xmax])
abundance_axis.set_xticks([])

polymorphism_axis.set_ylim([1e1,1e05])
polymorphism_axis.fill_between([len(species_idxs)-0.5,len(species_idxs)+1],[1,1],[1e05, 1e05],color='0.9')
polymorphism_axis.fill_between([-2,-0.5],[1,1],[1e05, 1e05],color='0.9')

#fst_axis.set_xlim([xmin,xmax])
#fst_axis.set_xticks([])
#fst_axis.set_ylim([1e-01,1e05])

fixation_axis.set_ylim([7e-01,1e05])

fixation_axis.fill_between([len(species_idxs)-0.5,len(species_idxs)+1],[3e-01,3e-01],[3e04, 3e04],color='0.9')
fixation_axis.fill_between([-2, -0.5],[3e-01,3e-01],[3e04, 3e04],color='0.9')
retention_axis.set_ylim([-0.1,1.1])

retention_axis.fill_between([len(species_idxs)-0.5,len(species_idxs)+1],[-0.1,-0.1],[1.1, 1.1],color='0.9')
retention_axis.fill_between([-2, -0.5],[-0.1,-0.1],[1.1, 1.1],color='0.9')

polymorphism_axis.set_xlim([xmin,xmax])
polymorphism_axis.set_xticks([])
polymorphism_axis.semilogy([0],[1e-08],'k.')

fixation_axis.set_xlim([xmin,xmax])
fixation_axis.set_xticks([])
fixation_axis.semilogy([0],[1e-08],'k.')

initial_axis.set_ylim([-0.04,0.24])
initial_axis.fill_between([len(species_idxs)-0.5,len(species_idxs)+1],[-0.05,-0.05],[0.5, 0.5],color='0.9')
initial_axis.fill_between([-2,-0.5],[-0.05,-0.05],[0.5, 0.5],color='0.9')
initial_axis.set_xlim([xmin,xmax])
initial_axis.set_xticks([])

#abundance_axis.set_xticks(xticks)
#abundance_axis.set_xticklabels(xticklabels,rotation=90)
#abundance_axis.get_xaxis().tick_top()
#abundance_axis.xaxis.set_label_position('top') 
#abundance_axis.tick_params(axis='x', labelsize=5,direction='out',length=3,pad=1)
#retention_axis.set_xlim([xmin,xmax])
#retention_axis.set_xticks([])

retention_axis.set_xlim([xmin,xmax])
retention_axis.set_xticks(xticks)
retention_axis.set_xticklabels(xticklabels,rotation=90)
retention_axis.tick_params(axis='x', labelsize=5,direction='out',length=3,pad=1)
retention_axis.get_xaxis().tick_bottom()


fig = pylab.gcf()
fig.savefig('%s/cross_species_summary.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
#change_fig.savefig('%s/species_freq_change.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
#focal_fig.savefig('%s/species_focal_freq_timecourse.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')

print ""
print "FINAL SUMMARY"
print NUM_SPECIES_GENETIC_CHANGES_CHECKED, "species examined for genetic changes"
print NUM_SPECIES_GENETIC_CHANGES_ABX, "species with ABX genetic change"
print NUM_SPECIES_GENETIC_CHANGES_ABX*1.0/NUM_SPECIES_GENETIC_CHANGES_CHECKED

garudgoodetal_checked = 801
garudgoodetal_changes = 96

print "Performing a Fisher's exact test w/ 801 resident populations and 96 replacements + modifications (%g percent)" % (garudgoodetal_changes*1.0/garudgoodetal_checked)


from scipy.stats import fisher_exact

oddsratio, pvalue = fisher_exact([[NUM_SPECIES_GENETIC_CHANGES_CHECKED-NUM_SPECIES_GENETIC_CHANGES_ABX, garudgoodetal_checked-garudgoodetal_changes],[NUM_SPECIES_GENETIC_CHANGES_ABX, garudgoodetal_changes]])
print "Fisher exact pvalue = %g" % pvalue
    
