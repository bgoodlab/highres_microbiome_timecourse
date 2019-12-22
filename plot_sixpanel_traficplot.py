###############################
#
# Import tools
#
################################
import matplotlib  
matplotlib.use('Agg') 
import pylab
import numpy
import sys
from math import log10, floor, ceil,log, fabs
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from numpy.random import binomial, random_sample, shuffle,random
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
#import calculate_within_species_fixations
import calculate_barcode_within_species_fixations as calculate_within_species_fixations
import figure_utils
import config
import diversity_utils
import sfs_utils
import stats_utils

################################################################################
#
# Standard header to read in argument information
#
################################################################################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('settingsfilename', type=str, help="settings-filename")
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
args = parser.parse_args()

debug = args.debug
chunk_size = args.chunk_size
settings_filename = args.settingsfilename
################################################################################

#vmin=0
#vmax=1
vmin=log(3e-03/(1-3e-03))
vmax=log((1-3e-03)/3e-03)
cmap='RdYlBu'

idx_letter_map = {0:'a', 1:'b',2:'c', 3:'d', 4:'e', 5:'f'}

from matplotlib.colors import LinearSegmentedColormap

basic_cols=['#d7191c', '#984ea3', '#2b83ba']
my_cmap=LinearSegmentedColormap.from_list('mycmap', basic_cols)

jet = cm = pylab.get_cmap(cmap)
jet = cm = my_cmap 
cmap = 'mycmap'
cNorm  = colors.Normalize(vmin=vmin, vmax=vmax)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

def get_prevalence_color(prevalence):
    clipped_prevalence = max([prevalence,3e-03])
    clipped_prevalence = min([clipped_prevalence,1-3e-03])
    colorVal = scalarMap.to_rgba(log(clipped_prevalence/(1-clipped_prevalence)))
    return colorVal    
            

variant_type_colors = {'4D': '#b3de69', '2D': '#b15928', '3D': '#b15928', '1D': '#ff7f00'}
marker_colors = {'preserved': '0.7', 'disrupted': basic_cols[-1]}
       
#desired_samples = parse_timecourse_data.morteza_samples
desired_samples = parse_timecourse_data.all_samples
allowed_variant_types=set(['1D','2D','3D','4D'])  
# Mininum coverage for frequency estimation vs interpolation 
min_coverage = 6
negative_min_coverage = 6
COLORED_LINEWIDTH=1
barcodes_instead_of_reads = True
private_prevalence_threshold = 1e-06 # Not quite private, but close
min_time = -10000
plot_genes = False
   
mpl.rcParams['font.size'] = 6.0
mpl.rcParams['lines.linewidth'] = 0.25
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'
mpl.rcParams['axes.linewidth'] = 0.5

# load settings
settings_file = open(settings_filename,"r")
settings_string = "\n".join(settings_file.readlines())
settings_file.close()
exec settings_string    

if barcodes_instead_of_reads:
    max_d = config.cluster_distance_threshold_barcodes
    snps_directory=parse_midas_data.barcode_directory+"/barcode_snps/"
else:
    snps_directory=parse_midas_data.default_snps_directory
    max_d = config.cluster_distance_threshold_reads

    

desired_species = species_names

num_species = len(desired_species)

fig_width = 8
fig_height = 1.7*3

pylab.figure(2,figsize=(5, 3))
stats_fig = pylab.gcf()
outer_grid = gridspec.GridSpec(2,1,height_ratios=[0.03,0.9],hspace=0.45) 
stats_grid  = gridspec.GridSpecFromSubplotSpec(2,3,height_ratios=[1,1],width_ratios=[1,1,1],hspace=0.5, wspace=0.2, subplot_spec=outer_grid[1])
colorbar_grid  = gridspec.GridSpecFromSubplotSpec(1,3,width_ratios=[0.35, 0.3,0.25],wspace=0.1, subplot_spec=outer_grid[0])

colorbar_axis = plt.Subplot(stats_fig, colorbar_grid[1])
stats_fig.add_subplot(colorbar_axis)

dummy_axis = plt.Subplot(stats_fig, colorbar_grid[0])
stats_fig.add_subplot(dummy_axis)
dummy_axis.set_ylim([0,1])
dummy_axis.set_xlim([0,1])

dummy_axis.spines['top'].set_visible(False)
dummy_axis.spines['right'].set_visible(False)
dummy_axis.spines['left'].set_visible(False)
dummy_axis.spines['bottom'].set_visible(False)

dummy_axis.set_xticks([])
dummy_axis.set_yticks([])

dummy_axis.plot([-1],[-1],'s',markeredgewidth=0,color='w',label='Variant type:')
dummy_axis.plot([-1],[-1],'s',markeredgewidth=0,color=variant_type_colors['2D'],label='(2D & 3D)')
dummy_axis.plot([-1],[-1],'s',markeredgewidth=0,color=variant_type_colors['4D'],label='syn (4D)')
dummy_axis.plot([-1],[-1],'s',markeredgewidth=0,color=variant_type_colors['1D'],label='non (1D)')
dummy_axis.legend(loc='center left',frameon=False,fontsize=6,numpoints=1,ncol=2,handlelength=0.5)   

dummy_axis = plt.Subplot(stats_fig, colorbar_grid[2])
stats_fig.add_subplot(dummy_axis)
dummy_axis.set_ylim([0,1])
dummy_axis.set_xlim([0,1])

dummy_axis.spines['top'].set_visible(False)
dummy_axis.spines['right'].set_visible(False)
dummy_axis.spines['left'].set_visible(False)
dummy_axis.spines['bottom'].set_visible(False)

dummy_axis.set_xticks([])
dummy_axis.set_yticks([])

dummy_axis.plot([-1],[-1],'s',markeredgewidth=0,color=marker_colors['preserved'],label='preserved')
dummy_axis.plot([-1],[-1],'s',markeredgewidth=0,color=marker_colors['disrupted'],label='disrupted')
dummy_axis.legend(loc='center left',frameon=False,fontsize=6,numpoints=1,ncol=1,handlelength=1,title='Private marker status')   

pylab.figure(1,figsize=(fig_width, fig_height))
fig = pylab.gcf()
outer_grid  = gridspec.GridSpec(3,2,height_ratios=[1,1,1],width_ratios=[1,1], hspace=0.25, wspace=0.05)

# First 

species_coverage_matrix, samples, species = parse_midas_data.parse_global_marker_gene_coverages()
sample_time_map = parse_timecourse_data.parse_sample_time_map()
species_ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, samples,min_time=min_time)
samples = numpy.array(samples)[sample_idxs]
species_coverage_matrix = species_coverage_matrix[:,sample_idxs]
total_coverage = species_coverage_matrix.sum(axis=0)
species_freq_matrix = numpy.clip(species_coverage_matrix*1.0/total_coverage,0, 2)    

mint = species_ts[0]-2
maxt = species_ts[-1]+2

freq_xticks = [species_ts[0]] + [20*i for i in xrange(0,8)]
freq_xticklabels = ['-2yr'] + [str(20*i) for i in xrange(0,8)]

max_abundance = 0
min_abundance = 1

species_freq_map = {species[species_idx]: species_freq_matrix[species_idx,:] for species_idx in xrange(0,len(species))}

    
fig_idx = 1
freq_axis = None

output_items = []
 
epoch_sets = [parse_timecourse_data.initial_plus_previous_epoch_intervals]

previous_samples = parse_timecourse_data.epoch_samples['previous']
initial_samples = parse_timecourse_data.epoch_samples['initial'] #+parse_timecourse_data.epoch_samples['disease']
middle_samples = parse_timecourse_data.epoch_samples['antibiotic']+parse_timecourse_data.epoch_samples['postantibiotic']+parse_timecourse_data.epoch_samples['disease']
final_samples = parse_timecourse_data.epoch_samples['prefinal']+parse_timecourse_data.epoch_samples['final']

previous_ts = [sample_time_map[sample] for sample in previous_samples]
initial_ts = [sample_time_map[sample] for sample in initial_samples] 
middle_ts = [sample_time_map[sample] for sample in middle_samples]
final_ts = [sample_time_map[sample] for sample in final_samples]

##
#
# Helper function for flipping alt allele
#
##
def flip_allele(allele):
    if allele=='A':
        return 'R'
    else:
        return 'A'

from math import fabs
#####
#
# This function tells us how we should polarize SNV clusters
#
####
def should_polarize(ts,fs,lower=0.3,upper=0.7,df=0.6):
    
     
    for tf in middle_ts:
        if (ts==tf).sum()<0.5:
            continue
                
        ff = fs[ts==tf][0]
                
        for t0 in initial_ts: 
            if (ts==t0).sum()<0.5:
                continue
                    
            f0 = fs[ts==t0][0]
                
            if (ff>upper and f0<lower)*(fabs(ff-f0)>df):
                return False
            elif (ff<lower and f0>upper)*(fabs(ff-f0)>df):
                return True
            else:
                pass
        
    for tf in reversed(final_ts):
        if (ts==tf).sum()<0.5:
            continue
                
        ff = fs[ts==tf][0]
                
        for t0 in initial_ts: 
            if (ts==t0).sum()<0.5:
                continue
                    
            f0 = fs[ts==t0][0]
            
            if (ff>upper and f0<lower)*(fabs(ff-f0)>df):
                return False
            elif (ff<lower and f0>upper)*(fabs(ff-f0)>df):
                return True
            else:
                pass
        
    
    # nothing happened so far, 
    # polarize based on initial timepoint
    if fs[0]>0.5:
        return True
    else:
        return False

for species_idx in xrange(0,len(desired_species)):        

    print species_idx, desired_species[species_idx]

    species_name = desired_species[species_idx]
    
    sys.stderr.write("Processing %s...\n" % species_name)

    
    sample_coverage_map = parse_midas_data.parse_median_coverage_map(species_name,snps_directory=snps_directory)
    
    desired_samples = parse_timecourse_data.all_samples
    
    dummy_samples, non_sfs_map = parse_midas_data.parse_within_sample_sfs(species_name, allowed_variant_types=set(['1D']))
    
    dummy_samples, syn_sfs_map = parse_midas_data.parse_within_sample_sfs(species_name, allowed_variant_types=set(['4D']))
    
    non_opportunity_vector = []
    syn_opportunity_vector = []
    for sample in desired_samples:
        if sample not in dummy_samples:
            continue
        
        within_sites, between_sites, total_sites = sfs_utils.calculate_polymorphism_rates_from_sfs_map(non_sfs_map[sample])
        
        non_opportunity_vector.append(total_sites)

        within_sites, between_sites, total_sites = sfs_utils.calculate_polymorphism_rates_from_sfs_map(syn_sfs_map[sample])
        
        syn_opportunity_vector.append(total_sites)            
    
    non_opportunity_vector = numpy.array(non_opportunity_vector)
    syn_opportunity_vector = numpy.array(syn_opportunity_vector)
    
    expected_nonsyn_ratio = numpy.median(non_opportunity_vector*1.0/syn_opportunity_vector)
    
    total_opportunities = numpy.median(non_opportunity_vector+syn_opportunity_vector)
    non_opportunities = expected_nonsyn_ratio/(1.0+expected_nonsyn_ratio)*total_opportunities
    syn_opportunities = 1.0/(1.0+expected_nonsyn_ratio)*total_opportunities
    
    
    sys.stderr.write("Loading core genes...\n")
    core_genes = core_gene_utils.parse_core_genes(species_name)
    personal_core_genes = core_gene_utils.parse_personal_core_genes(species_name)
    #non_shared_genes = core_gene_utils.parse_non_shared_reference_genes(species_name)
    non_shared_genes = personal_core_genes
    shared_pangenome_genes = core_gene_utils.parse_shared_genes(species_name)
    sys.stderr.write("Done! %d core genes and %d shared genes and %d non-shared genes\n" % (len(core_genes), len(shared_pangenome_genes), len(non_shared_genes)))
    
    sys.stderr.write("Loading SNV prevalences...\n")
    preexisting_snps = calculate_snp_prevalences.parse_population_freqs(species_name,polarize_by_consensus=True)
    # Polarization is such that it should match the output of parse_midas_data.parse_snps
    sys.stderr.write("Done! Loaded database of %d snps.\n" % len(preexisting_snps))
  
    sys.stderr.write("Loading fixations...\n")
    fixations_by_epoch = []
    gene_fixations_by_epoch = []
    for epoch_idx in xrange(0,len(epoch_sets)):
        
        epoch_set = epoch_sets[epoch_idx]
        
        snp_changes, gene_changes = calculate_within_species_fixations.load_within_species_fixations(species_name,allowed_epochs=epoch_sets[epoch_idx])
        #print epoch_set, len(snp_changes), len(gene_changes)
        
        
        fixations = []
        for snp_change in snp_changes:
            fixations.append((snp_change[0], snp_change[1]))
        fixations = set(fixations)
        fixations_by_epoch.append(fixations)
        
        gene_fixations = []
        for gene_name in gene_changes:
            gene_fixations.append(gene_name)
        gene_fixations = set(gene_fixations)
        gene_fixations_by_epoch.append(gene_fixations)
        
    sys.stderr.write("Done!")
        
    # Load gene coverage information for species_name
    sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
    gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name,allowed_samples=desired_samples, disallowed_genes=shared_pangenome_genes)
    sys.stderr.write("Done!\n")
    if len(gene_names)==0:
        sys.stderr.write("No genes, continuing!\n")
        continue
    
    marker_coverage_times, marker_coverage_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, gene_samples,min_time=min_time)

    gene_copynum_matrix = gene_depth_matrix*1.0/(marker_coverages+(marker_coverages==0))
    marker_coverages = marker_coverages[marker_coverage_idxs]
    gene_copynum_matrix = gene_copynum_matrix[:,marker_coverage_idxs]                

    times = []
    alt_matrix = []
    depth_matrix = []
    snp_infos = []
    
    opportunity_vectors = {}
    
    final_line_number = 0
    while final_line_number >= 0:
    
        sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
        samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_variant_types=set(['1D','2D','3D','4D']),chunk_size=chunk_size,allowed_samples=desired_samples, initial_line_number=final_line_number, allowed_genes=non_shared_genes, snps_directory=snps_directory)
        sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
    
        if len(opportunity_vectors)==0:
            for variant_type in allowed_variant_types:
                opportunity_vectors[variant_type] = numpy.zeros(len(samples)) 
    
            
        sample_ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, samples,min_time=min_time)

        Dbars = []
        for sample in samples:
            if sample in sample_coverage_map:
                Dbar = sample_coverage_map[sample]
            else:
                Dbar = 0
            Dbars.append(Dbar)
            
        Dbars = numpy.array(Dbars)
            

        # Calculate fixation matrix
        sys.stderr.write("Calculating allele freqs...\n")
        chunk_alts, chunk_depths, chunk_snp_infos = timecourse_utils.calculate_read_count_matrix(allele_counts_map, passed_sites_map, allowed_variant_types=set(['1D','2D','3D','4D']))    
        sys.stderr.write("Done!\n")
    
    
        sys.stderr.write("Calculating opportunities...\n")
        chunk_opportunity_vectors = timecourse_utils.calculate_opportunity_vectors(allele_counts_map, passed_sites_map, allowed_variant_types=set(['1D','2D','3D','4D']))    
        sys.stderr.write("Done!\n")
    
        for variant_type in chunk_opportunity_vectors:
            opportunity_vectors[variant_type] += chunk_opportunity_vectors[variant_type]
    
        chunk_alts = chunk_alts[:,sample_idxs]
        chunk_depths = chunk_depths[:,sample_idxs]
        
        #desired_sites = numpy.logical_not(numpy.logical_or((chunk_alts<=(0.1*chunk_depths)).all(axis=1), (chunk_alts>=(0.9*chunk_depths)).all(axis=1)))*((chunk_depths>0).sum(axis=1)>2)
        desired_sites = numpy.logical_not((chunk_alts<=(0.1*chunk_depths)).all(axis=1))*((chunk_depths>0).sum(axis=1)>2)
            
        chunk_alts = chunk_alts[desired_sites,:]
        chunk_depths = chunk_depths[desired_sites,:]
        chunk_allele_freqs = chunk_alts*1.0/(chunk_depths+(chunk_depths==0))
    
        if len(times)==0:
            times = sample_ts
        
                      
        if desired_sites.sum()>0:
            alt_matrix.append(chunk_alts)
            depth_matrix.append(chunk_depths)
            desired_site_idxs = numpy.nonzero(desired_sites)[0]
            for idx in desired_site_idxs:
                snp_infos.append(chunk_snp_infos[idx])  
        
    if len(alt_matrix)>0:     
        alt_matrix = numpy.vstack(alt_matrix)
        depth_matrix = numpy.vstack(depth_matrix) 
    else:
        alt_matrix = numpy.array([])
        depth_matrix = numpy.array([])
    
    
    sys.stderr.write("Done loading SNPs!\n")
    

    # set up figure axis
    species_row_idx = (species_idx / 3)
    species_col_idx = (species_idx % 3)
    
    stats_species_outer_grid = gridspec.GridSpecFromSubplotSpec(1,2,width_ratios=[0.3,1],wspace=0.1, subplot_spec=stats_grid[species_row_idx, species_col_idx])

    dnds_axis = plt.Subplot(stats_fig, stats_species_outer_grid[0])
    stats_fig.add_subplot(dnds_axis)
    dnds_axis.spines['top'].set_visible(False)
    dnds_axis.spines['right'].set_visible(False)
    dnds_axis.spines['left'].set_visible(False)
    dnds_axis.get_xaxis().tick_bottom()
    
    dnds_axis.set_ylim([0,1])
    dnds_axis.set_yticks([])
    dnds_axis.set_xlim([-0.55,1.55])
    dnds_axis.set_xticks([0.5])
    dnds_axis.set_xticklabels(['n=%d' % 0]) #,rotation='vertical')        

    marker_axis = plt.Subplot(stats_fig, stats_species_outer_grid[1])
    stats_fig.add_subplot(marker_axis)
    marker_axis.set_aspect('equal')
    marker_axis.spines['top'].set_visible(False)
    marker_axis.spines['right'].set_visible(False)
    marker_axis.spines['left'].set_visible(False)
    marker_axis.spines['bottom'].set_visible(False)
    marker_axis.set_yticks([])
    marker_axis.set_xticks([0])
    marker_axis.set_xticklabels(['m=%d' % 0])        
    marker_axis.get_xaxis().tick_bottom()
    
    species_row_idx = (species_idx%3)
    species_col_idx = (species_idx/3)
    species_grid  = gridspec.GridSpecFromSubplotSpec(2,1,height_ratios=[0.5,1],hspace=0.05, subplot_spec=outer_grid[species_row_idx, species_col_idx])
        
    abundance_axis = plt.Subplot(fig, species_grid[0])
    fig.add_subplot(abundance_axis)
    
    title_text = ('(%s) %s' % (idx_letter_map[species_idx], figure_utils.get_pretty_species_name(species_name)))
    abundance_axis.set_title(title_text,loc='left',fontsize=6)
    dnds_axis.text(-0.55,1.05,title_text)
    
    abundance_axis.set_xlim([mint,maxt])   
   
    #abundance_axis.spines['top'].set_visible(False)
    #abundance_axis.spines['right'].set_visible(False)
    abundance_axis.get_xaxis().tick_bottom()
    abundance_axis.get_yaxis().tick_left()
    
    abundance_axis.set_xticks([])
    abundance_axis.set_xticklabels([])
    
    species_freqs = species_freq_map[species_name]
    
    #max_abundance = max([max_abundance, species_freq_matrix[species_idx,:].max()])
    #min_abundance = min([min_abundance, species_freq_matrix[species_idx,:].min()])
    
    abundance_axis.fill_between([parse_timecourse_data.antibiotic_start, parse_timecourse_data.antibiotic_end],[1e-04,1e-04],[1,1],color='#CB0317',linewidth=0,alpha=0.5)
     
    abundance_axis.fill_between(species_ts, numpy.ones_like(species_freqs)*1e-05, species_freqs,color='#A6AAA9')
    abundance_axis.semilogy(species_ts, species_freqs,'k.-',markersize=2, linewidth=0.25)
    

    max_freq = species_freqs.max()
    power = floor(log10(max_freq))

    boosted_freq = ceil(max_freq/(10**(power)))
    
    if boosted_freq > 9:
        power+=1
        boosted_freq = 1
    
    max_abundance = boosted_freq*(10**power) #10**(ceil(log10(species_freqs.max())))
    min_abundance = 0 #max([10**(floor(log10(species_freqs.min()))), 3e-04])
    
    min_abundance = 3e-04
    max_abundance = 1
    
    
    abundance_axis.set_ylim([min_abundance,max_abundance])
    #abundance_axis.set_yticks([max_abundance])
    freq_axis = plt.Subplot(fig, species_grid[1])
    fig.add_subplot(freq_axis)
    
    #freq_axis.spines['top'].set_visible(False)
    #freq_axis.spines['right'].set_visible(False)
    freq_axis.get_xaxis().tick_bottom()
    freq_axis.get_yaxis().tick_left()
    
    freq_axis.set_xlim([mint,maxt])
    freq_axis.set_xticks([])
    freq_axis.set_xticklabels([])
       
    freq_axis.set_ylim([0,1.1])
    
    freq_axis.fill_between([parse_timecourse_data.antibiotic_start, parse_timecourse_data.antibiotic_end],[0,0],[1.02,1.02],color='#CB0317',linewidth=0,alpha=0.5)
    
    
    if species_col_idx==0:
        freq_axis.set_ylabel('Allele\nfrequency')
        abundance_axis.set_ylabel('Relative\nabundance')    
    else:
        freq_axis.set_yticklabels([])
        abundance_axis.set_yticklabels([])
    
    if species_row_idx == (num_species/2)-1:
        freq_axis.set_xlabel('Time (days)')
        freq_axis.set_xticks(freq_xticks)
        freq_axis.set_xticklabels(freq_xticklabels)
        freq_axis.tick_params(axis='x', labelsize=6,direction='out',length=2,pad=1)
        
    
    abundance_axis.tick_params(axis='y', labelsize=6,direction='out',length=2,pad=1)
    abundance_axis.tick_params(axis='y', which='minor',labelsize=6,direction='out',length=1,pad=1)
    
    freq_axis.tick_params(axis='y', labelsize=6,direction='out',length=2,pad=1)
        
    
        
    num_colored_mutations = 0
    num_total_mutations = 0
    num_privates = 0
    num_high_freqs = 0
    num_markers = 0
    num_preserved_private_markers = 0
    num_disrupted_private_markers = 0
    
    # Calculate min coverage thresholds
    min_coverages = numpy.ones_like(times)*min_coverage
    min_coverages[times<0] = negative_min_coverage
            
    
    if len(alt_matrix)>0:
        
        cluster_As = []
        cluster_Ds = []
        cluster_prevalences = []
        cluster_infos = []
        variant_type_counts = {variant_type:0 for variant_type in allowed_variant_types}
        polarized_prevalences = []
        
        marker_As = []
        marker_Ds = []
        marker_infos = []
        
        high_freq_infos = []
        
        background_As = []
        background_Ds = []
        
        for mutation_idx in xrange(0,len(snp_infos)):
        
            num_total_mutations += 1
            
            chromosome, location, gene_name, variant_type, alt_allele = snp_infos[mutation_idx]
            alts = alt_matrix[mutation_idx,:]
            depths = depth_matrix[mutation_idx,:]
            good_idxs = (depths>=min_coverages)
            
            copynums = depths/Dbars
            
            num_highcoverage_samples = (Dbars>=min_coverages).sum()
            
            num_good_timepoints = (good_idxs).sum()
            if num_good_timepoints<5 or (num_good_timepoints<0.75*num_highcoverage_samples):
                continue
                    
            #alts = alts*(depths>=2)
            #depths = depths*(depths>=2)
            freqs = alts*1.0/(depths+(depths==0))
        
            masked_times = times[good_idxs]
            masked_freqs = freqs[good_idxs]
            masked_depths = depths[good_idxs]
        
            # Figure out whether it is a private SNV or not  
            private_snp = False
            if len(preexisting_snps)>0:
                if (chromosome,location) not in preexisting_snps:
                    private_snp = True
                else:
                    if preexisting_snps[(chromosome,location)] < private_prevalence_threshold:
                        private_snp = True
            # Since we polarized by consensus, all prevalences should be <= 0.5
            if private_snp:
                prevalence = -1e-06
            else:
                prevalence = preexisting_snps[(chromosome,location)] 
            
            always_high_freq = (masked_freqs>=0.5).all() and ((masked_freqs>=0.8).sum()*1.0/len(masked_freqs) >= 0.8)
            
            # One of the preserved markers           
            if private_snp and always_high_freq and (gene_name in core_genes):
                num_preserved_private_markers += 1
            
            high_freq_marker_snv = always_high_freq and (gene_name in core_genes)
            
            #if high_freq_marker_snv and (prevalence<0.1):
                # add to separate list
            #    high_freq_infos.append((chromosome, location, alt_allele, variant_type, gene_name, prevalence))
                
            # If it never makes it to intermediate frequencies, don't bother plotting it    
            if (masked_freqs<0.2).all() or (masked_freqs>0.8).all():
                continue
            
            # Then decide whether to color it:
            id = (gene_name, chromosome, location, variant_type)
            
            color_condition = False
            for epoch_idx in xrange(0,len(epoch_sets)):    
                if (chromosome, location) in fixations_by_epoch[epoch_idx]:
                    color_condition = True
                    color = 'r'
                    break   
            
               
            color_condition = color_condition or color_rule_in_condition(species_idx, chromosome, location, gene_name, variant_type, masked_times, masked_freqs, masked_depths)
            
            color_condition = color_condition and (not color_rule_out_condition(species_idx, chromosome, location, gene_name, variant_type, masked_times, masked_freqs, masked_depths))
            
            if color_condition:
            
                # One of the colored ones!
                num_colored_mutations+=1
                
                # get color
                
                cluster_As.append(alts)
                cluster_Ds.append(depths)
                cluster_infos.append((snp_infos[mutation_idx][0], snp_infos[mutation_idx][1], variant_type, gene_name, prevalence, alt_allele))
                
                cluster_prevalences.append(prevalence)
                
                variant_type_counts[variant_type] += 1
                
                
            else:
                # We don't color it 
              
                background_As.append(alts)
                background_Ds.append(depths)
                
                
        cluster_As = numpy.array(cluster_As)
        cluster_Ds = numpy.array(cluster_Ds)
        
        
        background_As = numpy.array(background_As)
        background_Ds = numpy.array(background_Ds)
        # Limit # of background SNVs to ~300
        max_num_backgrounds = 300
        if len(background_As)>max_num_backgrounds:
            background_idxs = numpy.arange(0,len(background_As))
            shuffle(background_idxs)
            background_idxs = background_idxs[0:max_num_backgrounds]
            background_As = background_As[background_idxs,:]
            background_Ds = background_Ds[background_idxs,:]
        
            
        output_items.append('# %s' % species_name)
        snp_output_items = []
        
        start_idx = 0    
             
        if len(cluster_As)==0:
            
            sys.stderr.write("No temporally variable SNVs, just plotting background...\n")
            # No temporally variable SNVs. Just print background
        
            background_cluster_map = cluster_utils.cluster_snps_by_distance(background_As, background_Ds,max_d=max_d)
            
            for cluster_id in background_cluster_map:
                
                avg_fs, total_Ds = background_cluster_map[cluster_id]['centroid']
            
                good_idxs = (total_Ds>=min_coverages)
            
                if good_idxs.sum()<1.5:
                    continue
            
                masked_times = times[good_idxs]
                masked_avg_fs = avg_fs[good_idxs]
                polarize = should_polarize(masked_times, masked_avg_fs)
                if polarize:
                    masked_avg_fs = 1-masked_avg_fs
            
                for snp_idx, flip in background_cluster_map[cluster_id]['snps']:
                        
                    alts = background_As[snp_idx,:]
                    depths = background_Ds[snp_idx,:]
                    freqs = alts/(depths+(depths==0))
                    good_idxs = (depths>=min_coverages)
                
                    masked_times = times[good_idxs]
                    masked_freqs = freqs[good_idxs]
                
                    if flip:
                        masked_freqs = 1-masked_freqs
                
                    # Now get it in line with what the SNV cluster is now. 
                    if polarize:
                        masked_freqs = 1-masked_freqs
                
                        
                    freq_axis.plot(masked_times, masked_freqs, '-', color='0.7', alpha=0.5, markersize=3,label=gene_name,linewidth=0.25,zorder=1)
     
            sys.stderr.write("Done!\n")
                           
        if len(cluster_As)>0:
        
            sys.stderr.write("Proceeding with %d temporally variable SNVs and %d background SNVs...\n" % (len(cluster_As), len(background_As)))
            
            
            cluster_map = cluster_utils.fast_cluster_snps_by_distance(cluster_As[:,start_idx:],cluster_Ds[:,start_idx:],max_d=max_d,min_coverage=min_coverage)
            if len(background_As)>0:
                background_cluster_map = cluster_utils.cluster_secondary_snps_by_distance(cluster_map, background_As[:,start_idx:], background_Ds[:,start_idx:],max_d=max_d,min_coverage=min_coverage)
            else:
                background_cluster_map = {}
            sys.stderr.write("Proceeding with %d temporally variable clusters, %d total clusters\n" % (len(cluster_map), len(set(cluster_map.keys())|set(background_cluster_map.keys()))))
            
            
            cluster_ids = cluster_map.keys()
            cluster_sizes = [len(cluster_map[cluster_id]['snps']) for cluster_id in cluster_ids]
            sorted_cluster_sizes, sorted_cluster_ids = zip(*sorted(zip(cluster_sizes, cluster_ids), reverse=True))
            total_snvs = sum(sorted_cluster_sizes)
            
            # Get some metadata for each main cluster
            cluster_polarization_map = {}
            cluster_color_map = {}
            cluster_size_map = {}
            cluster_fraction_map = {}
            for cluster_id in sorted_cluster_ids:
                
                color = None
                avg_fs, total_Ds = cluster_map[cluster_id]['centroid']
            
                good_idxs = (total_Ds>=min_coverages[start_idx:])
                if good_idxs.sum()<1.5:
                    continue
            
                # TODO: left here!
            
                masked_times = times[start_idx:][good_idxs]
                masked_avg_fs = avg_fs[good_idxs]
 
                polarize = should_polarize(masked_times, masked_avg_fs)
                cluster_polarization_map[cluster_id] = polarize
                if polarize:
                    masked_avg_fs = 1-masked_avg_fs
            
                # add offset to make sure it isn't plotted, but still assign a color
                line, = freq_axis.plot(masked_times, masked_avg_fs-5, '-o', alpha=1, markersize=3, markeredgecolor='none', zorder=5, linewidth=1)
                cluster_color = pylab.getp(line,'color')
                
                cluster_color_map[cluster_id] = cluster_color
                cluster_size_map[cluster_id] = len(cluster_map[cluster_id]['snps'])
                cluster_fraction_map[cluster_id] = cluster_size_map[cluster_id]*1.0/total_snvs 
            
            
            # Plot SNVs in background cluster map
            for cluster_id in background_cluster_map:
                
                avg_fs, total_Ds = background_cluster_map[cluster_id]['centroid']
                good_idxs = (total_Ds>=min_coverages[start_idx:])
            
                if good_idxs.sum()<1.5:
                    continue
            
                masked_times = times[start_idx:][good_idxs]
                masked_avg_fs = avg_fs[good_idxs]
                polarize = should_polarize(masked_times, masked_avg_fs)
                if polarize:
                    masked_avg_fs = 1-masked_avg_fs
            
                for snp_idx, flip in background_cluster_map[cluster_id]['snps']:
                        
                    alts = background_As[snp_idx,:]
                    depths = background_Ds[snp_idx,:]
                
                    freqs = alts/(depths+(depths==0))
                    good_idxs = depths>=min_coverages
                    masked_times = times[good_idxs]
                    masked_freqs = freqs[good_idxs]
                
                    if flip:
                        masked_freqs = 1-masked_freqs
                
                    # Now get it in line with what the SNV cluster is now. 
                    if polarize:
                        masked_freqs = 1-masked_freqs
                   
                    freq_axis.plot(masked_times, masked_freqs, '-', color='0.7', alpha=0.5, markersize=3,label=gene_name,linewidth=0.25,zorder=1)
            
            output_items.append('# Clusters:')
                
            # Now plot temporally variable SNVs
            total_plotted = 0    
            for cluster_id in sorted_cluster_ids:
                
                # For estimating the cluster freqs
                total_As = numpy.zeros_like(times)
                total_Ds = numpy.zeros_like(times)
                
                # Recall data from before...
                polarize = cluster_polarization_map[cluster_id]
                cluster_color = cluster_color_map[cluster_id]
                cluster_size = cluster_size_map[cluster_id]
                cluster_fraction = cluster_fraction_map[cluster_id]
                
                
                if total_plotted>300 and cluster_fraction<0.05:
                    plot_cluster=False
                else:
                    plot_cluster=True
                total_plotted += cluster_size
                
                # If cluster is big, we won't plot the whole cluster, just the top 300 ones
                # We won't plot the whole cluster. Just the top 300 random ones
                p = 300.0/cluster_size
                
                
                # Calculate cluster trajectory
                for snp_idx, flip in cluster_map[cluster_id]['snps']:
                
                    alts = cluster_As[snp_idx,:]
                    depths = cluster_Ds[snp_idx,:]
                 
                    # First get it in line with what the SNV cluster *was*
                    if flip:
                        alts = depths-alts
                        
                    # Now get it in line with what the SNV cluster is now. 
                    if polarize:
                        alts = depths-alts
                        
                    total_As += alts
                    total_Ds += depths
        
                output_str = '# %s' % ("; ".join([cluster_color]+["%d,%d,%d" % (t,A,D) for t,A,D in zip(times,total_As,total_Ds)]))
                output_items.append(output_str)
            
                cluster_freqs = total_As*1.0/(total_Ds+(total_Ds==0))
                good_cluster_idxs = (total_Ds>=min_coverages)
        
                dN = 0
                dS = 0
                for snp_idx, flip in cluster_map[cluster_id]['snps']:
                    
                    alts = cluster_As[snp_idx,:]
                    depths = cluster_Ds[snp_idx,:]
                    
                    good_idxs = (depths>=min_coverages)
                    bad_idxs = numpy.logical_not(good_idxs)
                    bad_and_good_cluster_idxs = numpy.logical_and(bad_idxs,good_cluster_idxs)
                
                    if good_idxs.sum()<3:
                        continue
                
                    chromosome, location, variant_type, gene_name, prevalence, alt_allele = cluster_infos[snp_idx]
                    
                    if variant_type=='1D':
                        dN+=1
                    elif variant_type=='4D':
                        dS+=1
                    else:
                        pass
                     
                    # First get it in line with what the SNV cluster *was*
                    if flip:
                        alts = depths-alts
                        prevalence = 1-prevalence
                        alt_allele = flip_allele(alt_allele)
                
                    # Now get it in line with what the SNV cluster is now. 
                    if polarize:
                        alts = depths-alts
                        prevalence=1-prevalence
                        alt_allele = flip_allele(alt_allele)
                
                    if prevalence>=1:
                        # A private SNV that is disrupted!
                        num_disrupted_private_markers += 1
                
                    polarized_prevalences.append(prevalence)
                    
                    # Decide whether to plot
                    if random()>p:
                        continue
                    
                    if not plot_cluster:
                        continue
                    
                    
                        
                    # Fill in undersampled points using cluster trajectory
                    if False and ((bad_and_good_cluster_idxs.sum())>0):
                        alts[bad_and_good_cluster_idxs] = binomial(numpy.array(min_coverages[bad_and_good_cluster_idxs],dtype=numpy.int32),cluster_freqs[bad_and_good_cluster_idxs])
                        depths[bad_and_good_cluster_idxs] = min_coverages[bad_and_good_cluster_idxs]
                    
                    good_idxs = (depths>=min_coverages)
                    masked_times = times[good_idxs]
                    masked_As = alts[good_idxs]
                    masked_Ds = depths[good_idxs]
                    masked_freqs = masked_As*1.0/(masked_Ds)
                
                    variable = (((masked_freqs<=0.2).any() and (masked_freqs>=0.75).any()) or ((masked_freqs<=0.25).any() and (masked_freqs>=0.8).any()))
                    if not variable:
                        pass
                        #print times
                        #print alts
                        #print depths
                        #print alts*1.0/(depths+(depths==0))
                    
                    
                    output_str = "\t".join([species_name, "%s|%d" % (chromosome, location), variant_type, gene_name, str(prevalence), alt_allele, cluster_color])
                    snp_output_items.append( output_str )
        
                    #print chromosome, location, gene_name, masked_freqs
                      
                    colorVal = get_prevalence_color(prevalence)
                    
                    if prevalence>=1:
                        linestyle = 'v-'
                        zorder=5
                    elif prevalence<=0:
                        linestyle = '^-'
                        zorder=5
                    else:
                        linestyle = '-'
                        zorder=4
                    
                    if cluster_size>30:
                        linewidth=0.25
                    else:
                        linewidth=0.75
                    
                        
                        
                    freq_axis.plot(masked_times, masked_freqs, linestyle, color=cluster_color, alpha=0.5, markersize=1, markeredgecolor='none', zorder=zorder, linewidth=linewidth)
                
                # Do dN/dS calculate
                if dS>0.5:
                    print "dN/dS=%g " % (dN*1.0/dS/expected_nonsyn_ratio),
                else:
                    print "dN/dS=inf ",
                
                dtot = dN+dS
                
                if dtot>0.5:
                    pvalue = stats_utils.calculate_binomial_LRT_pvalue(dN,dtot,non_opportunities-dN,total_opportunities-dtot) 
                    print " (nN=%d, nS=%d, p=%g)" % (dN,dS, pvalue)    
                    print expected_nonsyn_ratio
                    print expected_nonsyn_ratio*1.0/(1+expected_nonsyn_ratio)    
                
                # Plot the cluster trajectory if cluster is big enough
                if plot_cluster and cluster_size>30:
                    
                    freq_axis.plot(times[good_cluster_idxs], cluster_freqs[good_cluster_idxs], 'k-', alpha=1, zorder=zorder+1, linewidth=0.75)
                    freq_axis.plot(times[good_cluster_idxs], cluster_freqs[good_cluster_idxs], 'o', color=cluster_color, alpha=1, markersize=2, zorder=zorder+2, linewidth=0.5)
    
           
        sys.stderr.write("%d markers preserved, %d disrupted\n" % (num_preserved_private_markers, num_disrupted_private_markers))
        output_items.append('# Private marker SNVs:')
        output_items.append("# Preserved=%d; Disrupted=%d" % (num_preserved_private_markers, num_disrupted_private_markers))

        output_items.append('# dN/dS:')
        for variant_type in sorted(opportunity_vectors):
            output_str = "# %s; Observed=%d; Expected=%s" % (variant_type, variant_type_counts[variant_type], ", ".join([str(n) for n in opportunity_vectors[variant_type]]))
            output_items.append(output_str)
        
        total_markers = num_disrupted_private_markers + num_preserved_private_markers

        output_items.extend(snp_output_items)
        
        if total_markers==0:
            total_markers+=1
        
        marker_axis.pie([num_disrupted_private_markers, num_preserved_private_markers],colors=[marker_colors['disrupted'],marker_colors['preserved']],wedgeprops = {'linewidth': 0})
    
        marker_axis.set_xticks([0])
        marker_axis.set_xticklabels(['m=%d' % total_markers])        
        marker_axis.tick_params(axis='x', labelsize=6,direction='out',length=0,pad=0)
        

        total_variant_type_counts = sum(variant_type_counts.values())
        if total_variant_type_counts>0.5:
            upper = 1.0
            for variant_type in ['2D','3D','1D','4D']:
                dnds_axis.bar([-0.5],[upper], width=0.8,color=variant_type_colors[variant_type],edgecolor='none')
                upper-= variant_type_counts[variant_type]*1.0/total_variant_type_counts
            
            dnds_axis.set_xticklabels(['n=%d' % total_variant_type_counts]) #,rotation='vertical')        
        
        polarized_prevalences = numpy.array(polarized_prevalences)
        total_prevalences = len(polarized_prevalences)
        if total_prevalences>0:
        
            fraction_private = (polarized_prevalences<0).sum()*1.0/total_prevalences
            fraction_public = (polarized_prevalences>1).sum()*1.0/total_prevalences
            fraction_intermediate = 1-fraction_private-fraction_public
            
            intermediate_prevalences = polarized_prevalences[(polarized_prevalences>=0)*(polarized_prevalences<=1)]
            
            if len(intermediate_prevalences)>100:
                shuffle(intermediate_prevalences)
                intermediate_prevalences = intermediate_prevalences[0:100]
            
            intermediate_prevalences = numpy.sort(intermediate_prevalences)
            
            intermediate_dy = fraction_intermediate/len(intermediate_prevalences)
            
            
            upper = 1.0
            if fraction_public>0:
                dnds_axis.bar([0.7], [fraction_public], width=0.8,bottom=[upper-fraction_public],color=get_prevalence_color(1.1),linewidth=0,alpha=0.5)
            upper-= fraction_public
            
            for prevalence_idx in reversed(xrange(0,len(intermediate_prevalences))):
                prevalence = intermediate_prevalences[prevalence_idx]
                    
                dnds_axis.bar([0.7], [intermediate_dy], width=0.8,bottom=[upper-intermediate_dy],color=get_prevalence_color(prevalence),linewidth=0)
                upper-= intermediate_dy
            
            if fraction_private>0:
                dnds_axis.bar([0.7], [upper], width=0.8,color=get_prevalence_color(-0.1),linewidth=0,alpha=0.5)
            
m = freq_axis.scatter([200],[1],c=[0.5], vmin=vmin, vmax=vmax, cmap=jet, marker='^')



cbar = fig.colorbar(m,cax=colorbar_axis,orientation='horizontal')
cbar.set_label('Prevalence of allele across hosts',labelpad=-25)        
#cbar.set_ticks([0,0.2,0.4,0.6,0.8,1])

cbar.set_ticks([log(1e-02/(1-1e-02)),log(1e-01/(1-1e-01)),0,log((1-1e-01)/1e-01),log((1-1e-02)/1e-02)])
cbar.set_ticklabels(['1%','10%','50%','90%','99%'])

sys.stderr.write("Saving final PNG image...\t")
if filename.endswith('png'):
    fig.savefig(filename, bbox_inches='tight', dpi=300, transparent=True)
    stats_fig.savefig(stats_filename, bbox_inches='tight', transparent=True)
else:
    fig.savefig(filename, bbox_inches='tight', transparent=True)
    stats_fig.savefig(stats_filename, bbox_inches='tight', transparent=True)
    
pylab.close(fig)
sys.stderr.write("Done!\n")

if len(output_filename)!="":
    sys.stderr.write("Saving output file!\n")
    output_file = open(output_filename,"w")
    output_file.write("\n".join(output_items))
    output_file.close()
    sys.stderr.write("Done!\n")
