import matplotlib  
matplotlib.use('Agg') 
import parse_midas_data
import parse_timecourse_data
import bacterial_phylogeny_utils
import pylab
import sys
import numpy
from math import log10, fabs, log
import figure_utils

import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint

mpl.rcParams['font.size'] = 8
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

min_marker_coverage = 20
min_prevalence=5

bacteroides_color = '#084594'
alistipes_color = '#B10026'
rest_color = '0.7'

bacteroides_colors = ['#fff5eb', '#fee6ce', '#fdd0a2', '#fdae6b', '#fd8d3c', '#f16913', '#d94801', '#8c2d04']
bacteroides_colors = ['#fee6ce', '#fdd0a2', '#fdae6b', '#fd8d3c', '#f16913', '#d94801', '#a63603', '#7f2704']
bacteroides_colors = bacteroides_colors[::-1] # reverse order
current_bacteroides_color = 0

alistipes_colors = ['#fbb4b9', '#f768a1', '#ae017e']
alistipes_colors = alistipes_colors[::-1] # reverse order
current_alistipes_color = 0

eubacterium_colors = ['#bae4b3', '#74c476', '#238b45'] #['#e5f5e0', '#a1d99b', '#31a354']
eubacterium_colors = eubacterium_colors[::-1] # reverse order
current_eubacterium_color = 0

lachno_colors = ['#bdd7e7', '#6baed6', '#2171b5'] #['#deebf7', '#9ecae1', '#3182bd']
lachno_colors = lachno_colors[::-1] # reverse order

current_lachno_color = 0

other_colors = ['#e6ab02', '#ef3b2c', '#6a51a3']# ]
current_other_color = 0


species_coverage_matrix, samples, species = parse_midas_data.parse_global_marker_gene_coverages()

print samples

sample_time_map = parse_timecourse_data.parse_sample_time_map()

ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, samples, min_time=-1000)

max_dt = 2.0


dts = numpy.diff(ts)*1.0

left_dts = numpy.hstack([[500.0],dts/4])
right_dts = numpy.hstack([dts/4,[500.0]])

dts = numpy.fmin(left_dts, right_dts)
dts = numpy.fmin(dts,max_dt/2)

print "ts:", ts

print "dts: ", dts

#ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, parse_timecourse_data.highcoverage_samples)

species_coverage_matrix = species_coverage_matrix[:,sample_idxs]

total_coverage = species_coverage_matrix.sum(axis=0)
species_freq_matrix = numpy.clip(species_coverage_matrix*1.0/total_coverage,0, 2)    

species_max_freqs = species_freq_matrix.max(axis=1)

print total_coverage

shannon_diversity = -1*(species_freq_matrix*numpy.log2(species_freq_matrix+(species_freq_matrix==0))).sum(axis=0)

js_divergence_matrix = []
# Calculate JSDs
for t0 in xrange(0,len(shannon_diversity)):
    js_divergences = []
    for t in xrange(0,len(shannon_diversity)):
    
        initial_freqs = (species_freq_matrix[:,t0])
        current_freqs = species_freq_matrix[:,t]
        mixture_freqs = initial_freqs*0.5+current_freqs*0.5
    
        good_idxs = (mixture_freqs>0)
    
        initial_freqs = initial_freqs[good_idxs]
        current_freqs = current_freqs[good_idxs]
        mixture_freqs = mixture_freqs[good_idxs]
    
        js_divergence = 0.5*(current_freqs*numpy.log2((current_freqs+(current_freqs==0))/mixture_freqs)).sum() + 0.5*(initial_freqs*numpy.log2((initial_freqs+(initial_freqs==0))/mixture_freqs)).sum()
        js_divergences.append(js_divergence)
    js_divergence_matrix.append(js_divergences)

js_divergence_matrix = numpy.array(js_divergence_matrix)

# Set up figure
fig = plt.figure(figsize=(14, 2.7))

# Set up grids to hold figure panels
outer_grid = gridspec.GridSpec(1, 2, width_ratios=[4,2], wspace=0.1)

trajectory_grid = gridspec.GridSpecFromSubplotSpec(2, 1, height_ratios=[0.3,1],
                subplot_spec=outer_grid[0], hspace=0.05) #, hspace=0.08)

#alpha_axis = plt.Subplot(fig, trajectory_grid[0])
#fig.add_subplot(alpha_axis)

#alpha_axis.set_ylabel('$\\alpha$ diversity \n (Shannon entropy)')
 
#alpha_axis.spines['top'].set_visible(False)
#alpha_axis.spines['right'].set_visible(False)
#alpha_axis.get_xaxis().tick_bottom()
#alpha_axis.get_yaxis().tick_left()

beta_axis = plt.Subplot(fig, trajectory_grid[0])
fig.add_subplot(beta_axis)

beta_axis.set_ylabel('$\\beta$ diversity \n (JS distance)')
 
beta_axis.spines['top'].set_visible(False)
beta_axis.spines['right'].set_visible(False)
beta_axis.get_xaxis().tick_bottom()
beta_axis.get_yaxis().tick_left()


freq_axis = plt.Subplot(fig, trajectory_grid[1])
fig.add_subplot(freq_axis)

freq_axis.set_ylabel('Species composition (%)')
freq_axis.set_xlabel('Day')
#freq_axis.set_ylim([1e-04,1])
freq_axis.set_ylim([0,1])
freq_axis.set_yticks([])
 
freq_axis.spines['top'].set_visible(False)
freq_axis.spines['right'].set_visible(False)
freq_axis.get_xaxis().tick_bottom()
freq_axis.get_yaxis().tick_left()


# Fill in markers
beta_axis.plot([parse_timecourse_data.antibiotic_start,parse_timecourse_data.antibiotic_end],[0.41,0.41],'k|-',markersize=4)

beta_axis.plot([parse_timecourse_data.lyme_infection],[0.30],'k|',markersize=4)

beta_axis.plot([parse_timecourse_data.antibiotic_start+2],[0.30],'k>',markersize=2.5)
beta_axis.plot([parse_timecourse_data.lyme_infection,parse_timecourse_data.antibiotic_start+2],[0.30,0.3],'k-')


beta_axis.plot([parse_timecourse_data.hrv_infection-1, parse_timecourse_data.hrv_infection+1], [0.30,0.30], 'k|-',markersize=4)

beta_axis.text(parse_timecourse_data.lyme_infection+(parse_timecourse_data.antibiotic_start-parse_timecourse_data.lyme_infection)*0.5,0.36,'Lyme \n disease',ha='center',fontsize=7)
beta_axis.text(parse_timecourse_data.antibiotic_start,0.47,'Antibiotics\n(doxycycline)',ha='left',fontsize=7)
beta_axis.text(parse_timecourse_data.hrv_infection,0.36,'HRV',ha='center',fontsize=7)



beta_axis.set_ylim([0,0.45])
beta_axis.set_yticks([0,0.1,0.2,0.3,0.4])

#alpha_axis.set_xlim([ts[0]-1,160])
beta_axis.set_xlim([ts[0]-1,160])
freq_axis.set_xlim([ts[0]-1,160])

xticks=[ts[0]]+[20*i for i in xrange(0,9)]
xticklabels=['-2yr']+[str(20*i) for i in xrange(0,9)]

#alpha_axis.set_xticks(xticks)
#alpha_axis.set_xticklabels([])

beta_axis.set_xticks(xticks)
beta_axis.set_xticklabels([])

freq_axis.set_xticks(xticks)
freq_axis.set_xticklabels(xticklabels)


legend_axis = plt.Subplot(fig, outer_grid[1])
fig.add_subplot(legend_axis)

legend_axis.set_ylim([0,1])
legend_axis.set_xlim([0,1])

legend_axis.spines['top'].set_visible(False)
legend_axis.spines['right'].set_visible(False)
legend_axis.spines['left'].set_visible(False)
legend_axis.spines['bottom'].set_visible(False)

legend_axis.set_xticks([])
legend_axis.set_yticks([])  


good_species_list = parse_midas_data.parse_good_species_list()

good_species_pretty_list = figure_utils.get_pretty_species_names(good_species_list)

pretty_species_name_map = {species_name: pretty_species_name for species_name, pretty_species_name in zip(good_species_list, good_species_pretty_list)}

freq_threshold = 5e-02

family_freqs = {}
genus_family_map = bacterial_phylogeny_utils.get_genus_family_map()

phylum_freqs = {}
genus_phylum_map = bacterial_phylogeny_utils.get_genus_phylum_map()


display_freqs = []
display_idxs = []

rest_freqs = numpy.zeros_like(species_freq_matrix[0,:])

for i in xrange(0,len(species)):
        
    species_coverages = species_coverage_matrix[i,:]
    species_freqs = species_freq_matrix[i,:]
    species_name = species[i]
    
    family_name = bacterial_phylogeny_utils.get_family_name(species_name, genus_family_map)
    if family_name not in family_freqs:
        family_freqs[family_name] = numpy.zeros_like(rest_freqs)        
    family_freqs[family_name] += species_freqs
    
    phylum_name = bacterial_phylogeny_utils.get_phylum_name(species_name, genus_phylum_map)
    if phylum_name not in phylum_freqs:
        phylum_freqs[phylum_name] = numpy.zeros_like(rest_freqs)        
    phylum_freqs[phylum_name] += species_freqs
    
    
    if (species_freqs>freq_threshold).any() or species_name.startswith('Alistipes_finegoldii'):
        # it's a display freq!
        display_freqs.append(species_freqs)
        display_idxs.append(i)
        print species_name, species_freqs.max(), species_freqs.mean()
    
        
          
    else:
        rest_freqs += species_freqs

# sort families in ascending order of frequency
sorted_families = list(family_freqs.keys())
sorted_families = list(sorted(sorted_families, key = lambda f: family_freqs[f][0]))
family_order_map = {}
for idx in xrange(0,len(sorted_families)):
    family_order_map[sorted_families[idx]] = idx
    if (family_freqs[sorted_families[idx]]>freq_threshold).any():
        print sorted_families[idx]

# sort families in ascending order of frequency
sorted_phyla = list(phylum_freqs.keys())
sorted_phyla = list(sorted(sorted_phyla, key = lambda f: phylum_freqs[f][0]))
phylum_freq_matrix = numpy.array([phylum_freqs[phylum] for phylum in sorted_phyla])

phylum_order_map = {}
for idx in xrange(0,len(sorted_phyla)):
    phylum_order_map[sorted_phyla[idx]] = idx
    if (phylum_freqs[sorted_phyla[idx]]>freq_threshold).any():
        print sorted_phyla[idx]

#print family_order_map

# sort species in descending order of family, then in descending order of abundance        
sorted_idxs = range(0,len(display_idxs))
#sorted_idxs = list(sorted(sorted_idxs, key = lambda idx: (family_order_map[bacterial_phylogeny_utils.get_family_name(species[idx])], (display_freqs[idx][8])/2.0)))
#sorted_idxs = list(sorted(sorted_idxs, key = lambda idx: (family_order_map[bacterial_phylogeny_utils.get_family_name(species[idx])], (display_freqs[idx][0])/2.0)))
sorted_idxs = list(sorted(sorted_idxs, key = lambda idx: (phylum_order_map[bacterial_phylogeny_utils.get_phylum_name(species[display_idxs[idx]])], family_order_map[bacterial_phylogeny_utils.get_family_name(species[display_idxs[idx]])], (display_freqs[idx].mean()))))

lower = numpy.ones_like(rest_freqs)

for idx in reversed(sorted_idxs):
    species_name = species[display_idxs[idx]]
    family_name = bacterial_phylogeny_utils.get_family_name(species_name, genus_family_map)
    phylum_name = bacterial_phylogeny_utils.get_phylum_name(species_name, genus_phylum_map)
    
    print phylum_name, family_name, species_name, idx, display_freqs[idx].mean()
     
    upper = lower
    lower = upper-display_freqs[idx]
    
    if family_name=='Bacteroidaceae':
        colorVal=bacteroides_colors[current_bacteroides_color]
        current_bacteroides_color+=1
    elif family_name=='Rikenellaceae':
        colorVal=alistipes_colors[current_alistipes_color]
        current_alistipes_color+=1
    elif family_name=='Eubacteriaceae':
        colorVal=eubacterium_colors[current_eubacterium_color]
        current_eubacterium_color+=1
    elif family_name=='Lachnospiraceae':
        colorVal=lachno_colors[current_lachno_color]
        current_lachno_color+=1
    else:
        colorVal = other_colors[current_other_color]
        current_other_color+=1
        #line, = freq_axis.plot(ts,-1*numpy.ones_like(ts),'-')
        #colorVal = pylab.getp(line,'color')
    
    if species_name in pretty_species_name_map:
        species_label = pretty_species_name_map[species_name]
    else:
        species_label = figure_utils.get_pretty_species_name(species_name)
    
    legend_axis.plot([-2,-1],[-2,-1],'s',markersize=3,markeredgewidth=0.0, label=species_label,color=colorVal)
    
    # create interleaved version
    background_ts = []
    background_lowers = []
    background_uppers = []
    
    #freq_axis.fill_between(ts, lower, upper,color=colorVal,alpha=0.5,zorder=1)
    
    for t_idx in xrange(0,len(ts)):
        
        
        
        tlower = ts[t_idx]-dts[t_idx]
        tupper = ts[t_idx]+dts[t_idx]
        
        flower = lower[t_idx]
        fupper = upper[t_idx]
        
        background_ts.append(tlower)
        background_ts.append(tupper)
        background_lowers.append(flower)
        background_lowers.append(flower)
        background_uppers.append(fupper)
        background_uppers.append(fupper)
        
        freq_axis.fill_between([tlower,tupper], [flower,flower],[fupper,fupper],color=colorVal,zorder=2,linewidth=0.25)
    
    freq_axis.fill_between(background_ts, background_lowers, background_uppers,color=colorVal,alpha=0.5,zorder=1,linewidth=0.25)
    
    
upper = lower
lower = numpy.zeros_like(rest_freqs)
# create interleaved version
background_ts = []
background_lowers = []
background_uppers = []
    
#freq_axis.fill_between(ts, lower, upper,color=colorVal,alpha=0.5,zorder=1)
    
for t_idx in xrange(0,len(ts)):
        
    tlower = ts[t_idx]-dts[t_idx]
    tupper = ts[t_idx]+dts[t_idx]
        
    flower = lower[t_idx]
    fupper = upper[t_idx]
        
    background_ts.append(tlower)
    background_ts.append(tupper)
    background_lowers.append(flower)
    background_lowers.append(flower)
    background_uppers.append(fupper)
    background_uppers.append(fupper)
    freq_axis.fill_between([tlower,tupper], [flower,flower],[fupper,fupper],color='0.7',zorder=2,linewidth=0.25)
    
freq_axis.fill_between(background_ts, background_lowers, background_uppers,color='0.7',alpha=0.5,zorder=1,linewidth=0.25)
    
legend_axis.plot([-2,-1],[-2,-1],'s',markersize=3,markeredgewidth=0.0, color='0.7',label='Other')

white_ts = numpy.array([ts[0]+(ts[1]-ts[0])/3.0,ts[1]-(ts[1]-ts[0])/3.0])
white_lowers = numpy.array([0,0])
white_uppers = numpy.array([1,1])

print white_ts
print white_lowers
print white_uppers

freq_axis.fill_between(white_ts, white_lowers, white_uppers,color='w',zorder=3)
   
    
#alpha_axis.plot(ts, shannon_diversity,'b.-',linewidth=1.5)   


beta_trajectory = []
for idx in xrange(0,len(ts)):
	
	current_t = ts[idx]
	
	baseline_idxs = (ts<30)*(ts>-1)*(ts!=current_t)
	
	betavals = js_divergence_matrix[idx][baseline_idxs]
	
	beta_trajectory.append(betavals.mean())
	
	beta_axis.plot([current_t, current_t], [betavals.min(), betavals.max()],'r-')
	
beta_axis.plot(ts, beta_trajectory,'r.-',markersize=4,linewidth=1)
	
#beta_axis.plot(ts[js_divergence_matrix[1]>0], js_divergence_matrix[1][js_divergence_matrix[0]>0],'r.-',linewidth=1.5)  
#beta_axis.plot(ts[js_divergence_matrix[4]>0], js_divergence_matrix[4][js_divergence_matrix[3]>0],'r.-',linewidth=1.5,alpha=0.5)  

beta_axis.fill_between(white_ts, [0.01,0.01], [0.4,0.4],color='w',zorder=3)
beta_axis.plot(ts, beta_trajectory,'r:',markersize=4,linewidth=1,zorder=4)


legend_axis.legend(loc='upper center',frameon=False,fontsize=8,numpoints=1,ncol=2,handlelength=1)   



fig.savefig('%s/species_freq_timecourse.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')

# Now do extinction plot
# Set up figure

from cycler import cycler
default_cycler = cycler('color', ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2','#7f7f7f','#bcbd22','#17becf']) 
mpl.rc('axes', prop_cycle=default_cycler)
 

fig = plt.figure(figsize=(5, 5))

# Set up grids to hold figure panels
outer_grid = gridspec.GridSpec(2, 2, width_ratios=[1,0.5],height_ratios=[1,1], wspace=0.1, hspace=0.1)

initial_abx_axis = plt.Subplot(fig, outer_grid[0,0])
fig.add_subplot(initial_abx_axis)

initial_abx_axis.set_ylabel('Min ABX abundance')
initial_abx_axis.loglog([1e-06,1],[1e-06,1],'k-',linewidth=0.25)
initial_abx_axis.fill_between([1e-06,1],[1e-07,1e-07],[1e-07,1e-01],color='0.85') 

initial_abx_axis.set_ylim([7e-06,1])
initial_abx_axis.set_xlim([1e-04,1])
initial_abx_axis.set_xticklabels([])

initial_final_axis = plt.Subplot(fig, outer_grid[1,0])
fig.add_subplot(initial_final_axis)

initial_final_axis.set_ylabel('Max final abundance')
initial_final_axis.set_xlabel('Max baseline abundance')

initial_final_axis.loglog([1e-06,1],[1e-06,1],'k-',linewidth=0.25)
initial_final_axis.fill_between([1e-06,1],[1e-07,1e-07],[1e-07,1e-01],color='0.85') 
initial_final_axis.set_ylim([7e-06,1])
initial_final_axis.set_xlim([1e-04,1])

fig2 = plt.figure(figsize=(3.42, 2))

#final_abx_axis = plt.Subplot(fig, outer_grid[0,1])
#fig.add_subplot(final_abx_axis)
final_abx_axis = fig2.gca()

final_abx_axis.set_xlabel('Max final abundance')

final_abx_axis.loglog([1e-06,1],[1e-06,1],'k-',linewidth=0.25)
final_abx_axis.set_ylim([7e-06,1])
final_abx_axis.set_xlim([1e-04,1])
final_abx_axis.set_yticklabels([])

legend_axis = plt.Subplot(fig, outer_grid[0,1])
fig.add_subplot(legend_axis)

legend_axis.set_ylim([0,1])
legend_axis.set_xlim([0,1])

legend_axis.spines['top'].set_visible(False)
legend_axis.spines['right'].set_visible(False)
legend_axis.spines['left'].set_visible(False)
legend_axis.spines['bottom'].set_visible(False)

legend_axis.set_xticks([])
legend_axis.set_yticks([])  

genus_family_map = bacterial_phylogeny_utils.get_genus_family_map()
genus_phylum_map = bacterial_phylogeny_utils.get_genus_phylum_map()

family_colorsymbol_map = {}
used_colors = set()

sorted_idxs = range(0,len(species))
sorted_idxs = list(sorted(sorted_idxs, key = lambda idx: (species_max_freqs[idx]),reverse=True))

for i in xrange(0,len(species)):
        
    species_coverages = species_coverage_matrix[i,:]
    species_freqs = species_freq_matrix[i,:]
    species_name = species[i]
    
    max_initial = species_freqs[1:5].max()
    min_abx = species_freqs[9:15].min()
    max_final = species_freqs[17:].max()
    
    if max_initial<1e-04 and min_abx>=1e-04:
        continue
    
    if max_initial<1e-04 and max_final>1e-04:
        print "Invader:", species_name, max_initial, max_final
    
    #if species_freqs.max()<1e-04:
    #    continue
        
    if max([max_initial,min_abx,max_final])<1e-04:
        continue
            
    
    family_name = bacterial_phylogeny_utils.get_family_name(species_name, genus_family_map)
    if family_name not in family_colorsymbol_map:
        
        line, = initial_abx_axis.plot([1e-09,1e-09],[1e-09,1e-09],'.',markersize=3,markeredgewidth=0.0)
        colorVal = pylab.getp(line,'color')
        
        if colorVal not in used_colors:
            symbol='o'
            used_colors.add(colorVal)
        else:
            symbol='s'
        
        line, = legend_axis.plot([-2,-1],[-2,-1],symbol,markersize=3,markeredgewidth=0.0, label=family_name)
        
        family_colorsymbol_map[family_name] = (colorVal,symbol)
    
    colorVal,symbol = family_colorsymbol_map[family_name]
    
    
    clipped_max_initial=max([1e-04,max_initial])
    clipped_min_abx=max([1e-05,min_abx])
    clipped_max_final=max([1e-05,max_final])
    
    if min_abx<1e-05:
        markersize=2
        color='none'
        markeredgewidth=0.5
        markeredgecolor=colorVal
    else:
        markersize=2
        color=colorVal
        markeredgewidth=0
        markeredgecolor=colorVal
    
    initial_abx_axis.plot([clipped_max_initial],[clipped_min_abx],symbol,markersize=markersize,markeredgewidth=markeredgewidth, color=color,markeredgecolor=markeredgecolor)
    
    if max_final<1e-05:
        markersize=2
        color='none'
        markeredgewidth=0.5
        markeredgecolor=colorVal
        print species_name, max_initial, min_abx, max_final
    else:
        markersize=2
        color=colorVal
        markeredgewidth=0
        markeredgecolor=colorVal
    
 

    initial_final_axis.plot([clipped_max_initial],[clipped_max_final],symbol,markersize=markersize,markeredgewidth=markeredgewidth, color=color,markeredgecolor=markeredgecolor)
    
legend_axis.legend(loc='lower left',frameon=False,fontsize=5,numpoints=1,ncol=1,handlelength=1)   
    
fig.savefig('extinction_figure.pdf',bbox_inches='tight')