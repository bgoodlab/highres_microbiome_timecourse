import config
import parse_timecourse_data
import parse_midas_data
import gzip
import pylab
import numpy
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

mpl.rcParams['font.size'] = 7
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

antibiotic_directory = "/Users/bgood/highres_microbiome_timecourse_barcode_data/antibiotic_barcodes/"

samples = parse_timecourse_data.morteza_samples

species_coverage_matrix, samples, species = parse_midas_data.parse_global_marker_gene_coverages()
sample_time_map = parse_timecourse_data.parse_sample_time_map()
ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, samples)
samples = numpy.array(samples)[sample_idxs]

antibiotic_timecourse = {}

linked_genes = ['gi|511019447|ref|WP_016273783.1|', 'gi|1045851695|ref|WP_065538049.1|', 'AAM48122', 'AAM48121', 'AAM48125', 'AAM48124', 'U38243.1.gene1.p01', 'AY769933.1.gene1.p01', 'AF118110.1.gene1.p01', 'AF472622.2.gene1.p01', 'AAV37206', 'AAM48117', 'AAM48116', 'P30898', 'YP_098534']

for sample in samples:
    
    file = gzip.open("%s/%s_antibiotic_barcodes.txt.gz" % (antibiotic_directory, sample),"r")
    for line in file:
        items = line.split("\t")
        gene_name = items[0].strip()
        subtype = items[1].strip()
        antibiotic = subtype.split("__")[0]
        ppm = float(items[2])
        
        if gene_name in linked_genes:
            antibiotic = 'beta-lactam (tet-linked)'
        
        #antibiotic = subtype
        if antibiotic not in antibiotic_timecourse:
            antibiotic_timecourse[antibiotic] = {}
            
        if sample not in antibiotic_timecourse[antibiotic]:
            antibiotic_timecourse[antibiotic][sample] = 0
        
        antibiotic_timecourse[antibiotic][sample] += ppm
 
            
new_antibiotic_timecourse = {}
for antibiotic in antibiotic_timecourse:
    new_antibiotic_timecourse[antibiotic] = []
    for sample in samples:
        
        if sample in antibiotic_timecourse[antibiotic]:
            new_antibiotic_timecourse[antibiotic].append( antibiotic_timecourse[antibiotic][sample])
        else:
            new_antibiotic_timecourse[antibiotic].append(0)
    
    new_antibiotic_timecourse[antibiotic] = numpy.array(new_antibiotic_timecourse[antibiotic])    

antibiotic_timecourse = new_antibiotic_timecourse

pylab.figure(figsize=(7,3))
for antibiotic in antibiotic_timecourse:

    if antibiotic_timecourse[antibiotic].max() > 10 or (antibiotic.startswith('beta')):

        pylab.plot(ts, antibiotic_timecourse[antibiotic],'.-',label=antibiotic)

pylab.fill_between([parse_timecourse_data.antibiotic_start, parse_timecourse_data.antibiotic_end],[0,0],[600,600],color=parse_timecourse_data.antibiotics_color,linewidth=0)
        
pylab.ylim([0,550])
pylab.legend(frameon=False,loc='upper left')
pylab.ylabel('Reads (per million) mapped to ABX genes')
pylab.xlabel('Day')
pylab.savefig('antibiotic_ppm_timecourse.pdf',bbox_inches='tight')
