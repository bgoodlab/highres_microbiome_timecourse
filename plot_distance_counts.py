import pylab
import numpy
import config
import figure_utils
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint

mpl.rcParams['font.size'] = 7
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

#desired_speciess = ['Bacteroides_vulgatus_57955'] #,'Bacteroides_coprocola_61586','Bacteroides_uniformis_57318'] 

desired_speciess = ['Alistipes_finegoldii_56071','Alistipes_sp_60764','Alistipes_onderdonkii_55464','Eubacterium_eligens_61678','Phascolarctobacterium_sp_59817','Bacteroides_coprocola_61586']

pylab.figure(figsize=(2,1.25))

for desired_species in desired_speciess:

    filename = "%s/new_barcode_trajectories/%s_distance_counts.txt" % (config.barcode_directory, desired_species)
 
    file = open(filename,"r")
    line = file.readline()
    items = line.split(",")
    distances = numpy.array([float(item) for item in items])
    line = file.readline()
    items = line.split(",")
    Bs = numpy.array([float(item) for item in items])
    line = file.readline()
    items = line.split(",")
    Ss = numpy.array([float(item) for item in items])

    distance_bins = numpy.hstack([numpy.arange(0,100),     numpy.logspace(2,7,101)])

    distance_centers =     numpy.array(distance_bins[0:-1],copy=True)
    distance_bins-=0.5

    distance_idxs = numpy.digitize(distances,bins=distance_bins)

    binned_Bs = numpy.zeros_like(distance_centers)*1.0
    binned_Ss = numpy.zeros_like(binned_Bs)

    numpy.add.at(binned_Bs, distance_idxs, Bs)
    numpy.add.at(binned_Ss, distance_idxs, Ss)

    good_idxs = (binned_Bs>1000)
    distance_centers = distance_centers[good_idxs]
    binned_Ss = binned_Ss[good_idxs]
    binned_Bs = binned_Bs[good_idxs]


    pretty_species_name = figure_utils.get_pretty_species_name(desired_species)
    pylab.loglog(distance_centers, binned_Ss/binned_Bs,'-',label=pretty_species_name,color=figure_utils.example_species_color_map[desired_species])


pylab.xlim([1,1e05])
pylab.ylim([1e-02,1])
pylab.legend(loc=(0.9,0.2),frameon=False,numpoints=1)
pylab.gca().get_xaxis().tick_bottom()
pylab.gca().get_yaxis().tick_left()
pylab.gca().spines['top'].set_visible(False)
pylab.gca().spines['right'].set_visible(False)
pylab.xlabel('Coordinate distance, $\ell$ (bp)')
pylab.ylabel('Fraction shared read clouds')
pylab.savefig('barcode_sharing_distance_counts.pdf',bbox_inches='tight')