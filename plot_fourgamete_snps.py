import sys
import pylab
import numpy
from math import fabs,log10,log
import config
#####
#
# Load trajectories 
#
#####
# create output filename
species_name = None
debug=False
chunk_size=100000000

Pstar = 0.1


# Things with a lot of 4 gametes...
species_names = ['Faecalibacterium_cf_62236']
#species_names = ['Blautia_wexlerae_56130']
#species_names = [species_names[0]]
#species_names = ['Alistipes_sp_60764']
#species_names = ['Bacteroides_vulgatus_57955']
#speciess_names = ['Eubacterium_eligens_61678']
#species_names = ['Alistipes_finegoldii_56071']
species_trajectory_distances = {}

contig_id_map = {}
id_contig_map = []

filename = ('new_cross_species_snp_barcode_statistics_core.txt')

file = open(filename,"r")

linkage_data = {}
snp_coordinates = set()

for line in file:    
    items = line.split()
    species_name = items[0]
    
    if species_name not in species_names:
        continue
    
    testable = long(items[1])
    onegametes = long(items[2])
    twogametes = long(items[3])
    threegametes = long(items[4])
    testable_fourgametes = long(items[5])
    fourgametes = long(items[6])
    
    if len(items)<8:
        continue
        
    fourgamete_items = items[7:]
    for item in fourgamete_items:
        
        subitems = item.split(",")
        focal_snp_items = subitems[0].strip().split("|")
        target_snp_items = subitems[1].strip().split("|")
    
    
        focal_contig = focal_snp_items[0]
        focal_location = focal_snp_items[1]
        
        target_contig = target_snp_items[0]
        target_location = target_snp_items[1]
        
        if focal_contig!=target_contig:
            print "Weird thing happened!"
            continue
        
        if (species_name, focal_contig) not in contig_id_map:
            contig_id = len(id_contig_map)
            id_contig_map.append((species_name, focal_contig))
            contig_id_map[(species_name, focal_contig)] = contig_id
        
        contig_id = contig_id_map[(species_name, focal_contig)]
        
        snp_coordinates.add((focal_location, contig_id+1))
        snp_coordinates.add((target_location, contig_id+1))
        
        
        linkage_data[species_name] = (testable,onegametes,twogametes,threegametes,testable_fourgametes,fourgametes)
    

    species_list = linkage_data.keys()
    species_list.sort()
file.close()

sys.stderr.write("Plotting %d SNVs!\n" % len(snp_coordinates))

num_configs = len(id_contig_map)
pylab.figure(1,figsize=(21,2))
for x,y in snp_coordinates:
    pylab.plot([x],[y],'bo',markersize=3,markeredgewidth=0,alpha=0.5)

pylab.yticks(range(1,len(id_contig_map)+1),[contig[1] for contig in id_contig_map])    
pylab.ylim([0,len(id_contig_map)+1])  
pylab.savefig('fourgamete_snp_locations.pdf',bbox_inches='tight')

