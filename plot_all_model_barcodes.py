import sys
import cPickle
import numpy

#####
#
# Load precalculated fixation events
#
####


# create output filename
pickle_filename = sys.argv[1]

print "Created by plot_all_model_barcodes.py"
print "from", pickle_filename

d = cPickle.load( open( pickle_filename,"rb") )

# Compile site list
for feature in sorted(d):
    pass
    
for feature in sorted(d):
    
    focal_species = d[feature]['focal_species'] 
    focal_alleles = [a for a in d[feature]['focal_alleles']]
    allele_str = ", ".join(focal_alleles)
        
    all_barcodes = d[feature]['all']
    
    if all_barcodes.sum() < 2.5:
        continue
    
    print "Focal feature: %s, %s" % (focal_species, allele_str)
    print "Total barcodes: %s" % (", ".join([str(n) for n in all_barcodes]))
    
    # Sort things by # of barcodes
    other_species_list = d[feature]['species'].keys()
    
    species_gene_list = []
    species_gene_barcodes = []
    
    for other_species_name in sorted(d[feature]['species']):
        for other_gene_name in sorted(d[feature]['species'][other_species_name]):
                    
            barcodes, pvalues = d[feature]['species'][other_species_name][other_gene_name]
            total_barcodes = barcodes.sum()
            species_gene_list.append((other_species_name, other_gene_name))
            species_gene_barcodes.append(total_barcodes)
    
    # Sort in descending order
    if len(species_gene_list)>0:
        species_gene_barcodes, species_gene_list = (list(x) for x in zip(*sorted(zip(species_gene_barcodes, species_gene_list), key=lambda pair: pair[0], reverse=True)))
    
    # only look at the top 5 entries
    for other_idx in xrange(0,len(species_gene_list)):
        other_species_name, other_gene_name = species_gene_list[other_idx]
        
        is_snp = ((other_gene_name[-1]=='A') or (other_gene_name[-1]=='R'))
        
        if (not is_snp) and other_idx>=10:
            continue
        
        if other_species_name==focal_species:
            linkage_label = ( "%s->self" % (focal_species) )
        else:
            linkage_label = ("%s->%s" % (focal_species, other_species_name))
        
        barcodes, pvalues = d[feature]['species'][other_species_name][other_gene_name]
            
        if (barcodes>0).sum()<3:
            continue
            
        print "%s %s (%d): %s" % (linkage_label, other_gene_name, barcodes.sum(), ", ".join(["%d" % (n) for n,ntot,p in zip(barcodes, all_barcodes, pvalues)]))    
    
    print "-"        
    print ""
    print ""
                