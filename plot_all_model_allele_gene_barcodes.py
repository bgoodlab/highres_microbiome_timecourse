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

print "Created by plot_all_model_allele_gene_barcodes.py"
print "from", pickle_filename

d = cPickle.load( open( pickle_filename,"rb") )

data = {}     # map from allele_base_str to linked thing to allele to target.
    
for feature in sorted(d):
    
    focal_species = d[feature]['focal_species'] 
    focal_alleles = [a for a in d[feature]['focal_alleles']]
    total_allele_str = focal_alleles[0]

    if (total_allele_str[-1]=='A') or (total_allele_str[-1]=='R'):
        allele_str = "|".join(total_allele_str.split("|")[:-1])
        allele = total_allele_str[-1]
    else:
        allele_str = total_allele_str
        allele = ""
            
    all_barcodes = d[feature]['all']
    
    if all_barcodes.sum() < 2.5:
        continue
    
    source = (focal_species, allele_str)
    
    if source not in data:
        data[source] = {'total': {}, 'targets': {}}
        
    data[source]['total'][allele] = all_barcodes
    
    for target_species_name in sorted(d[feature]['species']):
        for target_total_allele_str in sorted(d[feature]['species'][target_species_name]):
            
            if (target_total_allele_str[-1]=='A') or (target_total_allele_str[-1]=='R'):
                target_allele_str = "|".join(target_total_allele_str.split("|")[:-1])
                target_allele = target_total_allele_str[-1]
            else:
                target_allele_str = target_total_allele_str
                target_allele = ""  
            
            if target_allele_str == allele_str:
                continue  
              
            target = (target_species_name, target_allele_str)  
            
            barcodes, pvalues = d[feature]['species'][target_species_name][target_total_allele_str]
            
            if target not in data[source]['targets']:
                data[source]['targets'][target] = {}
            
            data[source]['targets'][target][(allele, target_allele)] = barcodes
    
    
for source in sorted(data):

    # Sort targets by total
    target_list = []
    target_total_barcodes = []
    
    for target in sorted(data[source]['targets']):
        
        total_barcodes = 0
        for allele_tuple in data[source]['targets'][target]:
            total_barcodes += data[source]['targets'][target][allele_tuple].sum()
            
        target_list.append(target)
        target_total_barcodes.append(total_barcodes)
        
    if len(target_list)==0:
        continue
    
    target_total_barcodes, target_list = (list(x) for x in zip(*sorted(zip(target_total_barcodes, target_list), key=lambda pair: pair[0], reverse=True)))
    
    # Now print stuff!
    focal_species, focal_site = source
    print "Focal feature: %s, %s" % (focal_species, focal_site)
    print "Total barcodes:"
    for source_allele in sorted(data[source]['total']):
        print "%s: %s" % (source_allele, ", ".join([str(n) for n in data[source]['total'][source_allele]]))
    
    for target_idx in xrange(0,len(target_list)):
        
            
        if target_total_barcodes[target_idx] < 2.5:
            continue
        
        is_snp = ("|" in target_list[target_idx][1])
        # Only print out top 5 hits (+ all other SNVs)  
        if (not is_snp) and target_idx>=5:
            continue
        
        
        target = target_list[target_idx]
        target_species, target_allele_str = target_list[target_idx]
        
        if target_species==focal_species:
            linkage_label = ( "%s->self" % (focal_species) )
        else:
            linkage_label = ("%s->%s" % (focal_species, target_species))
        
        print "%s %s (%d)" % (linkage_label, target_allele_str, target_total_barcodes[target_idx])
        
        for allele_tuple in sorted(data[source]['targets'][target]):
            print "%s%s: %s" % (allele_tuple[0], allele_tuple[1], ", ".join(["%d" % (n) for n in data[source]['targets'][target][allele_tuple]]))        
        
    print "-"        
    print ""
    print ""
                