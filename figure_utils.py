
    

def get_pretty_species_name(species_name, include_number=False):
    
    items = species_name.split("_")
    
    pretty_name = "%s %s" % (items[0], items[1])
    
    if include_number:
        pretty_name += (" (%d)" % (items[2]))
        
    return pretty_name

def get_pretty_species_names(species_list):
    
    pretty_species_name_counts = {}
    pretty_species_name_idxs = {}
 
    pretty_species_names = []
    
    for species_name in species_list:
        
        pretty_species_name = get_pretty_species_name(species_name)
        
        pretty_species_names.append(pretty_species_name)
        
        if pretty_species_name not in pretty_species_name_counts:
            pretty_species_name_counts[pretty_species_name] = 0
            pretty_species_name_idxs[pretty_species_name] = 1
            
        pretty_species_name_counts[pretty_species_name] += 1
    
    
    for idx in xrange(0,len(pretty_species_names)):
        
        pretty_species_name = pretty_species_names[idx]
        
        if pretty_species_name_counts[pretty_species_name] > 1.5:
            
            new_pretty_species_name = pretty_species_name + (" %d" % pretty_species_name_idxs[pretty_species_name])

            pretty_species_name_idxs[pretty_species_name] += 1
        
        else:
            
            new_pretty_species_name = pretty_species_name
            
        pretty_species_names[idx] = new_pretty_species_name
        
    return pretty_species_names
        
def get_abbreviated_species_name(species_name):
    
    items = species_name.split("_")
    
    pretty_name = "%s. %s" % (items[0][0], items[1])
        
    return pretty_name

       

    
    
example_species_color_map = {
'Bacteroides_vulgatus_57955': '#000000',
'Alistipes_finegoldii_56071': '#2171b5',
'Alistipes_sp_60764': '#6baed6',
'Alistipes_onderdonkii_55464': '#bdd7e7',
'Eubacterium_eligens_61678': '#238b45',
'Phascolarctobacterium_sp_59817': '#fd8d3c',
'Bacteroides_coprocola_61586': '#a50f15'
}

