import config
import gzip 
import numpy
import os.path

gamete_idx_map = {('R','R'): 0, ('R','A'):1, ('A','R'):2, ('A','A'):3  }
allele_idx_map = {'R':0, 'A':1}

# Returns list of barcodes and what alleles they were found on.

base_idx_map = {('R','R'):0, ('A','A'):1, ('R','A'):2, ('A','R'):3}

def parse_tracked_snps(desired_speciess = set()):
	
	tracked_longsnps = set()
	
	filename = "%s/desired_snps.txt.gz" % (config.barcode_directory)
	file = gzip.GzipFile(filename,"r")
	file.readline() # header
	for line in file:
		items = line.split()
		species_name = items[0].strip()
		contig = items[1].strip()
		location = long(items[2])
		longsnp = (species_name, contig, location)
		if len(desired_speciess)==0 or (species_name in desired_speciess):
			tracked_longsnps.add(longsnp)
	file.close()   
	
	return tracked_longsnps

def calculate_linkage_score(gamete_vector, base_1='A', base_2='A'):
	
	n = gamete_vector.sum()
	
	linkage_score = (gamete_vector[0]+gamete_vector[1])*1.0/n - (gamete_vector[2]+gamete_vector[3])*1.0/n		 

	if base_1!=base_2:
		linkage_score *= -1
		
	return linkage_score

def minimum_gamete_fraction(gamete_vector):
	
	n = gamete_vector.sum()
	return gamete_vector.min()*1.0/n


def parse_all_mapped_barcodes(sample_name, corrected=True, min_depth=0):

	if corrected:
		barcode_filename = "%s%s/output/all_corrected_barcodes.gz" % (config.barcode_directory, sample_name)
	
	else:
		barcode_filename = "%s%s/output/all_barcodes.gz" % (config.barcode_directory, sample_name)
		
	barcode_depth_map = {} # map from barcode id->depth
	barcode_id_map = {} # map from (corrected) barcode str -> barcode_id
	
	if os.path.isfile(barcode_filename):
		barcode_file = gzip.GzipFile(barcode_filename,"r")
		barcode_file.readline() # header
	
		for line in barcode_file:
			items = line.split()
			barcode_id = long(items[0])
			barcode_str = items[1].strip()
			barcode_weight = long(items[2])
			
			if barcode_weight >= min_depth:
				barcode_depth_map[barcode_id] = barcode_weight
				barcode_id_map[barcode_str] = barcode_id
				
		barcode_file.close()
	
	return barcode_depth_map, barcode_id_map


def parse_total_barcode_map():
	
	barcode_filename = "%s/barcode_summary.txt" % (config.barcode_directory)

	total_barcode_map = {}

	file = open(barcode_filename,"r")
	file.readline() # header
	for line in file:
		items = line.split()
		sample_name = items[0].strip()
		unfiltered_total_barcodes = long(items[1])
		unfiltered_barcode_depth = long(items[2])
		total_barcodes = long(items[3])
		barcode_depth = long(items[4])
		
		total_barcode_map[sample_name] = (total_barcodes, barcode_depth)
		
	return total_barcode_map
		
def parse_barcode_depth_map(sample_name, corrected=True, min_depth=0):

	if corrected:
		barcode_filename = "%s%s/output/all_corrected_barcodes.gz" % (config.barcode_directory, sample_name)
	
	else:
		barcode_filename = "%s%s/output/all_barcodes.gz" % (config.barcode_directory, sample_name)
		
	barcode_depth_map = {}
	
	if os.path.isfile(barcode_filename):
		barcode_file = gzip.GzipFile(barcode_filename,"r")
		barcode_file.readline() # header
	
		for line in barcode_file:
			items = line.split()
			barcode_id = long(items[0])
			barcode_weight = long(items[2])
			
			if barcode_weight >= min_depth:
				barcode_depth_map[barcode_id] = barcode_weight
			
		barcode_file.close()
	
	return barcode_depth_map

def parse_raw_barcode_error_correction_map(sample_name):
	barcode_filename = "%s/barcode_error_maps/%s.error_map.gz" % (config.barcode_directory, sample_name)
	barcode_error_map = {}
	barcode_depth_map = {}
	
	if not os.path.isfile(barcode_filename):
		return barcode_error_map
	
	barcode_file = gzip.GzipFile(barcode_filename,"r")
	barcode_file.readline() # header
	
	for line in barcode_file:
		items = line.split()
		original_barcode_str = items[0].strip()
		corrected_barcode_str = items[1].strip()
		num_reads = long(items[2])
	
		if corrected_barcode_str not in barcode_depth_map:
			barcode_depth_map[corrected_barcode_str] = 0
		
		barcode_error_map[original_barcode_str] = corrected_barcode_str
		barcode_depth_map[corrected_barcode_str] += num_reads
	
	barcode_file.close()
	return barcode_error_map, barcode_depth_map	  

def calculate_barcode_blacklist(barcode_str_iterator):
	
	barcode_blacklist = set()
	
	for barcode_str in barcode_str_iterator:
		
		ntot = len(barcode_str)
		nA = barcode_str.count('A')
		nT = barcode_str.count('T')
		nC = barcode_str.count('C')
		nG = barcode_str.count('G')
		nN = barcode_str.count('N')
		
		
		
		
		
		#b = numpy.array(list(barcode_str))
		#ntot = len(b)
		#nA = (b=='A').sum()
		#nC = (b=='C').sum()
		#nT = (b=='T').sum()
		#nG = (b=='G').sum()
		#nN = (b=='N').sum()
		
		nmax = max([nA,nC,nT,nG])
		
		if nmax>=0.8*ntot:
			barcode_blacklist.add(barcode_str)
			
		if nN>=0.5*ntot:
			barcode_blacklist.add(barcode_str)
		
	return barcode_blacklist
		
		
def parse_barcode_error_correction_map(sample_name):

	barcode_filename = "%s%s/output/barcode_map.gz" % (config.barcode_directory, sample_name)
	barcode_error_map = {}
		
	if not os.path.isfile(barcode_filename):
		return barcode_error_map
	
	barcode_file = gzip.GzipFile(barcode_filename,"r")
	barcode_file.readline() # header
	
	for line in barcode_file:
		items = line.split()
		original_barcode = long(items[0])
		corrected_barcode = long(items[1])
		
		barcode_error_map[original_barcode] = corrected_barcode
		
	barcode_file.close()
	return barcode_error_map
	
	
def barcodes_exist(species_name, sample_name, corrected=True):
	if corrected:
		barcode_filename = "%s%s/output/%s.corrected_barcodes.gz" % (config.barcode_directory, sample_name, species_name)
	
	else:
		barcode_filename = "%s%s/output/%s.barcodes.gz" % (config.barcode_directory, sample_name, species_name)
		
	if os.path.isfile(barcode_filename):
		return True
	else:
		return False
		

def parse_uncorrected_allele_barcode_tuples(species_name, sample_name, corrected=True, bootstrapped=False):
	
	barcode_filename = "%s%s/output/%s" % (config.barcode_directory, sample_name, species_name)
	if corrected:
		barcode_filename += ".corrected_barcodes"
	else:
		barcode_filename += ".barcodes"
		
	if bootstrapped:
		barcode_filename += ".bootstrapped"
		
	barcode_filename+=".gz"
	
	barcode_file = gzip.GzipFile(barcode_filename,"r")
	barcode_file.readline()
	
	allele_barcode_map = {}
	
	for line in barcode_file:
		line = line.strip()
		items = line.split("\t")
		allele = items[0].strip()
		allele_barcode_map[allele] = []
		
		if len(items)>1:
		
			barcode_items = items[1].split(",")
			for barcode_item in barcode_items:
				barcode_subitems = barcode_item.split(":")
				barcode_id = long(barcode_subitems[0])
				barcode_weight = long(barcode_subitems[1])
				
				allele_barcode_map[allele].append((barcode_id, barcode_weight))
				
	return allele_barcode_map

# Like parse allele barcode tuples,
# but use multiple reads per barcode to correct sequencing errors for SNVs
def parse_allele_barcode_tuples(species_name, sample_name, corrected=True, bootstrapped=False):
	
	barcode_filename = "%s%s/output/%s" % (config.barcode_directory, sample_name, species_name)
	if corrected:
		barcode_filename += ".corrected_barcodes"
	else:
		barcode_filename += ".barcodes"
		
	if bootstrapped:
		barcode_filename += ".bootstrapped"
		
	barcode_filename+=".gz"
	
	barcode_file = gzip.GzipFile(barcode_filename,"r")
	barcode_file.readline()
	
	allele_barcode_map = {}
	allele_error_map = {}
	
	line = barcode_file.readline().strip()
	while line!="":
		
		items = line.split("\t")
		allele = items[0].strip()
		allele_barcode_map[allele] = []
		allele_error_map[allele] = []
		
		if allele[-1]!='R':
			# Loaded a gene proceed normally
			if len(items)>1:
		
				barcode_items = items[1].split(",")
				for barcode_item in barcode_items:
					barcode_subitems = barcode_item.split(":")
					barcode_id = long(barcode_subitems[0])
					barcode_weight = long(barcode_subitems[1])
				
					allele_barcode_map[allele].append((barcode_id, barcode_weight))
		
		else:
			# Loading a SNP. Parse two lines at once
			line2 = barcode_file.readline().strip()
			items2 = line2.split("\t")
			allele2 = items2[0].strip()
			allele_barcode_map[allele2] = []
			allele_error_map[allele2] = []
			
			items1 = items
			allele1 = allele
			
			barcode_weight_map = {}
			
			# first go through first allele (R)	
			if len(items1)>1:
					
				barcode_items = items1[1].split(",")
			
				for barcode_item in barcode_items:
					barcode_subitems = barcode_item.split(":")
					barcode_id = long(barcode_subitems[0])
					barcode_weight = long(barcode_subitems[1])
					barcode_weight_map[barcode_id] = barcode_weight
					
			# Now go through second allele (A)
			if len(items2)>1:
					
				barcode_items = items2[1].split(",")
			
				for barcode_item in barcode_items:
					barcode_subitems = barcode_item.split(":")
					barcode_id = long(barcode_subitems[0])
					barcode_weight = long(barcode_subitems[1])
					
					if barcode_id not in barcode_weight_map:
						# No collision with allele1. 
						# Add directly to output for allele2
						allele_barcode_map[allele2].append((barcode_id, barcode_weight))
					else:
						# Collision with allele1. 
						# Do error correction:
						if barcode_weight < barcode_weight_map[barcode_id]:
							# Correct to allele1
							error_weight = barcode_weight
							barcode_weight_map[barcode_id] += error_weight
							allele_error_map[allele1].append((barcode_id, error_weight))
						else:
							# Correct to allele2
							error_weight = barcode_weight_map[barcode_id]
							corrected_weight = barcode_weight+error_weight
							allele_barcode_map[allele2].append((barcode_id, corrected_weight))
							allele_error_map[allele2].append((barcode_id, error_weight))
							# Now remove from database
							del barcode_weight_map[barcode_id]
							
			# What's left now can be added to allele1											
			for barcode_id in barcode_weight_map:
				allele_barcode_map[allele1].append((barcode_id, barcode_weight_map[barcode_id]))
				
		line = barcode_file.readline().strip()	 
				
	return allele_barcode_map, allele_error_map
	
def parse_allele_barcodes(species_name, sample_name, corrected=True):
	
	if corrected:
		barcode_filename = "%s%s/output/%s.corrected_barcodes.gz" % (config.barcode_directory, sample_name, species_name)
	
	else:
		barcode_filename = "%s%s/output/%s.barcodes.gz" % (config.barcode_directory, sample_name, species_name)
	
	barcode_file = gzip.GzipFile(barcode_filename,"r")
	barcode_file.readline()
	
	allele_barcode_map = {}
	
	for line in barcode_file:
		line = line.strip()
		items = line.split("\t")
		allele = items[0].strip()
		allele_barcode_map[allele] = []
		
		if len(items)>1:
		
			barcode_items = items[1].split(",")
			for barcode_item in barcode_items:
				barcode_subitems = barcode_item.split(":")
				barcode_id = long(barcode_subitems[0])
				
				allele_barcode_map[allele].append(barcode_id)
				
	return allele_barcode_map
	
def calculate_barcode_allele_map(allele_barcode_map):

	barcode_allele_map = {}
	for allele in allele_barcode_map.keys():
		
		for barcode in allele_barcode_map[allele]:
			
			if barcode not in barcode_allele_map:
				barcode_allele_map[barcode] = []
			
			barcode_allele_map[barcode].append(allele)
	
	return barcode_allele_map

def calculate_num_shared_barcodes(allele_barcode_map, barcode_allele_map, desired_alleles):
	
	num_shared_barcodes = {allele: {} for allele in desired_alleles}
	for allele in desired_alleles:
		
		if allele not in allele_barcode_map:
			continue
			
		barcodes = allele_barcode_map[allele]
		
		if len(barcodes) == 0:
			continue
			
		for barcode in barcodes:
			for other_allele in barcode_allele_map[barcode]:
				if other_allele not in num_shared_barcodes[allele]:
					num_shared_barcodes[allele][other_allele] = 0
					
				num_shared_barcodes[allele][other_allele] += 1
				
				
	return num_shared_barcodes
	
def calculate_num_shared_barcodes_per_site(allele_barcode_map, barcode_allele_map, desired_sites=set([]), num_shared_barcodes_per_site = {}):
	
	desired_sites_set = set(desired_sites)
	
	allele_items_map = {}
	
	for allele_1 in allele_barcode_map.keys():
		
		if not ((allele_1[-1]=='A') or (allele_1[-1]=='R')):
			# not a SNP
			continue	
		
		if not allele_1 in allele_items_map:
			
			allele_items_1 = allele_1.split("|")
			contig_1 = allele_items_1[0]
			position_1 = long(allele_items_1[1]) 
			base_1 = allele_items_1[2].strip()
			site_1 = (contig_1, position_1)
			allele_items_map[allele_1] = (site_1, base_1)
			
		site_1, base_1 = allele_items_map[allele_1]
		
		if site_1 not in desired_sites_set:
			continue
		
		if site_1 not in num_shared_barcodes_per_site:
			num_shared_barcodes_per_site[site_1] = {}
	
		barcodes = allele_barcode_map[allele_1]
		
		if len(barcodes)==0:
			# only linked to itself!
			continue
			
		for barcode in barcodes:
		
			if len(barcode_allele_map[barcode]) < 2:
				# only linked to itself!
				continue
				
			for allele_2 in barcode_allele_map[barcode]:
				
				if not ((allele_2[-1]=='A') or (allele_2[-1]=='R')):
					# not a SNP
					continue
			
				if not allele_2 in allele_items_map:
		
					allele_items_2 = allele_2.split("|")
					contig_2 = allele_items_2[0]
					position_2 = long(allele_items_2[1]) 
					base_2 = allele_items_2[2].strip()
					site_2 = (contig_2, position_2)
					allele_items_map[allele_2] = (site_2, base_2)
			
				site_2, base_2 = allele_items_map[allele_2]
		 
				if site_2 not in num_shared_barcodes_per_site[site_1]:
					num_shared_barcodes_per_site[site_1][site_2] = numpy.zeros(len(base_idx_map))
		
				num_shared_barcodes_per_site[site_1][site_2][base_idx_map[(base_1, base_2)]] += 1
				
	return num_shared_barcodes_per_site

	
def calculate_linked_set(num_shared_barcodes, min_shared_barcodes=1):
	
	linked_set = set([])
	linked_set.update(num_shared_barcodes) 
	
	for allele in num_shared_barcodes:
		for other_allele in num_shared_barcodes[allele].keys():
			if num_shared_barcodes[allele][other_allele] >= min_shared_barcodes:
				linked_set.add(other_allele)
				
	return linked_set
	
#####
#
# Function used to calculate mean number of fragments from barcode distribution
#
#####	 
def calculate_M_from_D(D):
	
	if D<10:
		return D
	else:
		return 10

def parse_raw_G_function_map():
	 
	filename = "%s/barcode_G_function.txt" % (config.barcode_directory)
	 
	file = open(filename,"r")
	line = file.readline()
	items = line.split(",")
	Dbins = numpy.array([float(item) for item in items])
	 
	raw_G_function_map = {}
	 
	for line in file:
		items = line.split(",")
		sample_name = items[0].strip()
		 
		G_numerators = []
		G_denominators = []
		for item in items[1:]:
			subitems = item.split("|")
			G_numerators.append(float(subitems[0]))
			G_denominators.append(float(subitems[1]))
		
		raw_G_function_map[sample_name] = (numpy.array(G_numerators), numpy.array(G_denominators))
					 
	return Dbins, raw_G_function_map
  
def G_function_factory(Dbins, G_function):
	return lambda x: G_function[(x>=Dbins[0:-1])*(x<Dbins[1:])]
	
# Smooths it out for us
def parse_G_function_map():
	
	Dbins, raw_G_function_map = parse_raw_G_function_map()
	
	G_function_map = {}
	
	for sample_name in raw_G_function_map:

		numerators, denominators = raw_G_function_map[sample_name]

		max_idx = denominators.argmax()
		
		G_function = numerators*1.0/(denominators+(denominators==0))
		
		left_side = (numpy.arange(0,len(denominators))<max_idx)
		
		right_side = (numpy.arange(0,len(denominators))>max_idx)

		left_bad_idxs = numpy.nonzero(left_side*(denominators<10000))[0]
		
		if len(left_bad_idxs)==0:
			left_good_idx = 0
		else:
			left_good_idx = left_bad_idxs[-1]+1
			
		right_bad_idxs = numpy.nonzero(right_side*(denominators<10000))[0]
		
		if len(right_bad_idxs)==0:
			right_good_idx = len(denominators)-1
		else:
			right_good_idx = right_bad_idxs[0]-1
		
		left_good_G = G_function[left_good_idx]
		right_good_G = G_function[right_good_idx]
		
		if left_good_idx>0:
			G_function[:left_good_idx] = left_good_G
		
		if right_good_idx<(len(denominators)-1):
			G_function[right_good_idx:] = right_good_G
		
		
		functional_G_function = G_function_factory(Dbins,G_function)
		
		G_function_map[sample_name] = functional_G_function
 
	return G_function_map