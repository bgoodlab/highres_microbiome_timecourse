import sys
import cPickle
import numpy
from math import fabs,log10,log
import config
import parse_midas_data
import parse_timecourse_data
from numpy.random import random_sample
from scipy.stats import fisher_exact, chi2_contingency
import os.path
from stats_utils import calculate_poisson_logP
import core_gene_utils
import statsmodels.api as sm
import barcode_utils
from math import exp
   

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

species_names = parse_midas_data.parse_good_species_list()

# Things with a lot of 4 gametes...
#species_names = ['Faecalibacterium_cf_62236']
#species_names = ['Blautia_wexlerae_56130']
#species_names = [species_names[0]]
#species_names = ['Alistipes_sp_60764']
species_trajectory_distances = {}

output_file = open("new_cross_species_snp_barcode_statistics_nonshared.txt","w")

sys.stderr.write("Loading tracked SNVs...\n")
tracked_shortsnps = barcode_utils.parse_tracked_snps()
sys.stderr.write("Done!\n")

inconsistency_pvalue = 0.05
complete_FDR = 0.05
fourgamete_FDR = 0.05

#species_names = ['Alistipes_sp_60764']

for species_name in species_names:

	pickle_filename = ("%s/new_barcode_trajectories/%s_snp_barcodes.p" % (config.barcode_directory, species_name))
		
	if not os.path.isfile(pickle_filename):
		continue
	
	sys.stderr.write("Processing %s...\n" % species_name) 
	#keep to actual core genome	  
	#core_genes = core_gene_utils.parse_core_genes(species_name)
	# Keep to resident genome
	personal_core_genes = core_gene_utils.parse_personal_core_genes(species_name)
	non_shared_genes = core_gene_utils.parse_non_shared_reference_genes(species_name)
	core_genes = (non_shared_genes & personal_core_genes)
	core_gene_location_map = core_gene_utils.calculate_core_gene_location_map(species_name, core_genes)

	sys.stderr.write("Loading pickle from disk...\n")
	d = cPickle.load( open( pickle_filename,"rb") )
	sys.stderr.write("Done!\n")

	id_shortsnp_map = d['id_longsnp_map']
	allele_idx_map = d['allele_idx_map']
	B_matrix = d['B']
	S_matrix = d['S']
	
	gamete_idx_map = {('A','A'): (allele_idx_map['A'], allele_idx_map['A']),  ('A','R'): (allele_idx_map['A'], allele_idx_map['R']), ('R','A'): (allele_idx_map['R'], allele_idx_map['A']), ('R','R'): (allele_idx_map['R'], allele_idx_map['R'])}
	
	all_gametes = gamete_idx_map.values()
	on_diagonal_gametes = [gamete_idx_map[('A','A')], gamete_idx_map[('R','R')]]
	off_diagonal_gametes = [gamete_idx_map[('A','R')], gamete_idx_map[('R','A')]]
	
	twogamete_haplotypes = [on_diagonal_gametes, off_diagonal_gametes]

	# First loop through and find SNPs that themselves fail 2 gamete test 
	sys.stderr.write("Preprocessing SNVs...\n")
	maf_matrix = {}
	blacklisted_shortsnp_ids = set()
	total_snps = len(S_matrix)
	#corrected_logPstar = log(Pstar/total_snps)
	# We want to be conservative at ruling these out
	corrected_logPstar = log(inconsistency_pvalue)
	for focal_shortsnp_id in S_matrix:
	
		Bmin = B_matrix[focal_shortsnp_id].min()
		Btot = B_matrix[focal_shortsnp_id].sum()
		maf_matrix[focal_shortsnp_id] = Bmin*1.0/Btot
	
		self_S = S_matrix[focal_shortsnp_id][focal_shortsnp_id][0]
		self_S0 = S_matrix[focal_shortsnp_id][focal_shortsnp_id][0]
		
		S = long(self_S[off_diagonal_gametes].sum())
		S0 = self_S0[off_diagonal_gametes].sum()
		
		logP = calculate_poisson_logP(S,S0)
			
		if logP<corrected_logPstar:
			blacklisted_shortsnp_ids.add(focal_shortsnp_id)	 
	
		# Restrict to core genome
		focal_species_name, focal_contig, focal_location = id_shortsnp_map[focal_shortsnp_id]
		if not core_gene_utils.in_core_genome(core_gene_location_map, focal_contig, focal_location):
			blacklisted_shortsnp_ids.add(focal_shortsnp_id)
	
		# Restricted to tracked snps
		if id_shortsnp_map[focal_shortsnp_id] not in tracked_shortsnps:
			blacklisted_shortsnp_ids.add(focal_shortsnp_id)
		
	sys.stderr.write("Done!\n")
	
	focal_locations = set()
	for lmin,lmax in [(-1,5e04),(0,200),(200,2000),(2000,5e04)]:
		sys.stderr.write("Processing lmin=%d, lmax=%d...\n" % (lmin,lmax))
		sys.stderr.write("Finding threegamete testable pairs...\n")
		# Now go through and find pairs that are three gametes or higher
		testable_pairs = []
		for focal_shortsnp_id in S_matrix:
	
			if focal_shortsnp_id in blacklisted_shortsnp_ids:
				continue
	
			focal_species_name, focal_contig, focal_location = id_shortsnp_map[focal_shortsnp_id]
			
			focal_locations.add(focal_location)
		
			#print focal_contig, focal_location
			for target_shortsnp_id in S_matrix[focal_shortsnp_id]:
			
				if target_shortsnp_id in blacklisted_shortsnp_ids:
					continue
			
				
				target_species_name, target_contig, target_location = id_shortsnp_map[target_shortsnp_id] 
			
				#if focal_location==3779633:
				#	print focal_location, "->", target_location, S_matrix[focal_shortsnp_id][target_shortsnp_id][0].sum()
			    
				if target_shortsnp_id == focal_shortsnp_id:
					continue
			
				if target_contig!=focal_contig:
					continue
				
				if fabs(target_location-focal_location)>5e04:
					continue
				
				if fabs(target_location-focal_location)<lmin:
					#print "Should be here sometimes!"
					continue
				
				if fabs(target_location-focal_location)>lmax:
					#print "Should be here sometimes!"
					continue
				
				# restricted to long distances
				#if fabs(target_location-focal_location)<1e03:
				#	 continue	 
			
				Ss = S_matrix[focal_shortsnp_id][target_shortsnp_id][0]
				Stot = Ss.sum()
			
				focal_mac = Ss.sum(axis=1).min()
			
				if focal_mac < 3.5:
					continue
				
				if focal_mac < 0.1*Stot:
					continue
			
				target_mac = Ss.sum(axis=0).min()
				if target_mac < 3.5:
					continue
			
				if target_mac < 0.1*Stot:
					continue	   
			
				testable_pairs.append((focal_shortsnp_id, target_shortsnp_id))
		sys.stderr.write("Done! Found %d pairs!\n" % len(testable_pairs))
	
		#print sorted(focal_locations)
	
		if len(testable_pairs)==0:
			continue	
	
		sys.stderr.write("Calculating threegamete pvalues...\n")
		logPs = []
		for focal_shortsnp_id, target_shortsnp_id in testable_pairs:
					
			Ss = S_matrix[focal_shortsnp_id][target_shortsnp_id][0]
			S0s = S_matrix[focal_shortsnp_id][target_shortsnp_id][1]
			
			logP = max([calculate_poisson_logP(Ss[twogamete_haplotypes[i]].sum(),S0s[twogamete_haplotypes[i]].sum()) for i in xrange(0,len(twogamete_haplotypes))])
			logPs.append(logP)
		
		sys.stderr.write("Calculating qvalue threshold...\n")
		logPs = numpy.array(logPs)
		ntot = len(logPs)
		log_ntot = log(ntot)
		n_lower = 0
		for logP in sorted(logPs,reverse=True):
			n_lower += 1
			#n_lower =	(logPs <= logP).sum()
			if logP + log_ntot	< log(n_lower*complete_FDR):
				break
		
		# this is now the global threshold pvalue
		log_corrected_Pstar = logP
		sys.stderr.write("logp(FDR) = %g\n" % log_corrected_Pstar)

		
		sys.stderr.write("Calculating threegamete pairs...\n")
		# Now calculate which ones are three gamete
		total_pairs = len(testable_pairs)
		log_corrected_Pstar = log(2*Pstar/total_pairs/total_pairs)
		threegamete_pairs = []
		twogamete_pairs = []
		onegamete_pairs = []
		perfect_pairs = []
		extra_perfect_pairs = []
		for pair_idx in xrange(0,len(testable_pairs)):
	
			focal_shortsnp_id, target_shortsnp_id = testable_pairs[pair_idx]
			logP = logPs[pair_idx]
		
			if logP < log_corrected_Pstar:
				threegamete_pairs.append((focal_shortsnp_id, target_shortsnp_id))
			else:
				if logP < 0:
					# something consistent with perfect linkage
					perfect_pairs.append((focal_shortsnp_id, target_shortsnp_id))
				else:
					extra_perfect_pairs.append((focal_shortsnp_id, target_shortsnp_id))
					if (Ss>0).sum()==1:
						onegamete_pairs.append((focal_shortsnp_id, target_shortsnp_id))
					else:
						twogamete_pairs.append((focal_shortsnp_id, target_shortsnp_id))
		sys.stderr.write("Done! Found %d pairs!\n" % len(threegamete_pairs))
	
		sys.stderr.write("Calculating fourgamete testable pairs...\n")		  
		# Now get 4 gamete testable SNVS
		fourgamete_testable_pairs = []
		for focal_shortsnp_id, target_shortsnp_id in threegamete_pairs:
		
			f = maf_matrix[focal_shortsnp_id]
			g = maf_matrix[target_shortsnp_id]
		
			Stot = S_matrix[focal_shortsnp_id][target_shortsnp_id][0].sum()
		
			if Stot*f*g>=4:
				fourgamete_testable_pairs.append( (focal_shortsnp_id, target_shortsnp_id))
	
		fourgamete_testable_pairs = threegamete_pairs
		sys.stderr.write("Done! Found %d pairs!\n" % len(fourgamete_testable_pairs))

	
		sys.stderr.write("Calculating fourgamete pairs...\n")
		predictors = []
		responses = []
		fourgamete_pairs = []
		complete_pairs = []
		total_pairs = len(fourgamete_testable_pairs)
	
		if total_pairs>0:
	
			logPs = []
			Smins = []
			for focal_shortsnp_id, target_shortsnp_id in fourgamete_testable_pairs:
		
				focal_species_name, focal_contig, focal_location = id_shortsnp_map[focal_shortsnp_id]
				target_species_name, target_contig, target_location = id_shortsnp_map[target_shortsnp_id]
				distance = fabs(focal_location-target_location)
		
				predictors.append(log10(distance))
					
				# Testable! 
				Ss = S_matrix[focal_shortsnp_id][target_shortsnp_id][0]
				S0s = S_matrix[focal_shortsnp_id][target_shortsnp_id][1]
			
				logP = max([calculate_poisson_logP(Ss[all_gametes[i]],S0s[all_gametes[i]]) for i in xrange(0,len(all_gametes))])
				logPs.append(logP)
			
				Smins.append(Ss.min())
	
			logPs = numpy.array(logPs)
			Smins = numpy.array(Smins)
		
			sys.stderr.write("Calculating qvalue threshold...\n")
			ntot = len(logPs)
			log_ntot = log(ntot)
			n_lower = 0
			for logP in sorted(logPs,reverse=True):
				n_lower += 1
				if logP + log_ntot	< log(n_lower*fourgamete_FDR):
					break
	  
			# this is now the global threshold pvalue
			log_corrected_Pstar = logP
			sys.stderr.write("logp(FDR) = %g\n" % log_corrected_Pstar)
			# Bonferonni version
			#log_corrected_Pstar = log(2*Pstar/total_pairs/total_pairs)
			for pair_idx in xrange(0,len(fourgamete_testable_pairs)):
				focal_shortsnp_id, target_shortsnp_id = fourgamete_testable_pairs[pair_idx]
				logP = logPs[pair_idx]
				Smin = Smins[pair_idx]
			
				focal_species_name, focal_contig, focal_location = id_shortsnp_map[focal_shortsnp_id]
				target_species_name, target_contig, target_location = id_shortsnp_map[target_shortsnp_id]
				distance = fabs(focal_location-target_location)
		
				predictors.append(log10(distance))	   
			
				if logP < log_corrected_Pstar and (Smin>1.5):
					fourgamete_pairs.append((focal_shortsnp_id, target_shortsnp_id))		
					responses.append(1)
				else:
					responses.append(0)
					complete_pairs.append((focal_shortsnp_id, target_shortsnp_id))
			
		sys.stderr.write("Done! Found %d pairs!\n" % len(fourgamete_pairs))
	
	
		same_contig = 0
		diff_contig = 0
		within_50kb = 0
		involved_snvs = set([])
	
		output_items = [species_name, str(lmin), str(len(testable_pairs)), str(len(onegamete_pairs)), str(len(twogamete_pairs)), str(len(threegamete_pairs)), str(len(fourgamete_testable_pairs)), str(len(fourgamete_pairs))]
	
		for focal_shortsnp_id, target_shortsnp_id in fourgamete_pairs:
	
			focal_species_name, focal_contig, focal_location = id_shortsnp_map[focal_shortsnp_id]
			target_species_name, target_contig, target_location = id_shortsnp_map[target_shortsnp_id]
			 
			if target_contig==focal_contig:
			
				same_contig += 1
			
				if fabs(focal_location-target_location) < 5e04:
					within_50kb += 1
			else:
				diff_contig += 1
		
			involved_snvs.add(focal_shortsnp_id)
			involved_snvs.add(target_shortsnp_id)
		
			output_subitem = "%s|%d,%s|%d,%d,%d,%d,%d" % (focal_contig, focal_location, target_contig, target_location, S_matrix[focal_shortsnp_id][target_shortsnp_id][0][0,0], S_matrix[focal_shortsnp_id][target_shortsnp_id][0][1,0], S_matrix[focal_shortsnp_id][target_shortsnp_id][0][0,1], S_matrix[focal_shortsnp_id][target_shortsnp_id][0][1,1])
			output_items.append(output_subitem)
		
		output_file.write(" ".join(output_items))
		output_file.write("\n")		
					
	#print species_name, len(testable_pairs), len(onegamete_pairs), len(twogamete_pairs), len(threegamete_pairs), len(fourgamete_testable_pairs), len(fourgamete_pairs)
	
	#print diff_contig, same_contig, within_50kb
	
	#print len(involved_snvs)
output_file.close()	   