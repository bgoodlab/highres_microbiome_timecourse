import numpy
import parse_midas_data

filename = "selection_test_output_072419.txt"
file = open(filename,"r")

line = file.readline() # blank
line = file.readline()

data = {}

output_filename = parse_midas_data.analysis_directory+"table_S3.csv"

output_file = open(output_filename,"w")

output_items = ["SNV Cluster", "Ne all", "Ne post", "Ne pre", "Ne pre+post", "Pvalue all post", "Pvalue all pre", "Pvalue all both", "Pvalue joint post", "Pvalue joint pre", "Pvalue joint both", "Pvalue separate post", "Pvalue separate pre", "Pvalue separate both"]	
output_file.write(", ".join([str(item) for item in output_items]))
output_file.write("\n")
	

while line!="":
	
	species_name = line.strip()
	while not line.startswith('Effective generations per day:'):
		line = file.readline()
		
	gensperday = float(line.split(":")[1])
	
	file.readline() # Nhat
	
	while not line.startswith('All Ne:'):
		line = file.readline()
		
	all_Ne = float(line.split(":")[1])
	line = file.readline() # post Ne
	post_Ne = float(line.split(":")[1])
	line = file.readline() # pre Ne
	pre_Ne = float(line.split(":")[1])
	line = file.readline() # postpre Ne
	postpre_Ne = float(line.split(":")[1])
	
	file.readline() # "All"
	line = file.readline()
	pvalue_all_post = float(line.split(":")[1])
	line = file.readline()
	pvalue_all_pre = float(line.split(":")[1])
	line = file.readline()
	pvalue_all_both = float(line.split(":")[1])
	
	file.readline() # "joint"
	line = file.readline()
	pvalue_joint_post = float(line.split(":")[1])
	line = file.readline()
	pvalue_joint_pre = float(line.split(":")[1])
	line = file.readline()
	pvalue_joint_both = float(line.split(":")[1])
	
	file.readline() # "separate"
	line = file.readline()
	pvalue_separate_post = float(line.split(":")[1])
	line = file.readline()
	pvalue_separate_pre = float(line.split(":")[1])
	line = file.readline()
	pvalue_separate_both = float(line.split(":")[1])
	
	line = file.readline() # blank
	line = file.readline() # new species
	
	output_items = [species_name, all_Ne/gensperday, post_Ne/gensperday, pre_Ne/gensperday, postpre_Ne/gensperday, pvalue_all_post, pvalue_all_pre, pvalue_all_both, pvalue_joint_post, pvalue_joint_pre, pvalue_joint_both, pvalue_separate_post, pvalue_separate_pre, pvalue_separate_both]
	
	output_file.write( ", ".join([str(item) for item in output_items]) )
	output_file.write("\n")
	
output_file.close()
file.close()
	 