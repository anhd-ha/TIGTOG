#!/usr/bin/env python

import sys, os, re, shlex, subprocess, itertools, argparse, time, glob
from collections import defaultdict
from Bio import SeqIO
from Bio import SeqUtils
from Bio.Seq import Seq
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
import joblib
import csv


# predict proteins from genome FNA file
def predict_proteins(genome_file, project, redo):

	seqdict = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
	if len(seqdict) < 1:
		raise (genome_file+" does not appear to be in FASTA format!")

	file_name = os.path.basename(genome_file)
	file_base = os.path.splitext(file_name)[0]

	nucl_length = float(0)
	for i in SeqIO.parse(genome_file, "fasta"):
			nucl_length += len(i.seq)

	if nucl_length > 5000000:
		print("Warning: "+file_name+" is a very large genome. Are you sure this is a viral genome or MAG?")
		large = True
	else:
		large = False
	
	path_base = os.path.splitext(project)[0]
	protein_file = os.path.join(path_base, file_base+".faa")
	
	cmd = "prodigal -p meta -i "+ genome_file +" -a "+ protein_file
	cmd2 = shlex.split(cmd)
	if not redo:
		subprocess.call(cmd2, stdout=open("out.txt", "w"), stderr=open("err.txt", "w"))


	prot_length = float(0)
	for j in SeqIO.parse(protein_file, "fasta"):
		prot_length += len(j.seq)

	density = 100 * (prot_length * 3) / nucl_length
	if density < 85:
		print("Warning: "+file_name+" has low coding density (<85%). Results may not be accurate.")
		low_dens = True
	else:
		low_dens = False


	return protein_file, file_base, large, nucl_length, low_dens, density


# run HMMER3
def run_hmmer(input_file, cpus, redo, evalue):

	output_file = re.sub(".faa", ".gvog.out", input_file)
	cmd = "hmmsearch --cpu "+ cpus +" --tblout "+ output_file +" db/merged_GVOGs.hmm "+ input_file
	cmd2 = shlex.split(cmd)

	if not redo:
		subprocess.call(cmd2, stdout=open("out.txt", "w"), stderr=open("err.txt", "w"))
		os.remove("out.txt")
	return output_file


# define function for parsing HMMER3 output
def parse_hmmout(hmmout, evalue):

	input = open(hmmout, "r")
	hit_dict = defaultdict(lambda:"NA")
	bit_dict = defaultdict(float)

	for i in input.readlines():
		line = i.rstrip()
		if line.startswith("#"):
			pass
		else:
			newline = re.sub("\s+", "\t", line)
			tabs = newline.split("\t")
			hit = tabs[2]
			protein = tabs[0]
			eval = float(tabs[4])
			score = float(tabs[5])

			if score > bit_dict[protein] and eval <= float(evalue):
				bit_dict[protein] = score
				hit_dict[protein] = hit
			else:
				pass

	list_g = list()
	g_names = open("db/names_imp_GVOGs_both_levels.csv","r")
	for line in g_names.readlines(): 
		g = line.rstrip()
		g2 = re.sub("hmm$", "trim", g)
		list_g.append(g2)

	final_gvog_dict = {gvg: 0 for gvg in list_g}  	
	for key, value in hit_dict.items():
		final_gvog_dict[value] +=1

	hmmout = os.path.split(hmmout)[1]  
	seq_name = re.sub(".gvog.out","",hmmout) 
	df = pd.DataFrame(final_gvog_dict, index = [seq_name])  
	
	count_tab = "GVOG_count_unseen.tsv"
	df.to_csv(count_tab, sep ="\t",index_label="Sequence", mode ="a", header=not os.path.exists(count_tab)) 

	return df


# add GC content of sequences
def get_gc(folder):

	file_list = os.listdir(folder)
	file_names = []
	gc_content = []

	for genome in file_list:
		file_name = os.path.splitext(genome)[0]
		file_names.append(file_name)

		genome_file = os.path.join(folder, genome)
		genome_seq = open(genome_file, "r")
		full_seq = Seq("")
		for record in SeqIO.parse(genome_seq, "fasta"):
			seq = record.seq.strip().strip('*')
			full_seq = full_seq + seq
		gc = SeqUtils.GC(full_seq)
		gc = round(gc, 2)
		gc_content.append(gc)

	dat = pd.DataFrame()
	dat["Sequence"] = file_names
	dat["GC_content"] = gc_content

	return dat


# predict taxonomic classification
def tax_predict(output_f):
	# Predict orders
	print("Generating taxonomy predictions")
	clf_ord = joblib.load("clf/clf_Order_final.joblib")
	d = pd.read_csv('feat_order.tsv',sep = '\t')
	ord_test = d[d.columns[1:152]]
	name = d['Sequence']
	ord_pred = clf_ord.predict(ord_test)

	# Confidence of the predictions
	confidences_predictions = clf_ord.predict_proba(ord_test)
	class_labels = clf_ord.classes_
	confidence = []
	for i, prediction in enumerate(confidences_predictions):
	    max_index = np.argmax(prediction)
	    confidence.append(prediction[max_index])

	# Predict families  
	clf_fam = joblib.load("clf/clf_Fam_final.joblib")
	fd = pd.read_csv('feat_fam.tsv',sep = '\t')
	fam_test = fd[fd.columns[1:209]]
	fname = fd['Sequence']
	fam_pred = clf_fam.predict(fam_test)

	# Confidence of the predictions
	confidences_predictions_fam = clf_fam.predict_proba(fam_test)
	class_labels_fam = clf_fam.classes_
	confidence_fam = []
	for i, prediction_fam in enumerate(confidences_predictions_fam):
	    max_index = np.argmax(prediction_fam)
	    confidence_fam.append(prediction_fam[max_index])

	# Print out predictions and confidence
	pred = pd.DataFrame({'Sequence':name,'Predicted_Order':ord_pred, 'Order_Predict_Probability': confidence,'Predicted_Family':fam_pred, 'Family_Predict_Probability': confidence_fam})

	# clean up feature matrices
	fp_ord = os.path.join(output_f, 'feat_order.tsv')
	fp_fam = os.path.join(output_f, 'feat_fam.tsv')
	os.rename('feat_order.tsv', fp_ord)
	os.rename('feat_fam.tsv', fp_fam)

	return pred


# run lastp for AAI calculation
def aai_lastp(input_folder, output_folder, cpus):

	length_dict = {}
	proteins_in_genome = defaultdict(list)

	# define lastout parsing function
	def parse_lastout(lastout):
			bit_dict = defaultdict(float)
			hit_dict = defaultdict(list)
			perc_dict = defaultdict(float)
			done = {}		
			handle = open(lastout, "r")
			for j in handle.readlines():
				if j.startswith("#"):
					pass
				else:
					line = j.rstrip()
					tabs = line.split("\t")
					query = tabs[0]
					hit = tabs[1]
					perc = float(tabs[2])
					evalue = float(tabs[10])
					bit = float(tabs[11])
					aln = float(tabs[3])

					done[query] = query
					if evalue < 1e-5:
						if bit >= bit_dict[query]:
							hit_dict[query].append(hit)
							bit_dict[query] = bit
							perc_dict[query] = perc
						else:
							pass						
			return hit_dict, bit_dict, perc_dict

	input_db = "db/AAI_db/"
	folders = os.listdir(input_db)
	for faa in folders:
		if faa.endswith(".faa"):
			prefix = re.sub(".faa", "", faa)
			db = os.path.join(input_db, prefix)
			filepath = os.path.join(input_db, faa)
			lastpath = os.path.join(input_db, prefix+".suf")
			cmd = "lastdb -p "+ db +" "+ filepath
			cmd2 = shlex.split(cmd)
			subprocess.call(cmd2, stdin=open("stdout.txt", "w"), stdout=open("stderr.txt", "w"))

	AAI_df = pd.DataFrame(columns=['MAG', 'Ref', 'AAI', 'num_hits', 'num_prots', 'AF'])

	for fna in os.listdir(input_folder):
		print("Running AAI calculation for " + fna)
		file_name = os.path.basename(fna)
		mag_prefix = os.path.splitext(file_name)[0]		
		path_base = os.path.splitext(output_folder)[0]
		faa2 = os.path.join(path_base, mag_prefix, mag_prefix +".faa")

		# get length of each protein seq
		for prot in SeqIO.parse(faa2, "fasta"):
			length_dict[prot.id] = len(prot.seq)
			proteins_in_genome[mag_prefix].append(prot.id)

		for faa in folders:
			if faa.endswith(".faa"):
				prefix = re.sub(".faa", "", faa)
				db = os.path.join("db/AAI_db/", prefix)

				output = os.path.join(output_folder,mag_prefix +"__"+ prefix +".lastout")
				cmd = "lastal -P " + cpus + " -m 500 -f BlastTab "+ db +" "+ faa2
				cmd2 = shlex.split(cmd)
				subprocess.call(cmd2, stdin=open("stderr.txt", "w"), stdout=open(output, "w"))

				hit_dict, bit_dict, perc_dict = parse_lastout(output)

				aai_dict = defaultdict(list)
				for query in hit_dict:
					hit_list = hit_dict[query]
					aai_dict[mag_prefix +"__"+ prefix].append(perc_dict[query])

				if len(aai_dict[mag_prefix +"__"+ prefix]) == 0:
					aai = 0
				else:
					aai = np.mean(aai_dict[mag_prefix +"__"+ prefix])

					num_hits = float(len(aai_dict[mag_prefix +"__"+ prefix]))
					num_prots = float(len(proteins_in_genome[mag_prefix]))
					af = 100*(num_hits / num_prots)

					new_aai = pd.DataFrame([{'MAG': mag_prefix, 'Ref': prefix, 'AAI': aai, 'num_hits': num_hits, 'num_prots': num_prots,'AF': af}])
					AAI_df = pd.concat([AAI_df, new_aai], ignore_index=True)

	AAI_df.to_csv("AAI.tsv", sep="\t", index=False)
	to_del = glob.glob(os.path.join(output_folder,'*.lastout'))
	for file_path in to_del:
		os.remove(file_path)	
	return AAI_df

# parse for best AAI hits
def parse_best_aai_hit(output_folder):

	aai = open("AAI.tsv", "r")
	hit_dict = defaultdict(lambda:"none")
	aai_dict = defaultdict(float)
	af_dict = defaultdict(float)
	for i in aai.readlines():
		line = i.rstrip()
		if line.startswith("MAG"):
			pass
		else:
			newline = re.sub("\s+", "\t", line)
			tabs = newline.split("\t")
			hit = tabs[1]
			mag = tabs[0]
			aai = float(tabs[2])
			af = float(tabs[5])
			if aai > aai_dict[mag] and af >= 20:
				aai_dict[mag] = aai
				hit_dict[mag] = hit
				af_dict[mag] = af
			else:
				pass
	best_aai = []
	for key, value in aai_dict.items():
		value2 = hit_dict[key]
		value3 = af_dict[key]
		best_aai.append((key, value, value2, value3))

	aai_hit = pd.DataFrame(best_aai,columns=['Sequence','AAI','Best_AAI_hit','AF'])
	aai_tax = pd.read_csv('db/GVDB_tax.csv')
	aai_fin = aai_hit.set_index('Best_AAI_hit').join(aai_tax.set_index('GV'))	
	aai_fin['Best_AAI_hit'] = aai_fin.index
	aai_fin = aai_fin[['Best_AAI_hit','Sequence','AAI','AF','Order','Family','Genus']]
	aai_fin.columns = ['Best_AAI_hit', 'Sequence', 'AAI','Alignment_fraction (%)','Best_hit_Order','Best_hit_Family', 'Best_hit_Genus']
	
	
	fp = os.path.join(output_folder, 'AAI.tsv')
	os.rename('AAI.tsv', fp)

	return aai_fin


####### main function that runs the program ######
def run_program(input, project, evalue, cpus, aai, redo):

	print("Processing sequences")
	file_list = os.listdir(input)
	process_info =[]
	for i in file_list:
		genome_folder = os.path.join(project, i)
		genome = os.path.join(input, i)

		# create output directories
		foldername = os.path.splitext(genome_folder)[0]
		if os.path.isdir(foldername):
			pass
		else:
			os.mkdir(foldername)

		relpath = os.path.split(project)[1]
		relpathbase = os.path.splitext(relpath)[0]
		base = os.path.splitext(project)[0]


		# predict proteins and evaluate sequence
		protein_file, name, too_large, size, low_density, density = predict_proteins(genome, genome_folder, redo)
		genome_info = {'Sequence': name, 'Warning: Large size':too_large, 'Genome_Size (bp)':size, 'Warning: Low coding density': low_density, 'Coding_density (%)': density}
		process_info.append(genome_info)
		
		# run HMMER3 searches, and parse outputs
		print("Running feature search on "+ i)
		vog_out  = run_hmmer(protein_file, cpus, redo, evalue)
		vog_hit = parse_hmmout(vog_out, evalue)

	seq_processing = pd.DataFrame(process_info)

	# print feature matrices
	table = pd.read_table("GVOG_count_unseen.tsv")
	GC = get_gc(input)		
	table2 = GC.set_index('Sequence').join(table.set_index('Sequence'))	
	
	ord_feat_list = pd.read_csv('db/names_imp_gvog_order.csv', header=None)[0].tolist()
	selected_feat_ord = table2[ord_feat_list]
	selected_feat_ord.to_csv("feat_order.tsv", sep ="\t")

	fam_feat_list = pd.read_csv('db/names_imp_gvog_fam.csv', header=None)[0].tolist()
	selected_feat_fam = table2[fam_feat_list]
	selected_feat_fam.to_csv("feat_fam.tsv", sep ="\t")

	os.remove("GVOG_count_unseen.tsv")


	# predict taxonomy
	taxo = tax_predict(project)
	taxo2 = seq_processing.set_index('Sequence').join(GC.set_index('Sequence'))	
	tax = taxo2.join(taxo.set_index('Sequence'))	

	# calculate AAI if requested:
	if aai:
		aai_calc = aai_lastp(input, project, cpus)
		aai_final = parse_best_aai_hit(project)

		# print result table
		result = tax.join(aai_final.set_index('Sequence'))	
		result.to_csv("prediction_result.tsv", sep ="\t")

	else:
		tax.to_csv("prediction_result.tsv", sep ="\t")
		
	fin_result = pd.read_csv("prediction_result.tsv", sep ="\t")
	warnings = ['Warning: Large size','Warning: Low coding density']
	remaining_cols = [col for col in fin_result.columns if col not in warnings]
	ordered_cols = remaining_cols + warnings
	fin_result = fin_result[ordered_cols]
	fin_result.to_csv("pred_result.tsv", sep='\t', index=False)

	print("done!")
	
	# clean up
	os.remove("err.txt")
	os.remove("prediction_result.tsv")
	os.rename("pred_result.tsv", project+"."+"prediction_result.tsv")


########################################################################
# #### use argparse to run through the command line options given #######
# #######################################################################
def main(argv=None):

	args_parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="TIGTOG: Taxonomic Information of Giant viruses using Trademark Orthologous Groups", epilog='*******************************************************************\n\n*******************************************************************')
	args_parser.add_argument('-i', '--input', required=True, help='Folder containing input FASTA file (ending in .fna)')
	args_parser.add_argument('-n', '--project', required=True, help='project name for outputs')
	args_parser.add_argument('-e', '--evalue', required=False, default=str(1e-10), help='e-value that is passed to HMMER3 for the GVOG hmmsearch (default=1e-10)')
	args_parser.add_argument('-t', '--cpus', required=False, default=str(1), help='number of cpus to use for the HMMER3 search')
	args_parser.add_argument('-a', '--aai', type=bool, default=False, const=True, nargs='?', help='run AAI calculation')
	args_parser.add_argument('-r', '--redo', type=bool, default=False, const=True, nargs='?', help='run without re-launching prodigal and HMMER3 (for quickly re-calculating outputs with different parameters if you have already run once)')
	args_parser = args_parser.parse_args()

	# set up object names for input/output folders
	input = args_parser.input
	project = args_parser.project
	evalue = str(args_parser.evalue)
	cpus = args_parser.cpus
	aai = args_parser.aai
	redo = args_parser.redo

	project = project.rstrip("/")
	if os.path.isdir(project):
		pass
	else:
		os.mkdir(project)
		
	run_program(input, project, evalue, cpus, aai, redo)

	return 0

if __name__ == '__main__':
	status = main()
	sys.exit(status)





