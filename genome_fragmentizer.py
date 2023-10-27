import os
import sys
import re
from Bio import SeqIO
import random
import numpy as np
import argparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
# from Bio.Alphabet import IUPAC

#Defining a function
def get_contigs(genome, genomes, completeness_cutoff, info_out):
	info_line = list()
	info_line.append(genome)
	genome_file = os.path.join(genomes+genome)
	output = genome_file.replace(".fna", "_"+str(completeness_cutoff)+"_fna.comp")
	output_name = genome.replace(".fna", "_"+str(completeness_cutoff)+"_fna.comp")
	info_line.append(output_name)
	info_line.append(str(completeness_cutoff))
		
	#Get sequence
	sequence = list() #This list will contain the full sequence
	for record in SeqIO.parse(genome_file, "fasta"):
		seq_id = record.description
		seq = list(record.seq)
		sequence.append(seq)
	sequence = [j for i in sequence for j in i] #Merging the sublists in case we have multiple entries in the genome.

	#Info about the genome
	genome_length = len(sequence) #Full genome length.
	final_length = int((completeness_cutoff*genome_length)/100) #Length of output genome based on completeness cutoff.
		
	#Getting random contigs length based on the “discrete uniform” distribution.
	data = np.random.randint(low=5000, high=50000, size=(1000), dtype='int')
	data = data.tolist()
		
	#Picking contigs until our full length is reached.
	contigs_list = list()
	while genome_length>=sum(contigs_list):
		random_value = random.choice(data)
		if random_value not in contigs_list: #All the contigs will have diff size
			contigs_list.append(random_value)

	#Checking that the picked contigs are not longer that the complete genome 
	if sum(contigs_list)>genome_length:
		contigs_list.pop() #If so we will delete the last contig 
		diff = genome_length-sum(contigs_list) #Adjusting to reach the genome length
		contigs_list.append(diff)

	info_line.append(str(len(contigs_list)))
	info_line.append(str(round(len(contigs_list*completeness_cutoff)/100)))
	info_line.append(str(genome_length))
	info_line.append(str(final_length))

	#Fragment the complete genome into the picked contigs
	contigs_seqs = list() #Fragmented genome
	while len(sequence) != 0:
		for contig in contigs_list:
			seq_conting = sequence[0:contig] #Picking one contig at the time
			del sequence[0:contig] #Deleting the picked contig from the sequence
			contigs_seqs.append(seq_conting) #Appending the contig
		
	#Create a list of numbers equal to the number of contigs in the complete genome
	list_num = []
	for i in range(0,len(contigs_list)):
		list_num.append(i)		

	#If the last contig is shorter than 5000 it's going to be deleted from the fragmented genome
	if len(contigs_seqs[-1]) < 5000:
		contigs_seqs.pop()
		list_num.pop()
		
	#Picking random numbers, which are the contigs that will be in the output genome
	contigs_needed = int(round(completeness_cutoff*len(contigs_list))/100) #Contigs needed
	contigs_to_pick = random.sample(list_num, k=contigs_needed) #Random contigs
	contigs_to_pick.sort() #Sorting contigs
 		
	total = 0
	final_records = list()
	contig_num = 0
	#Get the random contigs from the original sequence
	for i in contigs_to_pick:
		random_contig = contigs_seqs[i]
		seq_ready = "".join(random_contig)
		contig_num = contig_num+1
		total = total+len(seq_ready)
		contig_len=str(len(seq_ready))
		contig_id = "contig_"+str(contig_num)+"_"+contig_len+"_bp"
		newseq = Seq(seq_ready)
		newrecord = SeqRecord(newseq, id=contig_id, description = "")
		final_records.append(newrecord)

	#Get final sequences file
	SeqIO.write(final_records, output, "fasta")

	info_line.append(str(total))
	info_line.append(str(int(final_length-total)))
	real_comp = (total*100)/genome_length
	info_line.append(str(real_comp))
	info_line = "\t".join(info_line)
	info_out.write(info_line+"\n")

	print(info_line)
###

def run_program(genomes, completeness_cutoff, info_out):
	#Output info file#
	header = ["Input_file", "Outfile", "Completeness", "N_contigs", "Incomplete_contigs", "100_length", "New_length", "Output_length", "Diff", "Real_comp"]
	header = "\t".join(header)
	info_out.write(header+"\n")
	print(header)
	for genome in os.listdir(genomes):
		if genome.endswith(".fna"):
			get_contigs(genome, genomes, completeness_cutoff, info_out)
	info_out.close()

### Input/Output files ###
def main(argv=None):
	args_parser = argparse.ArgumentParser(description="Script that creates incomplete genomes based on a given cutoff.\nComplete genomes must be used.\nOutput files will be stored in the folder where the input files are.\nContigs will be 5000-50000 bp.")
	args_parser.add_argument('-i', '--input_folder', required=True, help='Input folder where the FNA sequence files are.')
	args_parser.add_argument('-c', '--cutoff', required=True, help='Completeness cutoff needed for the output genomes.')
	args_parser.add_argument('-o', '--output_info', required=True, help='Output file that will have information about the incomplete genomes (i.e., Output name, completeness picked, genome size, etc.')
	args_parser = args_parser.parse_args()

	genomes = args_parser.input_folder
	completeness_cutoff = int(args_parser.cutoff)
	info_out = open(args_parser.output_info, "w")

	run_program(genomes, completeness_cutoff, info_out)

if __name__ == '__main__':
	status = main()
	sys.exit(status)
#######
