#!/usr/bin/env python

import argparse, os, subprocess as sp
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pyfaidx import Fasta, wrap_sequence

'''This function is to extract CDS from the Genbank file'''
def extractCDS(input_path, output_path):
	
	input_pathname = os.path.abspath(input_path)
	with open(input_pathname) as input_handle:
	
		if not os.path.isdir(output_path):
			os.makedirs(output_path)
		
		genelist =[]
		for seq_record in SeqIO.parse(input_handle, "genbank"):
			gb_feature = seq_record.features
			for seq_feature in seq_record.features:
				if seq_feature.type == "CDS":
					start = seq_feature.location.nofuzzy_start
					end = seq_feature.location.nofuzzy_end
					newseq = seq_record.seq[start:end]
					
					if 'gene' in seq_feature.qualifiers: 
						if len(seq_feature.qualifiers['gene']):
							genename = seq_feature.qualifiers['gene'][0]
						else:
							genename = "No Name"
					else:
						genename = "No Name"

					gene_file_name = os.path.basename(genename) + ".fa"
					gene_file_path = os.path.join(output_path, gene_file_name)
					
					new_seq = SeqRecord(newseq, 
									id=genename)
					
					with open(gene_file_path, "w+") as output_handle:
						if genename not in genelist:
							genelist.append(genename)
							SeqIO.write(new_seq, output_handle, 'fasta')


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-input_path",
	                    help="input directory path",
	                    required=True)

	parser.add_argument("-output_path",
	                    help="output directory path" ,
	                    required=True)

	args = parser.parse_args()
	input_path = args.input_path
	output_path = args.output_path

	extractCDS(input_path, output_path)
	
if __name__ == "__main__":
	main()

