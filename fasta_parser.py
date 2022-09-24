#!/usr/bin/env python 
# 4-space indented, v1.1
# File name: fasta_parser.py
# Description: This script can be used to calculate various metrics from an input fasta file.
# Author: Robert Linder
# Date: 2022-09-23

import re
import pandas as pd
from seq_classes import Fasta
import argparse


def parse_args():

	parser = argparse.ArgumentParser(description="List a variety of metrics from a FASTA file")
	parser.add_argument("FASTA_file", type=str, help="FASTA file to retrieve information from")
	parser.add_argument("-rframe", type=int, default=1, help="Use this reading frame for reporting ORFs")
	parser.add_argument("-sublen", type=int, default=3, help="Input this substring length to find the most frequently occurring substring of this length")
	parser.add_argument("--seqid", type=str, default=False, help="Retrieve information about this particular sequence")
	args = parser.parse_args()
	return args

def unpack_dictionaries(fn):
	def wrapper(*args):
		counter = 0
		output = pd.DataFrame()
		for dictionary in args:
			counter += 1
			df_dictionary = fn(**dictionary)
			output = pd.concat([output, df_dictionary], ignore_index=True)
			if counter == len(args):
				return output
	return wrapper

@unpack_dictionaries
def create_orf_df(**kwargs):
	"""create a dataframe from a list of dictionaries (that must be unpacked using *) containing any number of key: value pairs, with each key representing a different parameter"""
	orf_df = pd.DataFrame([kwargs])
	orf_df['input'] = 'orf'
	return orf_df

def seq_info(seq_df, fun):
	"""Find the lengths and identifiers of sequences of interest from a dataframe of sequences"""
	find_info = seq_df[seq_df['Length(bp)'] == fun(seq_df['Length(bp)'])].copy()
	find_info['Type'] = fun.__name__
	print(f"This is the {fun.__name__} value from the {seq_df['input'].iloc[0]} input:")
	print(find_info)

def seq_retrieval(seq_df, seq_id):
	"""Retrieve information about a seqeunce of interest by its' identifier from a dataframe of sequences"""
	find_info = seq_df[seq_df['Identifier'] == seq_id]
	print(f"This has information about the {seq_id} sequence from the {seq_df['input'].iloc[0]} input:")
	print(find_info)

def repeat_profiler(substring_dict, fun):
	"""Retrieve information about a substring of interest by a user-defined function"""
	print(f"This is the {fun.__name__} repeat value:")
	print([(k, v) for k,v in substring_dict.items() if v == fun(substring_dict.values())])


def main():
	inputs = parse_args()
	fasta_file = Fasta.preprocess_fasta(inputs)
	fasta_df = fasta_file.create_fasta_df()
	all_orfs = fasta_file.find_orfs(read_frame = 3)
	all_repeats = fasta_file.find_substrings(substring_len = 5)
	orf_df = create_orf_df(*all_orfs)
	max_info = seq_info(orf_df, max)
	if args.seqid:
		sequence_info = seq_retrieval(fasta_df, seq_id = ">gi|142022655|gb|EQ086233.1|91 marine metagenome JCVI_SCAF_1096627390048 genomic scaffold, whole genome shotgun sequence")
	max_repeat_count = repeat_profiler(all_repeats, max)

if __name__ == "__main__":
	main()