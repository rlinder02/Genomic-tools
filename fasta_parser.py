#!/usr/bin/env python 
# 4-space indented, v1.0.0
# File name: fasta_parser.py
# Description: This script can be used to calculate various metrics from an input fasta file.
# Author: Robert Linder
# Date: 2022-09-23

import argparse
import pandas as pd
import re
from functools import wraps
from seq_classes import Fasta



def parse_args():
	"""this enables users to provide input as command line arguments to minimize the need for directly modifying the script; please enter '-h' to see the full list of options"""
	parser = argparse.ArgumentParser(description="List a variety of metrics from a FASTA file")
	parser.add_argument("FASTA_file", type=str, help="FASTA file to retrieve information from")
	parser.add_argument('-r', '--reading-frame', action='store', type=int, default=1, choices=range(1,6), help="Use this reading frame for reporting ORFs")
	parser.add_argument('-s', '--substring-length', type=int, default=3, help="Input this substring length to find the most frequently occurring substring of this length")
	parser.add_argument("--seq-id", type=str, help="Retrieve information about this particular sequence; ensure the sequence name is enclosed in quotes")
	parser.add_argument('-m', "--metric", type=str, choices=['max', 'min', 'median', 'mean'], help="Retrieve a metric of interest about this particular sequence")
	args = parser.parse_args()
	return args

def unpack_dictionaries(fn):
	@wraps(fn)
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
	fasta_file = Fasta.preprocess_fasta(inputs.FASTA_file)
	fasta_df = fasta_file.create_fasta_df()
	## If you are interested in finding ORFs in a reading frame other than the first, please specify using the -r flag
	all_orfs = fasta_file.find_orfs(inputs.reading_frame)
	## If you are interested in finding substrings of a length other than 3, please specify using the -s flag
	all_repeats = fasta_file.find_substrings(inputs.substring_length)
	orf_df = create_orf_df(*all_orfs)
	if inputs.metric and not orf_df.empty:
		orf_info = seq_info(orf_df, eval(inputs.metric))
	if inputs.seq_id:
		sequence_info = seq_retrieval(fasta_df, seq_id = inputs.seq_id)
	if inputs.metric and any(all_repeats):
		repeat_info = repeat_profiler(all_repeats, eval(inputs.metric))


if __name__ == "__main__":
	main()