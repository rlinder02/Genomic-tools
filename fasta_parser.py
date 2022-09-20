#!/usr/bin/env python3, 4-space indented, v1.0
# File name: fasta.parser.py
# Description: 
# Author: Robert Linder
# Date: 2022-07-31

from Bio import SeqIO as Seq
import pandas as pd


def create_fasta_df(fasta):
	"""Create a dataframe of the names of each fasta record and the corresponding lengths of each sequence"""
	header_names = list({name for name in fasta.keys()})
	record_seqs = list({len(record) for record in fasta.values()})
	combine_lists = list(zip(header_names, record_seqs))
	fasta_df = pd.DataFrame(combine_lists, columns = ['Identifiers', 'Sequence_lengths'])
	return fasta_df

def sequence_info(fasta_df):
	"""Find the values and identifiers of the longest and shortest sequences"""
	max_info = fasta_df[fasta_df['Sequence_lengths'] == fasta_df['Sequence_lengths'].max()]
	max_info['Type'] = 'Max'
	min_info= fasta_df[fasta_df['Sequence_lengths'] == fasta_df['Sequence_lengths'].min()]
	min_info['Type'] = "Min"
	all_info = pd.concat([max_info, min_info])
	return(all_info)

def find_orfs(fasta, frame_start):
	"""Find all ORFS in all sequences in a fasta file"""
	orf_list = []
	for name, value in fasta.items():
		frame = value.seq[frame_start-1:].upper()
		for i in range(0, len(frame), 3):
			codon = frame[i:(i+3)]
			if codon == start:
				starting = i
				for j in range(starting, len(frame), 3):
					codon2 = frame[j:(j+3)]
					if codon2 in stop:
						ending = j
						sequence = frame[starting:(ending+3)]
						orf_list.append([name, sequence, starting+frame_start])
						break
	#print(max([len(orf[1]) for orf in orf_list]))
	return(orf_list)

def orf_info(orf_list, fasta):
	orf_names = [orf[0] for orf in orf_list]
	seqs = [orf[1] for orf in orf_list]
	orf_starts = [orf[2] for orf in orf_list]
	orf_lens = [len(orf[1]) for orf in orf_list]
	combined_orf_info = list(zip(orf_names, seqs, orf_starts, orf_lens))
	orf_df = pd.DataFrame(combined_orf_info, columns = ['Identifier', 'Sequence', 'Start_idx', 'Orf_length'])
	return(orf_df)


def find_max(orf_df, seq_name):
	max_orf = orf_df.loc[orf_df['Orf_length'] == max(orf_df['Orf_length'])]
	#print(max_orf)
	an_orf = orf_df.loc[orf_df['Identifier'] == seq_name]
	#print(an_orf)
	max_an_orf = an_orf.loc[an_orf['Orf_length'] == max(an_orf['Orf_length'])]
	#print(max_an_orf)
	return(max_an_orf)

def find_repeats(fasta, n):
	repeat_list = []
	for value in fasta.values():
		sequence = value.seq[:].upper()
		for i in range(0, len(sequence)):
			if len(sequence[i:(i+n)]) == n:
				repeat_list.append(sequence[i:(i+n)])
	unique_repeats = set(repeat_list)
	repeat_dict = {x: repeat_list.count(x) for x in unique_repeats}
	return(repeat_dict)


def main(frame_start, seq_name, repeat_len, fasta):
	fasta_df = create_fasta_df(fasta)
	seq_info = sequence_info(fasta_df)
	all_orfs = find_orfs(fasta_seq, frame_start)
	orf_df = orf_info(all_orfs, fasta)
	max_info = find_max(orf_df, seq_name)
	finding_repeats = find_repeats(fasta, repeat_len)
	mx = max(finding_repeats.values())
	print([k for k,v in finding_repeats.items() if v == mx])


if __name__ == "__main__":
	fasta_seq = Seq.index("dna2.fasta", "fasta")
	start = 'ATG'
	stop = ['TAA', 'TAG', 'TGA']
	frame_start = 3
	seq_name = "gi|142022655|gb|EQ086233.1|16"
	repeat_len = 7
	main(frame_start, seq_name, repeat_len, fasta_seq)