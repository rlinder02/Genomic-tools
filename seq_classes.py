#!/usr/bin/env python 
# 4-space indented, v1.0
# File name: seq_classes.py
# Description: This contains a list of types of sequences commonly encountered during genomics analyses and creates  a distinct class for each. Nucleic acids is the parent class that all other sequence classes inherit from.
# Author: Robert Linder
# Date: 2022-09-20

import re


class NucleicAcid:
	"""define a nucliec acid, including unnatural sequences that contain Ts and Us; can be passed in as a single string or a list of strings; the output will be a list"""

	_nuc = re.compile('[ATCGUN]+')
	seq_counter = 0

	def __init__(self, sequence):
		if not isinstance(sequence, list) and not isinstance(sequence, str):
			raise TypeError("Please enter a list of sequence strings or a single sequence string")
		if isinstance(sequence, list):
			self.sequence = [seq.upper().replace(" ", "").replace("-", "") for seq in sequence]
		else:
			self.sequence = [sequence.upper().replace(" ", "").replace("-", "")]
		for seq in self.sequence:
			NucleicAcid.seq_counter += 1
			if not NucleicAcid._nuc.fullmatch(seq):
				if len(self.sequence) > 1:
					raise ValueError(f"Sequence {NucleicAcid.seq_counter} is not a valid nucleic acid sequence")
				else:
					raise ValueError("Please use a valid nucleic acid sequence")

	def __repr__(self):
		if len(self.sequence) > 1:
			return f"This is a list of {NucleicAcid.seq_counter} nucleic acid sequences"
		else:
			return f"This is a nucleic acid sequence of length {len(self.sequence[0])}bp"

	def reverse(self):
		"""reverses a nucleic acid sequence, agnostic of type"""
		return [seq[::-1] for seq in self.sequence]
	

class DNA(NucleicAcid):
	"""define DNA, how the bases are paired, and how DNA can be converted into RNA"""
	
	_dna = re.compile('[ATCGN]+')
	_dna_pairings = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N': 'N'}
	_dna_to_rna = {'A': 'A', 'C': 'C', 'T': 'U', 'G': 'G', 'N': 'N'}
	dna_counter = 0

	def __init__(self, sequence):
		super().__init__(sequence)
		for seq in self.sequence:
			DNA.dna_counter += 1
			if not DNA._dna.fullmatch(seq):
				if len(self.sequence) > 1:
					raise ValueError(f"Sequence {DNA.dna_counter} is not a valid DNA sequence")
				else:
					raise ValueError("Please use a valid DNA sequence")
	
	def __repr__(self):
		if len(self.sequence) > 1:
			return f"This is a list of {DNA.dna_counter} DNA sequences"
		else:
			return f"This is a DNA sequence of length {len(self.sequence[0])}bp"

	def complement(self):
		"""complements a sequence of DNA"""
		comp_lst = []
		for seq in self.sequence:
			complement = [DNA._dna_pairings.get(key) for key in seq]
			complement = ''.join(complement)
			comp_lst.append(complement)
		return comp_lst

	def reverse_complement(self):
		"""reverse complements a sequence of DNA"""
		rc_lst = []
		for seq in self.sequence:
			rc = list(reversed([DNA._dna_pairings.get(key) for key in seq]))
			rc = ''.join(rc)
			rc_lst.append(rc)
		return rc_lst

	def transcribe(self):
		"""transcribe a sequence of DNA"""
		ts_lst = []
		for seq in self.sequence:
			transcribe = [DNA._dna_to_rna.get(key) for key in seq]
			transcribe = ''.join(transcribe)
			ts_lst.append(transcribe)
		return ts_lst


class RNA(NucleicAcid):
	"""define RNA, how the bases are paired, and how RNA can be converted back into DNA"""
	
	_rna = re.compile('[AUCGN]+')
	_rna_pairings = {'A': 'U', 'C': 'G', 'U': 'A', 'G': 'C', 'N': 'N'}
	_rna_to_dna = {'A': 'A', 'C': 'C', 'U': 'T', 'G': 'G', 'N': 'N'}

	rna_counter = 0

	def __init__(self, sequence):
		super().__init__(sequence)
		for seq in self.sequence:
			RNA.rna_counter += 1
			if not RNA._rna.fullmatch(seq):
				if len(self.sequence) > 1:
					raise ValueError(f"Sequence {RNA.rna_counter} is not a valid RNA sequence")
				else:
					raise ValueError("Please use a valid RNA sequence")
	
	def __repr__(self):
		if len(self.sequence) > 1:
			return f"This is a list of {RNA.rna_counter} RNA sequences"
		else:
			return f"This is a RNA sequence of length {len(self.sequence[0])}bp"

	def complement(self):
		"""complements a sequence of DNA"""
		comp_lst = []
		for seq in self.sequence:
			complement = [RNA._rna_pairings.get(key) for key in seq]
			complement = ''.join(complement)
			comp_lst.append(complement)
		return comp_lst

	def reverse_complement(self):
		"""reverse complements a sequence of DNA"""
		rc_lst = []
		for seq in self.sequence:
			rc = list(reversed([RNA._rna_pairings.get(key) for key in seq]))
			rc = ''.join(rc)
			rc_lst.append(rc)
		return rc_lst

	def reverse_transcribe(self):
		"""reverse transcribe a sequence of RNA"""
		rts_lst = []
		for seq in self.sequence:
			rev_transcribe = [RNA._rna_to_dna.get(key) for key in seq]
			rev_transcribe = ''.join(rev_transcribe)
			rts_lst.append(rev_transcribe)
		return rts_lst


class Orf(DNA):
	"""simplified DNA version that does not deal with introns; effectively the CDS; keeps track of the number of potential ORFs found, which can be reset"""

	orf_counter = 0

	_start = 'ATG'
	_stop = ['TAA', 'TAG', 'TGA']

	def __init__(self, sequence, reset = False):
		super().__init__(sequence)
		self.reset = reset
		if self.reset == True:
			self.resetter()
		for seq in self.sequence:
			Orf.orf_counter += 1
			# an ORF must have a start and stop codon as well as at least one codon in between
			if not seq.startswith(Orf._start) or not seq[-3:] in Orf._stop or not len(seq) % 3 == 0 or not len(seq) > 6:
				if len(self.sequence) > 1:
					raise ValueError(f"Sequence {Orf.orf_counter} is not a valid ORF")
				else:
					raise ValueError("Please use a valid ORF sequence")	
	
	def __repr__(self):
		if len(self.sequence) > 1:
			return f"This is a list of {Orf.orf_counter} ORFs"
		else:
			return f"This is an ORF of length {len(self.sequence[0])} bases"

	def resetter(self):
		orf_counter = 0

	def orf_lens(self):
		orf_lengths = [len(orf) for orf in self.sequence]
		print(f"There are {Orf.orf_counter} ORFs that are {orf_lengths} bases long")
		return orf_lengths


class Fasta(DNA):

	sequences_found = 0
	
	@classmethod
	def display_sequences(cls):
		return f"There are {cls.sequences_found} sequences in this FASTA file"

	@classmethod
	def preprocess_fasta(cls, fasta):
		fasta_dict = {}
		with open(fasta, 'r') as fa:
			for line in fa:
				if line.startswith('>'):
					key = line.rstrip()
					fasta_dict[key] = ''
				else:
					fasta_dict[key] = fasta_dict[key] + line.rstrip()
		return cls(fasta_dict)


		# keys_values = list(zip(fasta_keys, fasta_values))
		# fasta_dict.update(keys_values)
		# return cls(fasta_dict.keys(), fasta_dict.values())

	def __init__(self, sequence, seq_id):
		super().__init__(sequence)
		self.seq_id = seq_id
		Fasta.sequences_found += 1


# Adding this in the future
#class Fastq(DNA):

seq= Orf(['ATGAAATAA', 'ATGAAATTTTAA', 'ATGAAGGGATTTTAA'])
seq.orf_lens()