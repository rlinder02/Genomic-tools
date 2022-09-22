#!/usr/bin/env python 
# 4-space indented, v1.0
# File name: seq_classes.py
# Description: This contains a list of types of sequences commonly encountered during genomics analyses and creates  a distinct class for each. Nucleic acids is the parent class that all other sequence classes inherit from.
# Author: Robert Linder
# Date: 2022-09-20

import re


class NucleicAcid:
	"""define a nucliec acid, including unnatural sequences that contain Ts and Us; can be passed in as a single string or a list of strings"""

	_nuc = re.compile('[ATCGUN]+')

	def __init__(self, sequence):
		if not isinstance(sequence, list) and not isinstance(sequence, str):
			raise TypeError("Please enter a list of sequence strings or a single sequence string")
		self.sequence = [seq.upper().replace(" ", "").replace("-", "") for seq in sequence]
		seq_counter = 0
		for seq in self.sequence:
			seq_counter += 1
			if not NucleicAcid._nuc.fullmatch(self.seq):
				if isinstance(self.sequence, list):
					raise ValueError(f"Sequence {seq_counter} is not a valid nucleic acid sequence")
				else:
					raise ValueError("Please use a valid nucleic acid sequence")

		self.sequence = sequence.upper().replace(" ", "").replace("-", "")
		if not NucleicAcid._nuc.fullmatch(self.sequence):
			raise ValueError("Please use a valid nucleic acid sequence")

	def __repr__(self):
		return f"{self.sequence} is a nucleic acid sequence of length {len(self.sequence)}bp"

	def reverse(self):
		"""reverses a nucleic acid sequence, agnostic of type"""
		return self.sequence[::-1]

class DNA(NucleicAcid):
	"""define DNA, how the bases are paired, and how DNA can be converted into RNA"""
	
	_dna = re.compile('[ATCGN]+')
	_dna_pairings = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N': 'N'}
	_dna_to_rna = {'A': 'A', 'C': 'C', 'T': 'U', 'G': 'G', 'N': 'N'}

	def __init__(self, sequence):
		super().__init__(sequence)
		if not DNA._dna.fullmatch(self.sequence):
			raise ValueError("Please use a valid DNA sequence")
	
	def __repr__(self):
		return f"{self.sequence} is a DNA sequence of length {len(self.sequence)}bp"		

	def complement(self):
		"""complements a sequence of DNA"""
		complement = [DNA._dna_pairings.get(key) for key in self.sequence]
		return ''.join(complement)

	def reverse_complement(self):
		"""reverse complements a sequence of DNA"""
		rc = list(reversed([DNA._dna_pairings.get(key) for key in self.sequence]))
		return ''.join(rc)

	def transcribe(self):
		"""transcribe a sequence of DNA"""
		ts = [DNA._dna_to_rna.get(key) for key in self.sequence]
		return ''.join(ts)


class RNA(NucleicAcid):
	"""define RNA, how the bases are paired, and how RNA can be converted back into DNA"""
	
	_rna = re.compile('[AUCGN]+')
	_rna_pairings = {'A': 'U', 'C': 'G', 'U': 'A', 'G': 'C', 'N': 'N'}
	_rna_to_dna = {'A': 'A', 'C': 'C', 'U': 'T', 'G': 'G', 'N': 'N'}

	def __init__(self, sequence):
		super().__init__(sequence)
		if not RNA._rna.fullmatch(self.sequence):
			raise ValueError("Please use a valid RNA sequence")
	
	def __repr__(self):
		return f"{self.sequence} is a RNA sequence of length {len(self.sequence)}bp"		

	def complement(self):
		"""complements a sequence of DNA"""
		complement = [RNA._rna_pairings.get(key) for key in self.sequence]
		return ''.join(complement)

	def reverse_complement(self):
		"""reverse complements a sequence of DNA"""
		rc = list(reversed([RNA._rna_pairings.get(key) for key in self.sequence]))
		return ''.join(rc)

	def reverse_transcribe(self):
		"""reverse transcribe a sequence of RNA"""
		rts = [RNA._rna_to_dna.get(key) for key in self.sequence]
		return ''.join(rts)


class Orf(DNA):
	"""simplified version that does not deal with introns; effectively the CDS; keeps track of the number of potential ORFs found, which can be reset"""

	orfs_found = 0

	_start = 'ATG'
	_stop = ['TAA', 'TAG', 'TGA']

	@classmethod
	def display_potential_orfs(cls):
		return f"There are {cls.orfs_found} potential ORFs"

	def __init__(self, sequence, reset = False):
		super().__init__(sequence)
		self.reset = reset
		if not self.sequence.startswith(Orf._start) or not self.sequence[-3:] in Orf._stop or not len(self.sequence) % 3 == 0:
			raise ValueError("Please use a valid ORF sequence")
		if self.reset == True:
			self.resetter()
		else:
			Orf.orfs_found += 1

	def __repr__(self):
		return f"{self.sequence} is an ORF of length {len(self.sequence)} bases"

	def resetter(self):
		orfs_found = 0


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

		cls(fasta_seqs, fasta_headers)


		# keys_values = list(zip(fasta_keys, fasta_values))
		# fasta_dict.update(keys_values)
		# return cls(fasta_dict.keys(), fasta_dict.values())

	def __init__(self, sequence, seq_id):
		super().__init__(sequence)
		self.seq_id = seq_id
		Fasta.sequences_found += 1


# Adding this in the future
#class Fastq(DNA):

orf = Orf('ATGUUUTAA')
print(orf)

