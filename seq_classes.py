#!/usr/bin/env python 
# 4-space indented, v1.0
# File name: seq_classes.py
# Description: This contains a list of types of sequences commonly encountered during genomics analyses and creates  a distinct class for each. Nucleic acids is the parent class that all other sequence classes inherit from.
# Author: Robert Linder
# Date: 2022-09-20

import re


class NucleicAcid:

	_dna= re.compile('[ATCGN]+')
	_rna= re.compile('[AUCGN]+')
	_dna_pairings = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N': 'N'}
	_rna_pairings = {'A': 'U', 'C': 'G', 'U': 'A', 'G': 'C', 'N': 'N'}
	_dna_to_rna = {'A': 'A', 'C': 'C', 'T': 'U', 'G': 'G', 'N': 'N'}
	_rna_to_dna = {'A': 'A', 'C': 'C', 'U': 'T', 'G': 'G', 'N': 'N'}

	def __init__(self, sequence):
		self.sequence = sequence.upper().replace(" ", "").replace("-", "")
		if not NucleicAcid._dna.fullmatch(self.sequence) and not NucleicAcid._rna.fullmatch(self.sequence):
			raise ValueError("Please use a valid nucleic acid sequence")

	def __repr__(self):
		return f"{self.sequence} is a nucleic acid sequence of length {len(self.sequence)}bp"


	def reverse(self):
		"""reverses a nucleic acid sequence, agnostic of type"""
		return self.sequence[::-1]

class DNA(NucleicAcid):

	def __init__(self, sequence):
		super().__init__(sequence)
		if not NucleicAcid._dna.fullmatch(self.sequence):
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
		ts = [NucleicAcid._dna_to_rna.get(key) for key in self.sequence]
		return ''.join(ts)


class RNA(NucleicAcid):

	def __init__(self, sequence):
		super().__init__(sequence)
		if not NucleicAcid._rna.fullmatch(self.sequence):
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
		rts = [NucleicAcid._rna_to_dna.get(key) for key in self.sequence]
		return ''.join(rts)


	#def transcribe(self)


#class RNA(NucleicAcid):

	#def reverse_transcribe(self)


#class Orf(DNA):

#class 

#class Fasta(DNA):


#class Fastq(DNA):

#class Primer(DNA):

#class 


nuc = RNA(sequence = "acgauu")
print(nuc.reverse_transcribe())