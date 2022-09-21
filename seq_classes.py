#!/usr/bin/env python 
# 4-space indented, v1.0
# File name: seq_classes.py
# Description: This contains a list of types of sequences commonly encountered during genomics analyses and creates  a distinct class for each. Nucleic acids is the parent class that all other sequence classes inherit from.
# Author: Robert Linder
# Date: 2022-09-20


class NucleicAcid:


	def __init__(self, letters):
		self.letters = letters.upper().replace(" ", "").replace("-", "")

	def __repr__(self):
		return f"{self.letters} is a nucleic acid"

	@letters.setter
	def letter(self, )



#class DNA(NucleicAcid):

	#def transcribe(self)


#class RNA(NucleicAcid):

	#def reverse_transcribe(self)


#class Orf(DNA):

#class 

#class Fasta(DNA):


#class Fastq(DNA):

#class Primer(DNA):

#class 


nuc = NucleicAcid(letters = "act g")
print(nuc)