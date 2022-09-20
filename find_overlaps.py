#!/usr/bin/env python

import itertools

def readGenome(filename):
	genome = ''
	with open(filename, 'r') as f:
		for line in f:
			# ignore header line with genome information
			if not line[0] == '>':
				genome += line.rstrip()
	return genome

def readFastq(filename):
	sequences = []
	with open(filename) as fh:
		while True:
			fh.readline() # skip name line
			seq = fh.readline().rstrip() # read base sequence
			fh.readline() # skip placeholder line
			qual = fh.readline().rstrip() #base quality line
			if len(seq) == 0:
				break
			sequences.append(seq)
	return sequences

def editDistance(x, y):
	# Create distance matrix
	D = []
	for i in range(len(x)+1):
		D.append([0]*(len(y)+1))
	# Initialize first row and column of matrix
	for i in range(len(x)+1):
		D[i][0] = i
	for i in range(len(y)+1):
		D[0][i] = 0
	#print(D)
	# Fill in the rest of the matrix
	for i in range(1, len(x)+1):
		for j in range(1, len(y)+1):
			distHor = D[i][j-1] + 1
			distVer = D[i-1][j] + 1
			if x[i-1] == y[j-1]:
				distDiag = D[i-1][j-1]
			else:
				distDiag = D[i-1][j-1] + 1
			D[i][j] = min(distHor, distVer, distDiag)
	# Edit distance is the value in the bottom right corner of the matrix
	#print(D)
	return min(D[-1])

def checkSuffix(reads, k):
	all_overlaps = []
	read_dict = {}
	for read in reads:
		for i in range(len(read) - k + 1):
			kmer = read[i:i+k]
			if kmer in read_dict:
				read_dict[kmer].append(read)
			else:
				read_dict[kmer] = [read]
	unique_dict = {k:set(v) for k, v in read_dict.items()}
	for read_a in reads:
		suffix = read_a[-k:]
		reads_with = unique_dict[suffix]
		if len(reads_with) > 0:
			for read_b in reads_with:
				if read_a != read_b: 
					find_overlap = overlap(read_a, read_b, min_length = k)
					if not isinstance(find_overlap, int):
						all_overlaps.append(find_overlap)
					else:
						continue
				else:
					continue
	return all_overlaps

def overlap(a, b, min_length=3):
	""" Return length of longest suffix of 'a' matching
		a prefix of 'b' that is at least 'min_length'
		characters long.  If no such overlap exists,
		return 0. """
	start = 0  # start all the way at the left
	while True:
		start = a.find(b[:min_length], start)  # look for b's suffx in a
		if start == -1:  # no more occurrences to right
			return 0
		# found occurrence; check for full suffix/prefix match
		if b.startswith(a[start:]):
			return (a, b)
		start += 1  # move just past previous match

def scs(ss):
	""" Returns shortest common superstring of given
		strings, which must be the same length """
	shortest_sup = None
	all_sups = []
	for ssperm in itertools.permutations(ss):
		sup = ssperm[0]  # superstring starts as first string
		for i in range(len(ss)-1):
			# overlap adjacent strings A and B in the permutation
			olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
			# add non-overlapping portion of B to superstring
			sup += ssperm[i+1][olen:]
		all_sups.append(sup)
		if shortest_sup is None or len(sup) < len(shortest_sup):
			shortest_sup = sup # found shorter superstring
	smallest_sups = [s for s in all_sups if len(s) == len(shortest_sup)]
	return smallest_sups  # return shortest

if __name__ == '__main__':
	#genome = readGenome("chr1.GRCh38.excerpt.fasta")
	#p = "GATTTACCAGATTGAG"
	#t = "TATTGGCTATACGGTT"
	#practice = editDistance(p, genome)
	#print(practice)
	phiX = readFastq("ERR266411_1.for_asm.fastq")
	#reads = ['CGTACG', 'TACGTA', 'GTACGT', 'ACGTAC', 'GTACGA', 'TACGAT']
	find_overlaps = checkSuffix(phiX, 30)
	print(set(find_overlaps))
	print(len(set(find_overlaps)))
	all_suffs = [a[0] for a in find_overlaps]
	print(len(set(all_suffs)))
	#virus = readFastq("ads1_week4_reads.fq")
	#ss = ['GAT', 'TAG', 'TCG', 'TGC', 'AAT', 'ATA']
	#ss = ['CCT', 'CTT', 'TGC', 'TGG', 'GAT', 'ATT']
	#shortest = print(scs(ss))
