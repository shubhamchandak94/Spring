#!/usr/bin/python2
## computes noise_entropy by assuming per position

import sys
import os
import numpy as np
from tqdm import tqdm


###################################
def quality_to_prob(qual_string):
	_q = [ord(c) for c in qual_string]

	_q = (33.0 - np.array(_q))/10.0
	prob = np.power(10.0,_q)
	return prob


###################################
def get_M_readlen(inFile):
	M = 0
	N = 0
	with open(inFile,'r') as f:
	  for line in f:
	  	M += 1
	  	N = len(line) -1
	return (M,N)

###################################
def compute_entropy(probs):
	_e = -( probs*np.log(probs) + (1-probs)*np.log(1-probs) + probs*np.log(3.0) ) 
	_e = _e/np.log(2.0)
	entropy = np.sum(_e)
	return entropy,_e 

def get_genome_size(genomeFile_uncompressed):
	genome_size = 0
	with open(genomeFile_uncompressed,'r') as f:
		for i,line in enumerate(f):
			if i != 0:
				genome_size += len(line) - 1
	return genome_size
###################################
def main():
	inFile = sys.argv[1]
	genomeFile = sys.argv[2]
	genomeFile_uncompressed = sys.argv[3]
	outputFile = sys.argv[4]

	print("FASTQfile: ", inFile)
	M,readlen = get_M_readlen(inFile)
	N = get_genome_size(genomeFile_uncompressed)

	print "Read Length: ", readlen
	probs_sum = np.zeros(readlen)
	with open(inFile,'r') as f:
		for line in tqdm(f,total=M,desc="computing read probs ...",ascii=True):
			line = line.rstrip('\n')
			_p = quality_to_prob(line)
			probs_sum += _p
			# print line
			# print _p

	probs = probs_sum/(1.0*M)
	readlen = len(probs)
	entropy,_e = compute_entropy(probs)
	noise_bpb = entropy/readlen
	noise_size = entropy*M/8

	print "Error Probability: "
	print probs
	print "Entropy per position: "
	print _e

	print "Number of Reads: ", M
	print "Readlength: ", readlen
	print "Genome Size: ", N
	print "Noise bits per base:", noise_bpb
	print "Noise Size (in bytes): ", noise_size

	FASTA_statinfo = os.stat(genomeFile)
	FASTA_Size = FASTA_statinfo.st_size
	FASTA_bpb = FASTA_Size*8.0/N
	print "FASTA Size (in bytes): ", FASTA_Size
	print "FASTA bits per base: ", FASTA_bpb

	multinomial_size = ( (M+N)*np.log(M+N) - M*np.log(M) - N*np.log(N) )/(8*np.log(2))
	print "Multinomial Size (in bytes): ", multinomial_size

	total_size = multinomial_size + FASTA_Size + noise_size
	bits_per_base = total_size*8.0/(M*readlen)
	print "Total Size (in bytes): ", total_size
	print "Bits per base: ", bits_per_base

	output_statinfo = os.stat(outputFile)
	output_size = output_statinfo.st_size
	output_bpb = output_size*8.0/(M*readlen)

	print "Compressed Output Size (in bytes): ", output_size
	print "Compressed Output bits per base: ", output_bpb

###################################
if __name__ == '__main__':
    main()

