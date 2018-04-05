#Encoding stage for noisy reads. 
#infile - file containing reordered reads produced by reordernoisy.cpp 
#infile_flag - file containing flags produced by reordernoisy.cpp
#outfile_... - the encoded files. Compress these and the rev file produced in stage 1 using xz.
#readlen - length of reads (assumed constant)
#maxmatch - maxshift in paper (same as that for reordering stage)
#thresh - Hamming threshold - put it equal to that used in reordering stage.

import sys
import os

sys.path.append( os.path.abspath('.') )
from config import *

char2index = {'A':0,'C':1,'G':2,'T':3,'N':4}
index2char = {0:'A',1:'C',2:'G',3:'T',4:'N'}

# maxmatch = 30
# thresh = 24 # maximum number of mismatches allowed 

basename = sys.argv[1]
basedir = os.path.join(basename,"output") 
infile = os.path.join(basedir,"temp.dna")
infile_flag = os.path.join(basedir,"tempflag.txt")
infile_pos = os.path.join(basedir,"temppos.txt")
outfile_seq = os.path.join(basedir,"read_seq.txt")
outfile_meta = os.path.join(basedir,"read_meta.txt")
outfile_pos = os.path.join(basedir,"read_pos.txt")
outfile_noise = os.path.join(basedir,"read_noise.txt")
outfile_noisepos = os.path.join(basedir,"read_noisepos.txt")


in_flag = open(infile_flag,'r')
in_pos = open(infile_pos,'r')
f_seq = open(outfile_seq,'w')
f_meta = open(outfile_meta,'w')
f_pos = open(outfile_pos,'w')
f_noise = open(outfile_noise,'w')
f_noisepos = open(outfile_noisepos,'w')
k = 0

readlen = 0
with open(infile,'r') as f:
	for line in f:
		current = line.rstrip('\n')
		readlen = len(current)
		break;
print("readlen: ", readlen)



#inttoascii = {0:'a',1:'b',2:'c',3:'d',4:'e',5:'f',6:'g',7:'h',8:'i',9:'j',10:'k',11:'l',12:'m',13:'n',14:'o',15:'p',16:'q',17:'r',18:'s',19:'t',20:'u',21:'w',22:'x',23:'y',24:'z',25:'A',26:'B',27:'C',28:'D',29:'E',30:'F',31:'G',32:'H',33:'I',34:'J',35:'K',36:'L',37:'M',38:'N',39:'O',40:'P',readlen:'v'}

#noise encoding - uses substitution statistics from Minoche et al.
enc_noise = {
('A','C'):'0',('A','G'):'1',('A','T'):'2',('A','N'):'3',
('C','A'):'0',('C','G'):'1',('C','T'):'2',('C','N'):'3',
('G','T'):'0',('G','A'):'1',('G','C'):'2',('G','N'):'3',
('T','G'):'0',('T','C'):'1',('T','A'):'2',('T','N'):'3',
('N','A'):'0',('N','G'):'1',('N','C'):'2',('N','T'):'3',
}

def findmajority(count):
	maxcount = [max(s) for s in count]
	l = [index2char[s.index(maxcount[i])] for i,s in zip(range(len(count)),count)]
	return ''.join(l)

def buildcontig(reads,pos):
	if(len(reads) == 1): #singleton read
		return reads[0]
	count = [[0,0,0,0,0] for i in range(readlen)] #number of A's,C's,T's,G's,N's seen at each position in ref
	for i in range(readlen):
		count[i][char2index[reads[0][i]]] = 1
	prevpos = 0 
	for p,currentread in zip(pos[1:],reads[1:]):
		currentpos = prevpos + p
		count = count + [[0,0,0,0,0] for j in range(p)]
		for j in range(readlen):
			count[currentpos+j][char2index[currentread[j]]] += 1
		prevpos = currentpos
#		flag = 0
#		bestmatch = readlen
#		besti = 0
#		for i in range(maxmatch):
#			hammingdist = hamming(currentread[:(readlen-i)],prevread[i:])
#			if(hammingdist<=thresh):
#				pos.append(i+pos[-1])
#				count = count + [[0,0,0,0,0] for j in range(i)]
#				for j in range(readlen):
#					count[pos[-1]+j][char2index[currentread[j]]] += 1
#				flag = 1
#				break
#			if(hammingdist < bestmatch):
#				bestmatch = hammingdist
#				bestmatchpos = i
#		if flag == 0: #no match found due to some reason (this might happen because of matchsort9's
#		# handling of N's (if only N's are seen at a position, matchsort9 makes it A in the ref.
#			pos.append(bestmatchpos+pos[-1])
#			count = count + [[0,0,0,0,0] for j in range(bestmatchpos)]
#			for j in range(readlen):
#				count[pos[-1]+j][char2index[currentread[j]]] += 1
#		prevread = currentread	
	ref = findmajority(count)
	return ref

def writecontig(ref,pos,reads):
	f_seq.write(ref)
	if len(reads) == 1: #singleton read
		f_noise.write('\n')
		f_pos.write(str(unichr(pos[0])))
		return
	currentread = reads[0]
	prevj = 0
	for j in range(readlen):#first read (special case to handle pos[0] = readlen)
		if currentread[j]!=ref[j]:
			f_noise.write(enc_noise[(ref[j],currentread[j])])
			f_noisepos.write(str(unichr(j-prevj)))#delta encoding
			prevj = j
	f_noise.write('\n')
	f_pos.write(str(unichr(pos[0])))
	prevpos = 0

	for p,currentread in zip(pos[1:],reads[1:]):
		currentpos = prevpos + p	
		prevj = 0;	
		for j in range(readlen):
			if currentread[j]!=ref[currentpos+j]:
				f_noise.write(enc_noise[(ref[currentpos+j],currentread[j])])
				f_noisepos.write(str(unichr(j-prevj)))#delta encoding
				prevj = j
		f_noise.write('\n')
		f_pos.write(str(unichr(p)))
		prevpos = currentpos
	return		


reads = []
pos = []
with open(infile,'r') as f:
	for line in f:
		k = k + 1
		if k%1000000 == 0:
			print str(k//1000000)+'M done'
		current = line.rstrip('\n')
		c = in_flag.read(1)
		p = int(in_pos.readline())
		if c=='0':
			if len(reads) != 0:
				ref = buildcontig(reads,pos)
				writecontig(ref,pos,reads)
			reads = [current]
			pos = [p]
		else:
			reads.append(current)
			pos.append(p)

#last contig
ref = buildcontig(reads,pos)
writecontig(ref,pos,reads)

# write metadata
f_meta.write(str(readlen) + "\n")
f_seq.close()
f_meta.close()
f_pos.close()	
f_noise.close()
f_noisepos.close()
in_flag.close()
in_pos.close()
