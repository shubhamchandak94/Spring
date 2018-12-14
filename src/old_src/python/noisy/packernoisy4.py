#Similar to packernoisy2.py but encodes noise using the function encodenoise. The idea is that given the original base, the 
#noisy base has only 4 possibilities (e.g. if the original base was A, the new base can be C, G, T or N. Thus we can save
#some space after compression by using a cyclic encoding (look at the function encodenoise)
#The outfile_noise now contains numbers 1,2,3 or 4 instead of the bases.
import operator
from itertools import imap

#hamming2 function from http://code.activestate.com/recipes/499304-hamming-distance/
def hamming2(str1, str2):
    #assert len(str1) == len(str2)
    #ne = str.__ne__  ## this is surprisingly slow
    ne = operator.ne
    return sum(imap(ne, str1, str2))

def char2index(c):
	if c == 'A':
		return 0		
	if c == 'C':
		return 1
	if c == 'G':
		return 2
	if c == 'T':
		return 3
	if c == 'N':
		return 4
def findmajority(count):
	l = []
	for i in range(len(count[0])):
		s = [count[j][i] for j in range(5)]
		maxcount = max(s[0:4])
		if maxcount == 0: #only N's seen so far
			l.append('A')
			continue
		if s[0] == maxcount:
			l.append('A')
			continue
		if s[1] == maxcount:
			l.append('C')
			continue
		if s[2] == maxcount:
			l.append('G')
			continue
		if s[3] == maxcount:
			l.append('T')
			continue
	return ''.join(l)

def encodenoise(c1,c2):
	if c1 == 'A':
		if c2 == 'C':
			return 1
		if c2 == 'G':
			return 2
		if c2 == 'T':
			return 3
		if c2 == 'N':
			return 4
	if c1 == 'C':
		if c2 == 'A':
			return 1
		if c2 == 'G':
			return 2
		if c2 == 'T':
			return 3
		if c2 == 'N':
			return 4
	if c1 == 'G':
		if c2 == 'C':
			return 1
		if c2 == 'A':
			return 2
		if c2 == 'T':
			return 3
		if c2 == 'N':
			return 4
	if c1 == 'T':
		if c2 == 'C':
			return 1
		if c2 == 'G':
			return 2
		if c2 == 'A':
			return 3
		if c2 == 'N':
			return 4
	if c1 == 'N':
		if c2 == 'C':
			return 1
		if c2 == 'G':
			return 2
		if c2 == 'T':
			return 3
		if c2 == 'A':
			return 4




infile = "tempte8.dna"
outfile_seq = "read_seq32.txt"
outfile_flag = "read_flag32.txt"
outfile_noise = "read_noise32.txt"
outfile_noisepos = "read_noisepos32.txt"

readlen = 100
minmatch = 20
thresh = 20 # maximum number of mismatches allowed 

f_seq = open(outfile_seq,'w')
f_flag = open(outfile_flag,'w')
f_noise = open(outfile_noise,'w')
f_noisepos = open(outfile_noisepos,'w')
k = 0
ref = 'A'*readlen # ref is the reference which is constantly updated (introduced because matching a read to previous read leads to double noise than actual)
prev = 'A'*readlen
count = [[1]*100,[0]*100,[0]*100,[0]*100,[0]*100] #number of A's,C's,T's,G's and N's seen at each position in ref
#Note: N is never considered in the ref - we arbitrarily place an A if only N's are seen at some position
with open(infile,'r') as f:
	for line in f:
		k = k + 1
		if k%1000000 == 0:
			print str(k//1000000)+'M done'
		current = line.rstrip('\n')
		flag = 0
		for i in range(minmatch):
			if(hamming2(current[:(readlen-i)],ref[i:])<=thresh):
				if(hamming2(current[:(readlen-i)],ref[i:])<=hamming2(current[:(readlen-i)],prev[i:])):
					f_flag.write('r')
					f_seq.write(current[(readlen-i):]+'\n')
					prevj = 0;
					for j in range(readlen-i):
						count[char2index(current[j])][i+j] += 1		
						if current[j]!=ref[i+j]:
							f_noise.write(str(encodenoise(ref[i+j],current[j])))
							f_noisepos.write("%02d"%(j-prevj))#delta encoding
							prevj = j	
				else:
					f_flag.write('p')
					f_seq.write(current[(readlen-i):]+'\n')
					prevj = 0;
					for j in range(readlen-i):
						count[char2index(current[j])][i+j] += 1		
						if current[j]!=prev[i+j]:
							f_noise.write(str(encodenoise(prev[i+j],current[j])))
							f_noisepos.write("%02d"%(j-prevj))#delta encoding
							prevj = j	
				
				count = [count[j][i:]+[0]*i for j in range(5)]
				for j in range(readlen-i,readlen):
					count[char2index(current[j])][j] = 1
				
				ref = findmajority(count)	
				#ref = current#ref[i:]+current[readlen-i:]
				f_noise.write('\n')
				flag = 1
				break
		
		if flag == 0:
			f_flag.write('0')
			f_seq.write(current+'\n')
			count = [[0]*100 for j in range(5)]
			for j in range(readlen):
				count[char2index(current[j])][j] = 1
			ref = findmajority(count)
		prev = current						
f_seq.close()
f_flag.close()	
f_noise.close()
f_noisepos.close()

