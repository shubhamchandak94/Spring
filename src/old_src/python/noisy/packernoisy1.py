#Encoding for noisy reads. Encodes each read in terms of clean reference (maintained using findmajority function and count array) 
#Four output files - outfile_seq has the suffix (and newline), outfile_flag has flag (+ or 0)
#outfile_noise has the base in the current read and a newline after every read that was matched,
#noisepos file has the noise position delta encoded with 2 digits per position (no newline here)
#Note that the clean reference never has the N character.

from distance import hamming

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



infile = "tempte1.dna"
outfile_seq = "read_seq23.txt"
outfile_pos = "read_pos23.txt"
outfile_noise = "read_noise23.txt"
outfile_noisepos = "read_noisepos23.txt"

readlen = 100
minmatch = 20
thresh = 20 # maximum number of mismatches allowed 

f_seq = open(outfile_seq,'w')
f_pos = open(outfile_pos,'w')
f_noise = open(outfile_noise,'w')
f_noisepos = open(outfile_noisepos,'w')
k = 0
ref = 'A'*readlen # ref is the reference which is constantly updated (introduced because matching a read to previous read leads to double noise than actual)
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
				f_pos.write('+')
				f_seq.write(current[(readlen-i):]+'\n')
				prevj = 0;
				for j in range(readlen-i):
					count[char2index(current[j])][i+j] += 1		
					if current[j]!=ref[i+j]:
						f_noise.write(current[j])
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
			f_pos.write('0')
			f_seq.write(current+'\n')
			count = [[0]*100 for j in range(5)]
			for j in range(readlen):
				count[char2index(current[j])][j] = 1
			ref = findmajority(count)					
f_seq.close()
f_pos.close()	
f_noise.close()
f_noisepos.close()

