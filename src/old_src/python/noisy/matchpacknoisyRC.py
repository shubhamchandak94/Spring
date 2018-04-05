#implements v2 in the paper, does both the reordering and encoding steps for real data. Works with reads containing 'N'. ind contains the 
#indices in the dictionary - note that here the indices need not be symmetric as is the case in the C++ implementation.
#Also, this implementation does not use bitsets and hence the thresh parameter is the actual threshold. 

#Is much slower and memory-intensive as compared to the C++ implmentation.

#Produces 5 output files as mentioned in the report (encoding stage)

from distance import hamming
from Bio.Seq import Seq

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

infile = "SRR959239.dna"
outfile_seq = "read_seq32.txt"
outfile_flag = "read_flag32.txt"
outfile_rev = "read_rev32.txt"
outfile_noise = "read_noise32.txt"
outfile_noisepos = "read_noisepos32.txt"
readlen = 98
no_reads = 5372832
#matchlen = 80
maxmatch = 18
thresh = 10 #maximum number of mismatches allowed
num_dict = 1 # should divide matchlen

#ind = [[i for i in range(j,matchlen,num_dict)] for j in range(num_dict)]ind = [[i for i in range(30,70)]]
ind = [[i for i in range(35,65)]]
f = open(infile,'r')
lines = [f.readline().rstrip('\n') for i in range(no_reads)]
f.close()

print "Generating reverse complements"
revlines = [str(Seq(lines[i]).reverse_complement()) for i in range(no_reads)]

print "Constructing dictionaries"
d1 = [{} for i in range(num_dict)]
d2 = [{} for i in range(num_dict)] #reverse complement dictionary
for i in range(no_reads):
	l = [[lines[i][j] for j in ind[k]] for k in range(num_dict)]
	s = [''.join(l[k]) for k in range(num_dict)]
	for k in range(num_dict):
		if s[k] in d1[k]:
			d1[k][s[k]].append(i)
		else:
			d1[k][s[k]] = [i]

	l = [[revlines[i][j] for j in ind[k]] for k in range(num_dict)]
	s = [''.join(l[k]) for k in range(num_dict)]
	for k in range(num_dict):
		if s[k] in d2[k]:
			d2[k][s[k]].append(i)
		else:
			d2[k][s[k]] = [i]

print "Ordering reads and writing to file"
f_seq = open(outfile_seq,'w')
f_flag = open(outfile_flag,'w')
f_noise = open(outfile_noise,'w')
f_noisepos = open(outfile_noisepos,'w')
f_rev = open(outfile_rev,'w')

remainingreads = set([i for i in range(no_reads)])
current = 0
prev = lines[current]
count = [[0]*readlen for j in range(5)]
for j in range(readlen):
	count[char2index(lines[current][j])][j] = 1
f_flag.write('0')
f_seq.write(lines[current]+'\n')

ref = findmajority(count)

while True:
	flag = 0
	if len(remainingreads)%1000000 == 0:
		print str(len(remainingreads)//1000000)+'M reads remain'
	remainingreads.remove(current)
	if(len(remainingreads)==0):
		break
	for i in range(maxmatch):
		l = [[ref[j+i] for j in ind[k]] for k in range(num_dict)]
		s = [''.join(l[k]) for k in range(num_dict)]
		inter = set()
		for k in range(num_dict):	
			if s[k] in d1[k]:
				inter = inter.union(set(d1[k][s[k]]))
		inter = inter.intersection(remainingreads)
		if len(inter)>0:
			for j in inter:
				if(hamming(ref[i:],lines[j][:readlen-i])<=thresh):
					f_rev.write('d') # d for direct
					current = j
					currentseq = lines[j]
					flag = 1
					break

		if flag == 0:
			inter = set()
			for k in range(num_dict):	
				if s[k] in d2[k]:
					inter = inter.union(set(d2[k][s[k]]))
			inter = inter.intersection(remainingreads)
			if len(inter)>0:
				for j in inter:
					if(hamming(ref[i:],revlines[j][:readlen-i])<=thresh):
						f_rev.write('r') # r for reverse
						current = j
						currentseq = revlines[j]
						flag = 1
						break
		if flag == 1:
			if(hamming(currentseq[:(readlen-i)],ref[i:])<=hamming(currentseq[:(readlen-i)],prev[i:])):
				f_flag.write('r') #r for ref
				f_seq.write(currentseq[(readlen-i):]+'\n')
				prevj = 0;
				for j in range(readlen-i):
					count[char2index(currentseq[j])][i+j] += 1		
					if currentseq[j]!=ref[i+j]:
						f_noise.write(currentseq[j])
						f_noisepos.write("%02d"%(j-prevj))#delta encoding
						prevj = j	
			else:
				f_flag.write('p') #p for prev
				f_seq.write(currentseq[(readlen-i):]+'\n')
				prevj = 0;
				for j in range(readlen-i):
					count[char2index(currentseq[j])][i+j] += 1		
					if currentseq[j]!=prev[i+j]:
						f_noise.write(currentseq[j])
						f_noisepos.write("%02d"%(j-prevj))#delta encoding
						prevj = j	
			
			count = [count[j][i:]+[0]*i for j in range(5)]
			for j in range(readlen-i,readlen):
				count[char2index(currentseq[j])][j] = 1
			
			ref = findmajority(count)	
			f_noise.write('\n')
			break

	if flag == 1:
		prev = currentseq
		continue
	
	current = remainingreads.pop()
	remainingreads.add(current)
	prev = lines[current]
	count = [[0]*readlen for j in range(5)]
	for j in range(readlen):
		count[char2index(lines[current][j])][j] = 1
	f_flag.write('0')
	f_seq.write(lines[current]+'\n')

	ref = findmajority(count)

print "Done"
f_seq.close()
f_flag.close()	
f_rev.close()
f_noise.close()
f_noisepos.close()

