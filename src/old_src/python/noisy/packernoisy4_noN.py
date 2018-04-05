#Similar to packernoisy4.py but for files with no Ns. Is much faster as it uses better findmajority and also needs the flag
#file to be provided as input (see packernoisy2_noN.py for more on these points). 
#Noise is now encoded as 1,2 or 3. The mapping to 1,2,3 is based on the paper on Error Profiles. For each original base,
#the most likely changed base is encoded as a 1 and so on.

#Leads to small savings in the noise file. 
from distance import hamming

char2index = {'A':0,'C':1,'G':2,'T':3}
index2char = {0:'A',1:'C',2:'G',3:'T'}

def findmajority(count):
	maxcount = [max(s) for s in count]
	l = [index2char[s.index(maxcount[i])] for i,s in zip(range(readlen),count)]
	return ''.join(l)



infile = "temp8.dna"
infile_flag = "tempflag8.txt"
outfile_seq = "read_seq82.txt"
outfile_flag = "read_flag82.txt"
outfile_noise = "read_noise82.txt"
outfile_noisepos = "read_noisepos82.txt"

readlen = 100
maxmatch = 20
thresh = 24 # maximum number of mismatches allowed 

def encodenoise(c1,c2):
	if c1 == 'A':
		if c2 == 'C':
			return 1
		if c2 == 'G':
			return 2
		if c2 == 'T':
			return 3
	if c1 == 'C':
		if c2 == 'A':
			return 1
		if c2 == 'G':
			return 2
		if c2 == 'T':
			return 3
	if c1 == 'G':
		if c2 == 'T':
			return 1
		if c2 == 'C':
			return 2
		if c2 == 'A':
			return 3
	if c1 == 'T':
		if c2 == 'G':
			return 1
		if c2 == 'C':
			return 2
		if c2 == 'A':
			return 3

in_flag = open(infile_flag,'r')
f_seq = open(outfile_seq,'w')
f_flag = open(outfile_flag,'w')
f_noise = open(outfile_noise,'w')
f_noisepos = open(outfile_noisepos,'w')
k = 0
ref = 'A'*readlen # ref is the reference which is constantly updated (introduced because matching a read to previous read leads to double noise than actual)
prev = 'A'*readlen
count = [[1,0,0,0] for i in range(readlen)] #number of A's,C's,T's,G's seen at each position in ref

with open(infile,'r') as f:
	for line in f:
		k = k + 1
		if k%1000000 == 0:
			print str(k//1000000)+'M done'
		current = line.rstrip('\n')
		c = in_flag.read(1)
		flag = 0
		if c=='1':
			for i in range(maxmatch):
				if(hamming(current[:(readlen-i)],ref[i:])<=thresh or hamming(current[:(readlen-i)],prev[i:])<=thresh):
					if(hamming(current[:(readlen-i)],ref[i:])<=hamming(current[:(readlen-i)],prev[i:])):
						f_flag.write('r')
						f_seq.write(current[(readlen-i):]+'\n')
						prevj = 0;	
						for j in range(readlen-i):
							count[i+j][char2index[current[j]]] += 1		
							if current[j]!=ref[i+j]:
								f_noise.write(str(encodenoise(ref[i+j],current[j])))
								f_noisepos.write("%02d"%(j-prevj))#delta encoding
								prevj = j	
					else:
						f_flag.write('p')
						f_seq.write(current[(readlen-i):]+'\n')
						prevj = 0;
						
						for j in range(readlen-i):
							count[i+j][char2index[current[j]]] += 1		
							if current[j]!=prev[i+j]:
								f_noise.write(str(encodenoise(prev[i+j],current[j])))
								f_noisepos.write("%02d"%(j-prevj))#delta encoding
								prevj = j	
					
					count = count[i:]+[[0,0,0,0] for t in range(i)]
					for j in range(readlen-i,readlen):
						count[j][char2index[current[j]]] = 1
					ref = findmajority(count)	
					#ref = current#ref[i:]+current[readlen-i:]
					f_noise.write('\n')
					flag = 1
					break
			
			if flag == 0:
				f_flag.write('0')
				f_seq.write(current+'\n')
				count = [[0,0,0,0] for i in range(readlen)]
				for j in range(readlen):
					count[j][char2index[current[j]]] = 1
				ref = current
			prev = current						
		else:
			f_flag.write('0')
			f_seq.write(current+'\n')
			count = [[0,0,0,0] for i in range(readlen)]
			for j in range(readlen):
				count[j][char2index[current[j]]] = 1
			ref = current
			prev = current						
f_seq.close()
f_flag.close()	
f_noise.close()
f_noisepos.close()

