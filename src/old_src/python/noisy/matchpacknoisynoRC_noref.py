#similar to v1 in the paper but does not support RC, does both the reordering and encoding steps for real data.
#Works with reads containing 'N'. ind contains the indices in the dictionary - note that here the indices need not be 
#symmetric as is the case in the C++ implementation.
#Also, this implementation does not use bitsets and hence the thresh parameter is the actual threshold. 

#This was written because matchpacknoisyRC.py was unable to handle SRR065390 but soon we moved to C++

from distance import hamming

infile = "SRR959239.dna"
outfile_seq = "read_seq33.txt"
outfile_flag = "read_flag33.txt"
outfile_noise = "read_noise33.txt"
outfile_noisepos = "read_noisepos33.txt"
readlen = 98
no_reads = 5372832
matchlen = 80
maxmatch = 20
num_dict = 1 # should divide matchlen
thresh = 10

#ind = [[i for i in range(j,matchlen,num_dict)] for j in range(num_dict)]#
indstart = 20
indend = 60
print "Reading file"
f = open(infile,'r')
lines = [f.readline().rstrip('\n') for i in range(no_reads)]
f.close()

print "Constructing dictionaries"
d = {}
for i in range(no_reads):
	s = lines[i][indstart:indend]
	if s in d:
		d[s].append(i)
	else:
		d[s] = [i]

print "Ordering reads and writing to file"
f_seq = open(outfile_seq,'w')
f_flag = open(outfile_flag,'w')
f_noise = open(outfile_noise,'w')
f_noisepos = open(outfile_noisepos,'w')

remainingreads = set([i for i in range(no_reads)])
current = 0
prev = lines[current]
f_flag.write('0')
f_seq.write(lines[current]+'\n')

while True:
	flag = 0
	if len(remainingreads)%1000000 == 0:
		print str(len(remainingreads)//1000000)+'M reads remain'
	remainingreads.remove(current)
	s = lines[current][indstart:indend]
	d[s].remove(current)
	if len(d[s]) == 0:
		del d[s]
	if(len(remainingreads)==0):
		break
	for i in range(maxmatch):
		s = lines[current][indstart+i:indend+i]
		if s in d:
			inter = set(d[s])
			for j in inter:
				if(hamming(lines[current][i:],lines[j][:readlen-i])<=thresh):
					current = j
					currentseq = lines[j]
					flag = 1
					break

		if flag == 1:
			f_flag.write('p') #p for prev
			f_seq.write(currentseq[(readlen-i):]+'\n')
			prevj = 0;
			for j in range(readlen-i):
				if currentseq[j]!=prev[i+j]:
					f_noise.write(currentseq[j])
					f_noisepos.write("%02d"%(j-prevj))#delta encoding
					prevj = j	
			f_noise.write('\n')
			break

	if flag == 1:
		prev = currentseq
		continue
	
	current = remainingreads.pop()
	remainingreads.add(current)
	prev = lines[current]
	f_flag.write('0')
	f_seq.write(lines[current]+'\n')

print "Done"
f_seq.close()
f_flag.close()	
f_noise.close()
f_noisepos.close()

