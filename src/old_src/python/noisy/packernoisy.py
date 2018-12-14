#Encoding for noisy reads. Compares each read to the previous read (no clean reference maintained). 
#Four output files - outfile_seq has the suffix (and newline), outfile_flag has flag (+ or 0)
#outfile_noise has the base in the current read and a newline after every read that was matched,
#noisepos file has the noise position stored in 2 digits per position (no newline here)

from distance import hamming

infile = "tempte.dna"
outfile_seq = "read_seq16.txt"
outfile_flag = "read_flag16.txt"
outfile_noise = "read_noise16.txt"
outfile_noisepos = "read_noisepos16.txt"

readlen = 100
maxmatch = 28
thresh = 20 # maximum number of mismatches allowed 

f_seq = open(outfile_seq,'w')
f_flag = open(outfile_flag,'w')
f_noise = open(outfile_noise,'w')
f_noisepos = open(outfile_noisepos,'w')
k = 0
prev = 'A'*readlen
with open(infile,'r') as f:
	for line in f:
		k = k + 1
		if k%1000000 == 0:
			print str(k//1000000)+'M done'
		current = line.rstrip('\n')
		flag = 0
		for i in range(maxmatch):
			if(hamming(current[:(readlen-i)],prev[i:])<=thresh):
				f_flag.write('+')
				f_seq.write(current[(readlen-i):]+'\n')
				for j in range(readlen-i):
					if current[j]!=prev[i+j]:
						f_noise.write(current[j])
						f_noisepos.write("%02d"%j)
				prev = current
				f_noise.write('\n')
				flag = 1
				break
		
		if flag == 0:
			f_flag.write('0')
			f_seq.write(current+'\n')
			prev = current
f_seq.close()
f_flag.close()	
f_noise.close()
f_noisepos.close()

