#Encoding stage for noiseless reads. 
#infile - file containing reordered reads produced by reordernoiseless.cpp 
#outfile_seq - the encoded file contaning suffixes. Compress using xz.
#outfile_flag - the encoded file contaning flags. Compress using xz.
#readlen - length of reads (assumed constant)
#maxmatch - maxshift in paper (same as that for reordering stage)

infile = "temp3.dna"
outfile_seq = "read_seq31.txt"
outfile_flag = "read_flag31.txt"
readlen = 100
maxmatch = 20
f_seq = open(outfile_seq,'w')
f_flag = open(outfile_flag,'w')
prev = 'A'*readlen
k = 0
with open(infile,'r') as f:
	for line in f:
		k = k + 1
		if (k%1000000 == 0):
			print str(k//1000000)+"M reads done"
		current = line.rstrip('\n')
		flag = 0
		for i in range(maxmatch):
			if(current[:readlen-i]==prev[i:]):
				f_flag.write('+')
				f_seq.write(current[(readlen-i):]+'\n')
				prev = current
				flag = 1
				break
		
		if flag == 0:
			f_flag.write('0')
			f_seq.write(current+'\n')
			prev = current
f_seq.close()
f_flag.close()			
