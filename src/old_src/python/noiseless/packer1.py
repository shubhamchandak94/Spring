#Second stage encoding for noiseless reads. The variable neg decides whether we want to match the current read with previous 
#read shifted in both directions. 
#The reordered reads in infile are encoded into outfile_seq (which stores the suffix of the read/full read for unmatched read), 
#outfile_flag which stores 0 for unmatched, + for shifted right and - for shifted left, 
#Instead of \n in read_seq we have a separate file outfile_pos which stores the offset from previous read

#We expected some improvement in the compression, but the result is slightly worse than packer.py

infile = "temp2.dna"
outfile_seq = "readseq.txt"
outfile_match = "readmatch.txt"
outfile_pos = "readpos.txt"

readlen = 100
minmatch = 20
neg = False #whether we want to match the current read with the previous read shifted rightward 
f_seq = open(outfile_seq,'w')
f_match = open(outfile_match,'w')
f_pos = open(outfile_pos,'w')
prev = 'A'*readlen
with open(infile,'r') as f:
	for line in f:
		current = line.rstrip('\n')
		flag = 0
		for i in range(minmatch):
			if(current[:readlen-i]==prev[i:]):
				f_match.write('+')
				f_seq.write(current[(readlen-i):])
				f_pos.write("%02d"%i)
				prev = current
				flag = 1
				break
		if flag == 0 and neg == True:
			for i in range(minmatch):
				if(current[i:]==prev[:(readlen-i)]):
					f_match.write('-')
					f_seq.write(current[:i])
					f_pos.write("%02d"%i)
					prev = current
					flag = 1
					break
		
		if flag == 0:
			f_match.write('0')
			f_seq.write(current)
			prev = current
f_seq.close()
f_pos.close()	
f_match.close()		

