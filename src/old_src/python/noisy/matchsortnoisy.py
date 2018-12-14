#Reordering for noisy reads without RC. Uses single dictionary 0:matchlen but has a thresh while matching. 

#Does not do too well. Lots of Hreads.

from distance import hamming

infile = "chrom22_50x_noRC_noisy.dna"
outfile = "tempte.dna"
readlen = 100
no_reads = 17500000
matchlen = 15
maxmatch = 10
thresh = 8 #maximum number of mismatches allowed

print "Reading file"
f = open(infile,'r')
lines = [f.readline().rstrip('\n') for i in range(no_reads)]
f.close()

print "Constructing dictionary"
d = {}
for i in range(no_reads):
	if lines[i][0:matchlen] in d:
		d[lines[i][0:matchlen]].add(i)
	else:
		d[lines[i][0:matchlen]] = set([i])

print "Ordering reads and writing to file"
remainingreads = set([i for i in range(no_reads)])
current = 0
fout = open(outfile,'w')
while True:
	flag = 0
	if len(remainingreads)%1000000 == 0:
		print str(len(remainingreads)//1000000)+'M reads remain'
	fout.write(lines[current]+'\n')
	d[lines[current][0:matchlen]].remove(current)
	remainingreads.remove(current)
	if(len(remainingreads)==0):
		break
	if len(d[lines[current][0:matchlen]]) == 0:
		del d[lines[current][0:matchlen]]
	else:
		for i in d[lines[current][0:matchlen]]:
			if hamming(lines[current][matchlen:],lines[i][matchlen:]) <= thresh:
				current = i
				flag = 1
				break
	if flag == 1:
		continue
	for j in range(1,maxmatch):
		if lines[current][j:j+matchlen] in d:
			for i in d[lines[current][j:j+matchlen]]:
				if hamming(lines[current][j+matchlen:],lines[i][matchlen:readlen-j]) <= thresh:
					current = i
					flag = 1
					break
			if flag == 1:
					break
	if flag == 1:
		continue
	current = remainingreads.pop()
	remainingreads.add(current)

print "Done"
fout.close()
