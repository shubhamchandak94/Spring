#Reordering for noisy reads without RC. Uses 4 dictionaries (defined in ind) but does not have any thresh (basically means
#thresh is very large. Any read in any of the dictionaries is picked.


infile = "chrom22_50x_noRC_noisy.dna"
outfile = "tempte.dna"
readlen = 100
no_reads = 17500000
maxmatch = 28
ind = [[i for i in range(j,72,4)] for j in range(4)]

print "Reading file"
f = open(infile,'r')
lines = [f.readline().rstrip('\n') for i in range(no_reads)]
f.close()

print "Constructing dictionaries"
d = [{} for i in range(4)]
for i in range(no_reads):
	l = [[lines[i][j] for j in ind[k]] for k in range(4)]
	s = [''.join(l[k]) for k in range(4)]
	for k in range(4):
		if s[k] in d[k]:
			d[k][s[k]].add(i)
		else:
			d[k][s[k]] = set([i])

print "Ordering reads and writing to file"
remainingreads = set([i for i in range(no_reads)])
current = 0
fout = open(outfile,'w')
while True:
	flag = 0
	if len(remainingreads)%1000000 == 0:
		print str(len(remainingreads)//1000000)+'M reads remain'
	fout.write(lines[current]+'\n')
	remainingreads.remove(current)
	if(len(remainingreads)==0):
		break
	for i in range(maxmatch):
		l = [[lines[current][j+i] for j in ind[k]] for k in range(4)]
		s = [''.join(l[k]) for k in range(4)]
		for k in range(4):	
			if s[k] in d[k]:
				if len(d[k][s[k]].intersection(remainingreads))>0:
					current = (d[k][s[k]].intersection(remainingreads)).pop()
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
