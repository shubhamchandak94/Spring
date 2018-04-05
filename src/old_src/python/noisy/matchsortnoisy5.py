#Reordering for noisy reads without RC. Uses best match within bin for a given shift (if within thresh). 4 dictionaries (ind). 
#Sort of similar to matchsort4.cpp

from distance import hamming

infile = "chrom22_50x_noRC_noisy.dna"
outfile = "tempte1.dna"
readlen = 100
no_reads = 17500000
matchlen = 80
maxmatch = 20
thresh = 8 #maximum number of mismatches allowed in matchlen to readlen
ind = [[i for i in range(j,80,4)] for j in range(4)]

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
			d[k][s[k]].append(i)
		else:
			d[k][s[k]] = [i]

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
		bestmatch = 100
		inter = set()
		for k in range(4):	
			if s[k] in d[k]:
				inter = inter.union(set(d[k][s[k]]))
		inter = inter.intersection(remainingreads)
		if len(inter)>0:
			for j in inter:
				if(hamming(lines[current][i:],lines[j][:readlen-i])<bestmatch):
					bestmatch = hamming(lines[current][i:],lines[j][:readlen-i])
					bestmatchindex = j
		if bestmatch <= thresh:
			current = bestmatchindex
			flag = 1
			break
	if flag == 1:
		continue
	current = remainingreads.pop()
	remainingreads.add(current)

print "Done"
fout.close()
