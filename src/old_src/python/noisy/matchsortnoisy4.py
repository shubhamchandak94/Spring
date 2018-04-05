#Reordering for noisy reads without RC. Similar to matchsort2 but has 5 dictionaries instead of 4.

from distance import hamming
	
infile = "SRR959239.dna"
outfile = "temp3.dna"
readlen = 98
no_reads = 5372832
matchlen = 80
maxmatch = 4
thresh = 8 #maximum number of mismatches allowed
ind = [[i for i in range(j,80,5)] for j in range(5)]

print "Reading file"
f = open(infile,'r')
lines = [f.readline().rstrip('\n') for i in range(no_reads)]
f.close()

print "Constructing dictionaries"
d = [{} for i in range(5)]
for i in range(no_reads):
	l = [[lines[i][j] for j in ind[k]] for k in range(5)]
	s = [''.join(l[k]) for k in range(5)]
	for k in range(5):
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
		l = [[lines[current][j+i] for j in ind[k]] for k in range(5)]
		s = [''.join(l[k]) for k in range(5)]
		inter = set()
		for k in range(5):	
			if s[k] in d[k]:
				inter = set(d[k][s[k]]).union(inter)
		inter = inter.intersection(remainingreads)
		if len(inter)>0:
			for j in inter:
				if(hamming(lines[current][i:],lines[j][:readlen-i])<=thresh):
					current = j
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
