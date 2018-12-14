#Reorders reads from infile and writes reordered reads to outfile. For noiseless reads without RC. 
#Picks the first read found in the bin (does not check if next read matches shifted version of the current read.) 
#One needs to use a noisy packer to encode the reordered reads
#Depending on genome size and repeats, this can give same or worse compression than matchsort.py

infile = "SRR959239.dna"
outfile = "temp2.dna"
readlen = 98
no_reads = 5372832
matchlen = 80
maxmatch = 18

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
			if True:#lines[current][matchlen:] == lines[i][matchlen:]:
				current = i
				flag = 1
				break
	if flag == 1:
		continue
	for j in range(1,maxmatch):
		if lines[current][j:j+matchlen] in d:
			for i in d[lines[current][j:j+matchlen]]:
				if True:#lines[current][j+matchlen:] == lines[i][matchlen:readlen-j]:
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
	
		

