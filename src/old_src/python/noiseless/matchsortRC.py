#Reordering for noiseless reads with RC. Reads are written to outfile as they are and hence packer needs to check for RC.
#Otherwise similar to matchsort.py

import random
from Bio.Seq import Seq

infile = "temp.dna"
outfile = "temp3.dna"
readlen = 40
no_reads = 67617092
matchlen = 20
maxmatch = 20

print "Reading file"
f = open(infile,'r')
lines = [f.readline().rstrip('\n') for i in range(no_reads)]
f.close()

random.shuffle(lines)
print "Generating Reverse Complements"
revlines = [str(Seq(lines[i]).reverse_complement()) for i in range(no_reads)]

print "Constructing dictionaries"
d1 = {}
d2 = {}
for i in range(no_reads):
	if lines[i][0:matchlen] in d1:
		d1[lines[i][0:matchlen]].append(i)
	else:
		d1[lines[i][0:matchlen]] = [i]
	
	if revlines[i][0:matchlen] in d2:
		d2[revlines[i][0:matchlen]].append(i)
	else:
		d2[revlines[i][0:matchlen]] = [i]

print "Ordering reads and writing to file"
remainingreads = set([i for i in range(no_reads)])

current = 0
currentseq = lines[current]
unmatched = 0
fout = open(outfile,'w')

while True:
	flag = 0
	if len(remainingreads)%1000000 == 0:
		print str(len(remainingreads)//1000000)+'M reads remain'
	fout.write(lines[current]+'\n')
	d1[lines[current][0:matchlen]].remove(current)
	d2[revlines[current][0:matchlen]].remove(current)
	remainingreads.remove(current)
	if(len(remainingreads)==0):
		break
	if len(d1[lines[current][0:matchlen]]) == 0:
		del d1[lines[current][0:matchlen]]
	if len(d2[revlines[current][0:matchlen]]) == 0:
		del d2[revlines[current][0:matchlen]]
	for j in range(maxmatch):
		if currentseq[j:j+matchlen] in d1:
			for i in d1[currentseq[j:j+matchlen]]:
				if currentseq[j+matchlen:] == lines[i][matchlen:readlen-j]:
					current = i
					currentseq = lines[i]
					flag = 1
					break
			if flag == 1:
					break
		if currentseq[j:j+matchlen] in d2:
			for i in d2[currentseq[j:j+matchlen]]:
				if currentseq[j+matchlen:] == revlines[i][matchlen:readlen-j]:
					current = i
					currentseq = revlines[i]
					flag = 1
					break
			if flag == 1:
					break
	if flag == 1:
		continue
#	current = random.sample(remainingreads,1)[0]
	unmatched += 1
	current = remainingreads.pop()
	currentseq = lines[current]	
	remainingreads.add(current)

print "Done, unmatched reads = "+str(unmatched)
fout.close()
