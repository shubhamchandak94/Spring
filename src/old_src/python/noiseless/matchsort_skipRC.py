#Similar to matchsort_skip.py but can handle RC.

#Reordering for noiseless reads with RC. Instead of using first matchlen bases as dictionary, it uses the bases from 
#rangestart to rangeend. Also while comparing current read to reads in the bin, the bases from skipstart-readlen are not
#considered.

#The idea was to try to utilise the fact that the middle part of real reads is relatively noiseless and hence maybe we can
#use the noiseless reordering by ignoring the last few bases. However the middle part is not completely noiseless and this 
#does not work well.


import random

infile = "SRR959239.dna"
outfile = "temp3.dna"
readlen = 98
no_reads = 5372832
rangestart = 20
rangeend = 40
skipstart = 80
maxmatch = 18

print "Reading file"
f = open(infile,'r')
lines = [f.readline().rstrip('\n') for i in range(no_reads)]
f.close()

print "Constructing dictionary"
d = {}
for i in range(no_reads):
	if lines[i][rangestart:rangeend] in d:
		d[lines[i][rangestart:rangeend]].add(i)
	else:
		d[lines[i][rangestart:rangeend]] = set([i])

print "Ordering reads and writing to file"
remainingreads = set([i for i in range(no_reads)])
current = 1000000
fout = open(outfile,'w')
while True:
	flag = 0
	if len(remainingreads)%1000000 == 0:
		print str(len(remainingreads)//1000000)+'M reads remain'
	fout.write(lines[current]+'\n')
	d[lines[current][rangestart:rangeend]].remove(current)
	remainingreads.remove(current)
	if(len(remainingreads)==0):
		break
	if len(d[lines[current][rangestart:rangeend]]) == 0:
		del d[lines[current][rangestart:rangeend]]
	else:
		for i in d[lines[current][rangestart:rangeend]]:
			if lines[current][rangeend:skipstart] == lines[i][rangeend:skipstart]:
				current = i
				flag = 1
				break
	if flag == 1:
		continue
	for j in range(1,maxmatch):
		if lines[current][j+rangestart:j+rangeend] in d:
			for i in d[lines[current][j+rangestart:j+rangeend]]:
				if lines[current][j+rangeend:skipstart] == lines[i][rangeend:skipstart-j]:
					current = i
					flag = 1
					break
			if flag == 1:
					break
	if flag == 1:
		continue
#	current = random.sample(remainingreads,1)[0]
	current = remainingreads.pop()
	remainingreads.add(current)

print "Done"
fout.close()
	
		

