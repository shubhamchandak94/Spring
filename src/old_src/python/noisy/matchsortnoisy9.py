#Similar to matchsortnoisy8.py but tries to find best match for the given shift, if the best match is less than thresh then 
#pick it or else go ahead.

#The idea was to reduce the noise, but the decrease is very marginal.

from itertools import imap
import operator

# hamming2 function from http://code.activestate.com/recipes/499304-hamming-distance/
def hamming2(str1, str2):
    #assert len(str1) == len(str2)
    #ne = str.__ne__  ## this is surprisingly slow
    ne = operator.ne
    return sum(imap(ne, str1, str2))
#	return sum([str1[i]!=str2[i] for i in range(len(str1))])	

def char2index(c):
	if c == 'A':
		return 0		
	if c == 'C':
		return 1
	if c == 'G':
		return 2
	if c == 'T':
		return 3
	if c == 'N':
		return 4
def findmajority(count):
	l = []
	for i in range(len(count[0])):
		s = [count[j][i] for j in range(5)]
		maxcount = max(s[0:4])
		if maxcount == 0: #only N's seen so far
			l.append('A')
			continue
		if s[0] == maxcount:
			l.append('A')
			continue
		if s[1] == maxcount:
			l.append('C')
			continue
		if s[2] == maxcount:
			l.append('G')
			continue
		if s[3] == maxcount:
			l.append('T')
			continue
	return ''.join(l)

infile = "chrom22_50x_noRC_noisy.dna"
outfile = "tempte8.dna"
readlen = 100
no_reads = 17500000
matchlen = 80
maxmatch = 20
thresh = 4 #maximum number of mismatches allowed
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
count = [[0]*100 for j in range(5)]
for j in range(readlen):
	count[char2index(lines[current][j])][j] = 1
ref = findmajority(count)
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
		l = [[ref[j+i] for j in ind[k]] for k in range(5)]
		s = [''.join(l[k]) for k in range(5)]
		bestmatch = 100
		inter = set()
		for k in range(5):	
			if s[k] in d[k]:
				inter = inter.union(set(d[k][s[k]]))
		inter = inter.intersection(remainingreads)
		if len(inter)>0:
			for j in inter:
				if(hamming2(ref[i:],lines[j][:readlen-i])<bestmatch):
					bestmatch = hamming2(ref[i:],lines[j][:readlen-i])
					bestmatchindex = j
			if bestmatch <= thresh:
				current = bestmatchindex
				flag = 1
		if flag == 1:
			for j in range(readlen-i):
				count[char2index(lines[current][j])][i+j] += 1		
			count = [count[j][i:]+[0]*i for j in range(5)]
			for j in range(readlen-i,readlen):
				count[char2index(lines[current][j])][j] = 1
			ref = findmajority(count)	
			break

#if len(d[lines[current][0:matchlen]]) == 0:
#		del d[lines[current][0:matchlen]]
#	else:
#		for i in d[lines[current][0:matchlen]]:
#			if hamming2(lines[current][matchlen:],lines[i][matchlen:]) <= thresh:
#				current = i
#				flag = 1
#				break
#	if flag == 1:
#		continue
#	for j in range(1,maxmatch):
#		if lines[current][j:j+matchlen] in d:
#			for i in d[lines[current][j:j+matchlen]]:
#				if hamming2(lines[current][j+matchlen:],lines[i][matchlen:readlen-j]) <= thresh:
#					current = i
#					flag = 1
#					break
#			if flag == 1:
#					break
	if flag == 1:
		continue
	current = remainingreads.pop()
	remainingreads.add(current)
	count = [[0]*100 for j in range(5)]
	for j in range(readlen):
		count[char2index(lines[current][j])][j] = 1
	ref = findmajority(count)


print "Done"
fout.close()
