#Reordering for noisy reads without RC. Similar to v2 in report without RC. Uses five dictionaries (ind) and a majority
#based reference.
#Note that here reads with 'N' are supported.

#Worked really well for simulated noise without RC. Slightly worse with real data due to RC.

from distance import hamming

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

infile = "SRR959239.dna"
outfile = "temp1.dna"
readlen = 98
no_reads = 5372832
matchlen = 80
maxmatch = 18
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
count = [[0]*readlen for j in range(5)]
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
		inter = set()
		for k in range(5):	
			if s[k] in d[k]:
				inter = inter.union(set(d[k][s[k]]))
		inter = inter.intersection(remainingreads)
		if len(inter)>0:
			for j in inter:
				if(hamming(ref[i:],lines[j][:readlen-i])<=thresh):
					current = j
					flag = 1
					break
		if flag == 1:
			for j in range(readlen-i):
				count[char2index(lines[current][j])][i+j] += 1		
			count = [count[j][i:]+[0]*i for j in range(5)]
			for j in range(readlen-i,readlen):
				count[char2index(lines[current][j])][j] = 1
			ref = findmajority(count)	
			break

	if flag == 1:
		continue
	current = remainingreads.pop()
	remainingreads.add(current)
	count = [[0]*readlen for j in range(5)]
	for j in range(readlen):
		count[char2index(lines[current][j])][j] = 1
	ref = findmajority(count)


print "Done"
fout.close()
