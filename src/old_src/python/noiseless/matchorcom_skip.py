#Reordering for noiseless reads without RC which have already been reordered by some algorithm.
#Combines ideas from matchorcom.py and matchsort_skip.py
#While building the chunks, the last few bases (skipstart:readlen) are not checked. For dictionary construction, bases 
#rangestart:rangeend are used instead of 0:matchlen. 

#The idea was to use the fact that the middle bases in a read are relatively noiseless and hence maybe a noiseless reordering
#can work. However it doesn't seem to work well.
import random
infile = "temp1.dna"
outfile = "temp4.dna"
readlen = 100
no_reads = 5372832
rangestart = 10
rangeend = 20
skipstart = 80
maxmatch = 30

print "Reading file"
f = open(infile,'r')
lines = [f.readline().rstrip('\n') for i in range(no_reads)]
f.close()

print "Finding orcom chunks"
prev = lines[0]
numchunks = 1
chunkpos = []
rshift = 1 # 0 for false, 1 for True, 2 for don't know (when shift is 0)
chunkstart = True
chunkstartpos = 0
for k in range(1,no_reads):
	if(k%1000000 == 0):
		print str(k//1000000)+"M reads done"
	prev = lines[k-1]
	current = lines[k]
	flag = 0
	if chunkstart == True:
		for i in range(maxmatch):
			if(current[:skipstart-i]==prev[i:skipstart]):
				flag = 1
				if i == 0:
				    	rshift = 2
				else:
					rshift = 1
				chunkstart = False
				break
		if flag == 0:
			for i in range(1,maxmatch):
				if(current[i:skipstart]==prev[:(skipstart-i)]):
					flag = 1
					rshift = 0
					chunkstart = False
					break
		if flag == 0:
			chunkpos.append([k-1,k-1])
			chunkstartpos = k
			chunkstart = True
			numchunks+=1
	else:
		if rshift == 1:		
			for i in range(maxmatch):
				if(current[:skipstart-i]==prev[i:skipstart]):
					flag = 1
					break
			if flag == 0:
				chunkpos.append([chunkstartpos,k-1])
				chunkstartpos = k
				chunkstart = True
				numchunks+=1
	
		elif rshift == 0:
			for i in range(maxmatch):
				if(current[i:skipstart]==prev[:(skipstart-i)]):
					flag = 1
					break
			if flag == 0:
				chunkpos.append([k-1,chunkstartpos])
				chunkstartpos = k
				chunkstart = True
				numchunks+=1
		else:
			for i in range(maxmatch):
				if(current[:skipstart-i]==prev[i:skipstart]):
					flag = 1
					if i > 0:
						rshift = 1
					break
			if flag == 0:
				for i in range(1,maxmatch):
					if(current[i:skipstart]==prev[:(skipstart-i)]):
						flag = 1
						rshift = 0	
						break
			if flag == 0:
				chunkpos.append([k-1,chunkstartpos])
				chunkstartpos = k
				chunkstart = True
				numchunks+=1
if chunkstart == True:
	chunkpos.append([no_reads-1,no_reads-1])
else:
	if rshift == True:
		chunkpos.append([chunkstartpos,no_reads-1])
	else:
		chunkpos.append([no_reads-1,chunkstartpos])

print "Number of chunks found = "+str(numchunks)
random.shuffle(chunkpos)

print "Constructing dictionary"
d = {}
for i in range(numchunks):
 	j = chunkpos[i][0]	
	if lines[j][rangestart:rangeend] in d:
		d[lines[j][rangestart:rangeend]].add(i)
	else:
		d[lines[j][rangestart:rangeend]] = set([i])

print "Ordering chunks and writing to file"

remainingchunks = set([i for i in range(numchunks)])
current = 0
fout = open(outfile,'w')
while True:
	flag = 0
	if len(remainingchunks)%1000000 == 0:
		print str(len(remainingchunks)//1000000)+'M chunks remain'
	beg = chunkpos[current][0]
	end = chunkpos[current][1]
	if end >= beg: #right shifting
		for i in range(beg,end+1):
			fout.write(lines[i]+'\n')
	else:
		for i in range(beg,end-1,-1):
			fout.write(lines[i]+'\n')		

	d[lines[beg][rangestart:rangeend]].remove(current)
	remainingchunks.remove(current)
	if(len(remainingchunks)==0):
		break
	if len(d[lines[beg][rangestart:rangeend]]) == 0:
		del d[lines[beg][rangestart:rangeend]]
	
	for j in range(maxmatch):
		if lines[end][j+rangestart:j+rangeend] in d:
			for i in d[lines[end][j+rangestart:j+rangeend]]:
				if lines[end][j+rangeend:skipstart] == lines[chunkpos[i][0]][rangeend:skipstart-j]:
					current = i
					flag = 1
					break
			if flag == 1:
					break
	if flag == 1:
		continue
#	current = random.sample(remainingchunks,1)[0]
	current = remainingchunks.pop()
	remainingchunks.add(current)

print "Done"
fout.close()
	
		

