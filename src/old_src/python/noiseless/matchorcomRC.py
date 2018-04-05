#Similar to matchorcom.py but with support for RC.

#Reordering for noiseless reads with RC which have already been reordered by some algorithm.
#It initially tries to find chunks (note that we say orcom chunks but any other algorithm will also work), where the chunk can
#have the reads in the right shift or the left shift - the finding chunks part identifies this with the first few reads in the 
#chunk and then searches only for shift in that direction. For each chunk the index of the first and last read in the chunk is 
#stored (chunk is reversed if the shift is leftward).
#After the chunks are identified, the first read in each chunk is indexed by the first matchlen bases and the program tries to 
#reorder the chunks.

#The idea behind this was to reduce the memory consumption for the dictionary since only the first read of each chunk needs to
#put in a dictionary.

#Running this on orcom ordered reads (for noiseless RC reads) produces very good results - slightly better than matchsortRC.py.
#The orcom chunks contained around 10 reads/chunk on average.
#Also due to random shuffling of chunks, running it again and again reduces the number of hard reads up to some low value.

#Unfortunately this idea does not extend well to noisy reads because some reads are left out from inside the chunks and hence
#reordering the chunk may not help. Also the chunks produced in the noisy case have only2-3 reads on average so the memory 
#saving is also insignificant.

infile = "temp1.dna"
outfile = "temp2.dna"
readlen = 34
no_reads = 5372832
matchlen = 10
maxmatch = 10

from Bio.Seq import Seq
import random

print "Reading file"
f = open(infile,'r')
lines = [f.readline().rstrip('\n') for i in range(no_reads)]
f.close()

print "Generating Reverse Complements"
revlines = [str(Seq(lines[i]).reverse_complement()) for i in range(no_reads)]

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
	revcurrent = revlines[k]
	flag = 0
	if chunkstart == True:
		for i in range(maxmatch):
			if(current[:readlen-i]==prev[i:] or revcurrent[:readlen-i]==prev[i:]):
				flag = 1
				if i == 0:
				    	rshift = 2
				else:
					rshift = 1
				chunkstart = False
				break
		if flag == 0:
			for i in range(1,maxmatch):
				if(current[i:]==prev[:(readlen-i)] or revcurrent[i:]==prev[:readlen-i]):
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
				if(current[:readlen-i]==prev[i:] or revcurrent[:readlen-i]==prev[i:]):
					flag = 1
					break
			if flag == 0:
				chunkpos.append([chunkstartpos,k-1])
				chunkstartpos = k
				chunkstart = True
				numchunks+=1
	
		elif rshift == 0:
			for i in range(maxmatch):
				if(current[i:]==prev[:(readlen-i)] or revcurrent[i:]==prev[:readlen-i]):
					flag = 1
					break
			if flag == 0:
				chunkpos.append([k-1,chunkstartpos])
				chunkstartpos = k
				chunkstart = True
				numchunks+=1
		else:
			for i in range(maxmatch):
				if(current[:readlen-i]==prev[i:] or revcurrent[:readlen-i]==prev[i:]):
					flag = 1
					if i > 0:
						rshift = 1
					break
			if flag == 0:
				for i in range(1,maxmatch):
					if(current[i:]==prev[:(readlen-i)] or revcurrent[i:]==prev[:readlen-i]):
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
d1 = {}
d2 = {}
for i in range(numchunks):
 	j = chunkpos[i][0]	
	if lines[j][0:matchlen] in d1:
		d1[lines[j][0:matchlen]].add(i)
	else:
		d1[lines[j][0:matchlen]] = set([i])
	if revlines[j][0:matchlen] in d2:
		d2[revlines[j][0:matchlen]].add(i)
	else:
		d2[revlines[j][0:matchlen]] = set([i]) 

print "Ordering chunks and writing to file"
unmatched = 1
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

	d1[lines[beg][0:matchlen]].remove(current)
	d2[revlines[beg][0:matchlen]].remove(current)
	remainingchunks.remove(current)
	if(len(remainingchunks)==0):
		break
	if len(d1[lines[beg][0:matchlen]]) == 0:
		del d1[lines[beg][0:matchlen]]
	if len(d2[revlines[beg][0:matchlen]]) == 0:
		del d2[revlines[beg][0:matchlen]]
	
	for j in range(maxmatch):
		if lines[end][j:j+matchlen] in d1:
			for i in d1[lines[end][j:j+matchlen]]:
				if lines[end][j+matchlen:] == lines[chunkpos[i][0]][matchlen:readlen-j]:
					current = i
					flag = 1
					break
			if flag == 1:
					break
		if lines[end][j:j+matchlen] in d2:
			for i in d2[lines[end][j:j+matchlen]]:
				if lines[end][j+matchlen:] == revlines[chunkpos[i][0]][matchlen:readlen-j]:
					current = i
					flag = 1
					break
			if flag == 1:
					break
	if flag == 1:
		continue
#	current = random.sample(remainingreads,1)[0]
	current = remainingchunks.pop()
	remainingchunks.add(current)
	unmatched += 1
print "Done, unmatched chunks = "+str(unmatched)
fout.close()
