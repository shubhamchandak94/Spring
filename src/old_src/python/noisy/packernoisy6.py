#find best match with ref or prev - for files with no Ns
from distance import hamming

char2index = {'A':0,'C':1,'G':2,'T':3,'N':4}
index2char = {0:'A',1:'C',2:'G',3:'T',4:'N'}

def findmajority(count):
	maxcount = [max(s) for s in count]
	l = [index2char[s.index(maxcount[i])] for i,s in zip(range(len(count)),count)]
	return ''.join(l)


infile = "temp5.dna"
infile_flag = "tempflag5.txt"
outfile_seq = "read_seq140.txt"
outfile_pos = "read_pos140.txt"
outfile_noise = "read_noise140.txt"
outfile_noisepos = "read_noisepos140.txt"

readlen = 100
maxmatch = 30
thresh = 30 # maximum number of mismatches allowed 
inttoascii = {0:'a',1:'b',2:'c',3:'d',4:'e',5:'f',6:'g',7:'h',8:'i',9:'j',10:'k',11:'l',12:'m',13:'n',14:'o',15:'p',16:'q',17:'r',18:'s',19:'t',20:'u',21:'w',22:'x',23:'y',24:'z',25:'A',26:'B',27:'C',28:'D',29:'E',30:'F',31:'G',32:'H',33:'I',34:'J',35:'K',36:'L',37:'M',38:'N',39:'O',40:'P',readlen:'v'}


def buildcontig(reads):
	if(len(reads) == 1): #singleton read
		return [reads[0],[0]]
	count = [[0,0,0,0,0] for i in range(readlen)] #number of A's,C's,T's,G's,N's seen at each position in ref
	pos = [0]
	for i in range(readlen):
		count[i][char2index[reads[0][i]]] = 1
	prevread = reads[0]
	for currentread in reads[1:]:
		flag = 0
		bestmatch = readlen
		besti = 0
		for i in range(maxmatch):
			hammingdist = hamming(currentread[:(readlen-i)],prevread[i:])
			if(hammingdist<=thresh):
				pos.append(i+pos[-1])
				count = count + [[0,0,0,0,0] for j in range(i)]
				for j in range(readlen):
					count[pos[-1]+j][char2index[currentread[j]]] += 1
				flag = 1
				break
			if(hammingdist < bestmatch):
				bestmatch = hammingdist
				bestmatchpos = i
		if flag == 0: #no match found due to some reason (this might happen because of matchsort9's
		# handling of N's (if only N's are seen at a position, matchsort9 makes it A in the ref.
			pos.append(bestmatchpos+pos[-1])
			count = count + [[0,0,0,0,0] for j in range(bestmatchpos)]
			for j in range(readlen):
				count[pos[-1]+j][char2index[currentread[j]]] += 1
		prevread = currentread	
	ref = findmajority(count)
	return [ref,pos]

def writecontig(ref,pos,reads):
	f_seq.write(ref+'\n')
	if pos[0]!=0:
		print "BAD"
	if len(reads) == 1: #singleton read
		f_noise.write('\n')
		f_pos.write(inttoascii[pos[0]]+'\n')
		return
	prevpos = 0
	for currentpos,currentread in zip(pos,reads):
		prevj = 0;	
		for j in range(readlen):
			if currentread[j]!=ref[currentpos+j]:
				f_noise.write(currentread[j])
				f_noisepos.write("%02d"%(j-prevj))#delta encoding
				prevj = j
		f_noise.write('\n')
		f_pos.write(inttoascii[(currentpos-prevpos)])
		prevpos = currentpos
	f_pos.write('\n')
	return		

in_flag = open(infile_flag,'r')
f_seq = open(outfile_seq,'w')
f_pos = open(outfile_pos,'w')
f_noise = open(outfile_noise,'w')
f_noisepos = open(outfile_noisepos,'w')
k = 0

reads = []

with open(infile,'r') as f:
	for line in f:
		k = k + 1
		if k%1000000 == 0:
			print str(k//1000000)+'M done'
		current = line.rstrip('\n')
		c = in_flag.read(1)
		if c=='0':
			if len(reads) != 0:
				[ref,pos] = buildcontig(reads)
				writecontig(ref,pos,reads)
			reads = [current]
		else:
			reads.append(current)

#last contig
[ref,pos] = buildcontig(reads)
writecontig(ref,pos,reads)
f_seq.close()
f_pos.close()	
f_noise.close()
f_noisepos.close()

