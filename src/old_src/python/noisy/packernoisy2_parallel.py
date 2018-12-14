#Same as packernoisy2 but uses multiple threads, the number of threads is given as a parameter numthreads.
#Supports reads with 'N' and hence is quite slow as compared to packernoisy_noN.py
#Divides the infile into numthreads parts which are then encoded separately, in all the output file, a blank line is placed to 
#demarcate between the outputs of different threads.
#It makes sense to keep numthreads quite large, maybe larger than the actual number of cores. This is because it allows for 
#the variable time which might be needed to encode different parts of the reordered file.


from joblib import Parallel, delayed
from distance import hamming

infile = "temp1.dna"
outfile_seq = "read_seq33.txt"
outfile_flag = "read_flag33.txt"
outfile_noise = "read_noise33.txt"
outfile_noisepos = "read_noisepos33.txt"

no_reads = 5128790
readlen = 98
maxmatch = 18
thresh = 20 # maximum number of mismatches allowed 
numthreads = 20

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

def pack(lines):
  s = ['','','','']
  ref = 'A'*readlen # ref is the reference which is constantly updated (introduced because matching a read to previous read leads to double noise than actual)
  prev = 'A'*readlen
  count = [[1]*readlen,[0]*readlen,[0]*readlen,[0]*readlen,[0]*readlen] #number of A's,C's,T's,G's and N's seen at each position in ref
  #Note: N is never considered in the ref - we arbitrarily place an A if only N's are seen at some position
  for current in lines:
      flag = 0
      for i in range(maxmatch):
        if(hamming(current[:(readlen-i)],ref[i:])<=thresh):
          if(hamming(current[:(readlen-i)],ref[i:])<=hamming(current[:(readlen-i)],prev[i:])):
            s[1]+='r'
            s[0]+=(current[(readlen-i):]+'\n')
            prevj = 0;
            for j in range(readlen-i):
              count[char2index(current[j])][i+j] += 1		
              if current[j]!=ref[i+j]:
                s[2]+=(current[j])
                s[3]+=("%02d"%(j-prevj))#delta encoding
                prevj = j	
          else:
            s[1]+=('p')
            s[0]+=(current[(readlen-i):]+'\n')
            prevj = 0;
            for j in range(readlen-i):
              count[char2index(current[j])][i+j] += 1		
              if current[j]!=prev[i+j]:
                s[2]+=(current[j])
                s[3]+=("%02d"%(j-prevj))#delta encoding
                prevj = j	

          count = [count[j][i:]+[0]*i for j in range(5)]
          for j in range(readlen-i,readlen):
            count[char2index(current[j])][j] = 1

          ref = findmajority(count)	
          #ref = current#ref[i:]+current[readlen-i:]
          s[2]+=('\n')
          flag = 1
          break

      if flag == 0:
        s[1]+=('0')
        s[0]+=(current+'\n')
        count = [[0]*readlen for j in range(5)]
        for j in range(readlen):
          count[char2index(current[j])][j] = 1
        ref = findmajority(count)
      prev = current
  return s						
 	
f_seq = open(outfile_seq,'w')
f_flag = open(outfile_flag,'w')
f_noise = open(outfile_noise,'w')
f_noisepos = open(outfile_noisepos,'w')

i = no_reads/numthread
startline = [i*j for j in range(numthreads)]
endline = [startline[j+1]-1 for j in range(numthreads-1)]
endline.append(no_reads)
f = open(infile,'r')
lines = [f.readline().rstrip('\n') for i in range(no_reads)]
s = Parallel(n_jobs=-1)(delayed(pack)(lines[startline[i]:endline[i]+1]) for i in range(numthreads))
print "Packing done, writing to files"
for i in range(numthreads):
  f_seq.write(s[i][0]+'\n')
  f_flag.write(s[i][1]+'\n')
  f_noise.write(s[i][2]+'\n')
  f_noisepos.write(s[i][3]+'\n')

f_seq.close()
f_flag.close()	
f_noise.close()
f_noisepos.close()
