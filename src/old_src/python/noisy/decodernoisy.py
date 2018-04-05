#Decode the reordered reads. readlen is the length of read (assumed constant). the infile's are the 5 encoded files and
#outfile is the file to which the reordered reads will be written.

from Bio.Seq import Seq
import sys
import os
#sys.path.append("..")
#import config.ini

basename = sys.argv[1]
basedir = os.path.join(basename,"output") 
outfile = os.path.join(basename,"output.dna")
infile_meta = os.path.join(basedir,"read_meta.txt")
infile_seq = os.path.join(basedir,"read_seq.txt")
infile_pos = os.path.join(basedir,"read_pos.txt")
infile_noise = os.path.join(basedir,"read_noise.txt")
infile_noisepos = os.path.join(basedir,"read_noisepos.txt")
infile_rev = os.path.join(basedir,"read_rev.txt")
infile_N = os.path.join(basedir,"input_N.dna")

f_meta = open(infile_meta,'r')
f_seq = open(infile_seq,'r')
f_pos = open(infile_pos,'r')
f_noise = open(infile_noise,'r')
f_noisepos = open(infile_noisepos,'r')
f_rev = open(infile_rev,'r')
f_out = open(outfile,'w')

readlen = 0 
for line in f_meta:
	current = line.rstrip('\n')
	readlen = int(current)
	break
print("readlen: ", readlen)

#asciitoint = {'a':0,'b':1,'c':2,'d':3,'e':4,'f':5,'g':6,'h':7,'i':8,'j':9,'k':10,'l':11,'m':12,'n':13,'o':14,'p':15,'q':16,'r':17,'s':18,'t':19,'u':20,'w':21,'x':22,'y':23,'z':24,'A':25,'B':26,'C':27,'D':28,'E':29,'F':30,'G':31,'H':32,'I':33,'J':34,'K':35,'L':36,'M':37,'N':38,'O':39,'P':40,'v':readlen}

dec_noise = {
('A','0'):'C',('A','1'):'G',('A','2'):'T',('A','3'):'N',
('C','0'):'A',('C','1'):'G',('C','2'):'T',('C','3'):'N',
('G','0'):'T',('G','1'):'A',('G','2'):'C',('G','3'):'N',
('T','0'):'G',('T','1'):'C',('T','2'):'A',('T','3'):'N',
('N','0'):'A',('N','1'):'G',('N','2'):'C',('N','3'):'T',
}

ref = ''
while True:
	p = f_pos.read(1)
	if not p:
		break
	p = ord(p)
	ref = ref[p:]+f_seq.read(p)
	currentread = ref
	noise = f_noise.readline().rstrip('\n')
	prevnoisepos = 0
	for n in noise:
		noisepos = ord(f_noisepos.read(1))+prevnoisepos
		currentread = currentread[:noisepos]+dec_noise[(ref[noisepos],n)]+currentread[noisepos+1:]
		prevnoisepos = noisepos
	rev = f_rev.read(1)
	if rev == 'd':
		f_out.write(currentread+'\n')	
	else:
		f_out.write(str(Seq(currentread).reverse_complement())+'\n')

with open(infile_N,'r') as f_N:
	for line in f_N:
		f_out.write(line)

f_out.close()
f_pos.close()	
f_noise.close()
f_noisepos.close()
f_rev.close()
