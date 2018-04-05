# SPRING

Compression tool for Fastq files. For paired end reads, both files should have same number of reads. Achieves near-optimal compression ratios and fast decompression. Supports upto 4.29 Billion fixed-length reads with lengths at most 256. The algorithm requires C++11 and g++ compiler and works on Linux. p7zip should already be installed ( on linux, run: sudo apt-get install p7zip-full).

### Installation
```bash
git clone https://github.com/shubhamchandak94/SPRING.git
cd SPRING
./install.sh
```

### Usage
##### Compression - compresses FASTQ files.
```bash
./spring -c -1 Fastq_file_1 [-2 Fastq_file_2] [-p] [-t num_threads] [-q mode] [-r qvz_ratio] [-i] -o outputfile
-2 second file for paired end reads
-p Preserve order of reads
-t num_threads - Default 8
-q Retain quality values. Possible modes:
	qvz - qvz specify bits/quality ratio using -r flag (default 8.0 lossless)
	bsc - use bsc compressor. 
	illumina_binning_bsc - bin into 8 levels and use bsc
	illumina_binning_qvz - bin into 8 levels and use qvz lossless
-r bits/quality ratio if -q qvz used [default 8.0 lossless]
-i Retain read IDs, if not specified fake ids will be generated during decompression
-o Output file name
```

##### Decompression - decompress compressed archive.
```bash
./spring -d compressed_file -o outputfile [-t num_threads]
-o outputfile name, if compressed with -P flag, two files created: outputfile.1 and outputfile.2. If quality is not retained, FASTA file is produced.
-t num_threads - Default 8
```

##### Help (this message)
```bash
./spring -h
```

### Example Usage of SPRING
For compressing file_1.fastq and file_2.fastq losslessly using default 8 threads and qvz for qualities.
```bash
./spring -c -1 file_1.fastq -2 file_2.fastq -q qvz -i -p -o outputname
```
Using BSC instead of QVZ, 16 threads.
```bash
./spring -c -1 file_1.fastq -2 file_2.fastq -q bsc -i -p -o outputname -t 16
```
Compressing with only paired end info preserved, ids not stored.
```bash
./spring -c -1 file_1.fastq -2 file_2.fastq -q qvz -o outputname
```
Compressing with Illumina binning followed by qvz.
```bash
./spring -c -1 file_1.fastq -2 file_2.fastq -q illumina_binning_qvz -i -p -o outputname
```
Compressing with Illumina binning followed by BSC.
```bash
./spring -c -1 file_1.fastq -2 file_2.fastq -q illumina_binning_bsc -i -p -o outputname
```
Compressing with QVZ lossy compression with rate 1.0 bits/quality value.
```bash
./spring -c -1 file_1.fastq -2 file_2.fastq -q qvz -r 1.0 -i -p -o outputname
```
Compressing only reads and ids.
```bash
./spring -c -1 file_1.fastq -2 file_2.fastq -i -p -o outputname
```
Compressing only reads losslessly.
```bash
./spring -c -1 file_1.fastq -2 file_2.fastq -p -o outputname
```
Compressing only reads with paired end info preserved.
```bash
./spring -c -1 file_1.fastq -2 file_2.fastq -o outputname
```
For single end file, compressing without order preserved.
```bash
./spring -c -1 file_1.fastq -q qvz -i -o outputname
```
For single end file, compressing with order preserved (lossless).
```bash
./spring -c -1 file_1.fastq -q qvz -i -p -o outputname
```
Decompression.
```bash
./spring -d compressedfilename -o uncompressedfilename
```
<!--
### Downloading datasets
###### Genomic sequencing reads
```bash
wget -b ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR554/SRR554369/SRR554369_1.fastq.gz
wget -b ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR554/SRR554369/SRR554369_2.fastq.gz
gunzip SRR554369_1.fastq.gz SRR554369_2.fastq.gz
cat SRR554369_1.fastq SRR554369_2.fastq > SRR554369.fastq
```

For some datasets (e.g. SRR327342 and SRR870667), the two fastq files may have reads of different lengths. To use HARC on such datasets, compress the two fastq files separately.

###### Metagenomics data
```bash
wget -b http://public.genomics.org.cn/BGI/gutmeta/High_quality_reads/MH0001/081026/MH0001_081026_clean.1.fq.gz
wget -b http://public.genomics.org.cn/BGI/gutmeta/High_quality_reads/MH0001/081026/MH0001_081026_clean.2.fq.gz
gunzip MH0001_081026_clean.1.fq.gz MH0001_081026_clean.2.fq.gz
cat MH0001_081026_clean.1.fq MH0001_081026_clean.2.fq > MH0001_081026_clean.fq
```

###### Human genome (hg19 - for generating simulated reads)
```bash
wget -b ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
```

Inside this file different chromosomes are demarcated.

##### Generating reads using gen_fastq (from orcom repo)
###### Error-free reads (without reverse complementation) - 35M reads of length 100 from chrom 22
```bash
cd util/gen_fastq_noRC
make
./gen_fastq_noRC 35000000 100 PATH/chrom22.fasta PATH/chrom22_reads.fastq
```

###### Reads with 1% uniform substitution rate (each base is equally likely to be changed to any of the 4 possibilities e.g. A - C,T,G,N) (without reverse complementation) - 35M reads of length 100 from chrom 22
```bash
./gen_fastq_noRC 35000000 100 PATH/chrom22.fasta PATH/chrom22_reads.fastq -e
```

##### Typical fastq format
```
@seq id
read
+
quality score
```

### Computing noise entropy
The directory util/ contains quality_counts.cpp and noise_entropy.py, which can be used to compute the noise entropy upper bound using the method described in the Supplementary Data (https://github.com/shubhamchandak94/HARC/blob/master/supplementary-data.pdf). To use these, 
1. Write the quality values (every fourth line of the FASTQ file) to a separate file, e.g., by using 
```bash
sed -n '4~4p' file.fastq > file.quality
```
2. Modify the read length in the header of quality_counts.cpp and compile it by running
```bash
g++ util/quality_counts.cpp -O3 -march=native -std=c++11 -o util/quality_counts.out
```
3. Generate the quality counts by running
```bash
./util/quality_counts.out file.quality file.quality_counts
```
4. Update the read length and input file in the header of noise_entropy.py and run 
```bash
python util/noise_entropy.py
``` 
to get the noise entropies for models of different orders.

### Other compressors (for evaluation)

##### Installing & Running orcom (boost should be installed)
```bash
git clone https://github.com/lrog/orcom.git
cd orcom
make boost
cd bin
./orcom_bin e -iPATH/SRR065390.fastq -oPATH/SRR065390.bin
./orcom_pack e -iPATH/SRR065390.bin -oPATH/SRR065390.orcom
```

##### Getting orcom ordered file
```bash
./orcom_pack d -iPATH/SRR065390.orcom -oPATH/SRR065390.dna
```

The dna file is of the form:
```
read1
read2
..
```

##### Installing and running Leon (seq-only mode)
```bash
git clone https://github.com/GATB/leon.git
cd leon
sh INSTALL
./leon -file PATH_TO_FASTQ -c -seq-only -nb-cores NUM_THR
```

##### Leon Decompression (order preserved)
```bash
./leon -file PATH_TO_LEON -d -nb-cores NUM_THR
```

The fasta.d file is of the form:
```
>1
read1
>2
read2
..
```
-->
<!---

##### Some python packages
```
sudo pip install distance
sudo pip install biopython
sudo pip install joblib
```

##### Running Proposed tool (noiseless)
First set the parameters at the top of the files.
```
g++ reordernoiseless.cpp -std=c++11 -o a.out
./a.out
python encodernoiseless.py
xz -k outfiles
```

##### Running Proposed tool (noisy)
First set the parameters at the top of the files.
```
g++ reordernoisy.cpp -std=c++11 -o a.out
./a.out
python encodernoisy.py
xz -k outfiles
```

##### Decompressing (noisy)
First set the parameters at the top of the files.
```
python decodernoisy.py
```
##### Using Google sparsehashmap
###### Installing
```
git clone https://github.com/sparsehash/sparsehash.git
cd sparsehash
./configure
sudo make install
```

###### Using
In the code replace 
```cpp
#include<unordered_map>
```
by
```cpp
#include<sparsehash/sparse_hash_map>
```
and replace the data type
```cpp
std::unordered_map<>
```
by
```cpp
google::sparse_hash_map<>
```
There should be no need to change anything else. However note that sparsehashmap does not seem to work with bitset as the index. One way to get around this (if the index length is not too large) is to use the bitset::to_ulong function.
-->

<!--

##### Running Python code
```
python matchsortnoisy8.py
python packernoisy2.py
```
##### Calculating Number of Hard reads
```
grep 0 -o read_flag.txt | wc -l
```

##### Converting Fastq file to dna file (needed for our code)
```
sed -n '2~4p' SRR065390.fastq > SRR065390.dna &
```



##### Removing reads containing 'N' (python code)
```python
fout = open('SRR065390_clean.dna','w')
with open('SRR065390.dna','r') as f:
  for line in f:
     if 'N' not in line:
          fout.write(line)
```


##### Converting dna file to fastq file (with fake seq id and quality scores) (python code)
```python
fout = open('SRR065390_clean.fastq','w')
with open('SRR065390_clean.dna','r') as f:
  for line in f:
    fout.write('@\n'+line+'+\n'+line)
```

##### Calculating Number of Singleton reads (i.e. 00's in the flag file) 
```
tr -cs 0 '\012' < read_flag.txt | awk '/00/{n += length - 1}; END {print n+0}'
```
##### Generating reads using gen_fastq (from orcom repo)
###### Error-free reads (without reverse complementation) - 35M reads of length 100 from chrom 22
```
cd gen_fastq_noRC
make
./gen_fastq_noRC 35000000 100 PATH/chrom22.fasta PATH/chrom22_reads.fastq
```

###### Reads with 1% uniform substitution rate (each base is equally likely to be changed to any of the 4 possibilities e.g. A - C,T,G,N) (without reverse complementation) - 35M reads of length 100 from chrom 22
```
./gen_fastq_noRC 35000000 100 PATH/chrom22.fasta PATH/chrom22_reads.fastq -e
```


##### C++ files
These are the more important files in the C++ folder. For the other files, see comments on top of those files. Note that the noisy files are currently unable to handle reads with N. Also read length is assumed to be constant for all codes.
###### Noiseless and no RC
1. matchsort2.cpp - reordering, parameters - matchlen, maxmatch. Generates outfile with reordered reads.

###### Noisy
1. matchsort3.cpp - v1 reordering described in the report. Tries to find mathces to the current read (no clean reference). Parameters - numdict, dictionary indices, maxmatch, thresh. Generates three files: outfile which has the reordered reads (some of which reverse complemented), outfileRC which has flags (0/1) to tell if the read has been reverse complemented, outfileflag which has the flags (0/1) to tell if the read is matched or not. This information helps packernoisy2_noN.py and packernoisy4_noN.py
2. matchsort6.cpp - has a recovery step after the reordering which tries to place the singleton reads before a matching read. thresh2 is the threshold for this second stage process. Other parameters and output files are same as matchsort3.
3. matchsort7.cpp - v2 reordering described in the report. Uses majority-based reference read for the reordering. Parameters and output files are same as matchsort3

##### Python files
These are the more important files in the python folder. For the other files, see comments on top of those files. Note that the some of the noisy files can handle reads with N. Also read length is assumed to be constant for all codes.
###### Noiseless and no RC
1. matchsort.py - reordering, parameters - matchlen, maxmatch. Generates outfile with reordered reads.
2. packer.py - encoding, takes infile with reordered reads and generates outfile_seq and outfile_flag. Parameter - maxmatch, neg (whether we want to look for matches with opposite shift as well.

###### Noisy for real data
1. matchpacknoisyRC.py - v2 reordering described in the report along with encoding. Uses majority-based reference read for the reordering. Parameters - numdict, dictionary indices, maxmatch, thresh. Directly generates encoding in the form of five files - outfile_seq, outfile_flag, outfile_rev, outfile_noise and outfile_noisepos. Works with reads containing N.
3. packernoisy2.py - encoding for reordered reads as described in report (uses reference). Input - infile containing reordered reads. Parameters - thresh, maxmatch. The reordering should generate the reverse complement flags separately, this does not consider RC. Produces 4 files - outfile_seq, outfile_flag, outfile_noise and outfile_noisepos. The current implementation is slow (due to findmajority function), see packernoisy2_noN.py for faster implementation.
4. packernoisy2_parallel.py - Parallel implementation of packernoisy2.py using joblib library - see comments on top of file. Leaves a blank line in each output file after each thread. Extra parameter - numthreads. Still slower than packernoisy2_noN.py.
5. packernoisy4.py - Similar to packernoisy2.py except to encoding of noise. Instead of storing the noisy base, we store cyclic shift using 1,2,3 or 4 (note that N is also a possibility). See function encodenoise for the exact encoding.

###### Encoders for C++ generated reordering (no reads with N)
1. packernoisy2_noN.py - Similar to packernoisy2.py but can't handle reads with N. Also an extra input file is needed - infile_flag which contains the flags (0/1 - unmatched/matched) produced by the C++ codes (matchsort3,6,7.cpp). Much faster than packernoisy2.py due to the flag file and the better written findmajority function.
2. packernoisy4_noN.py - Similar to packernoisy4.py but can't handle reads with N. Also an extra input file is needed - infile_flag which contains the flags (0/1 - unmatched/matched) produced by the C++ codes (matchsort3,6,7.cpp). Much faster than packernoisy4.py due to the flag file and the better written findmajority function.

###### Reordering for noisy simulated data without RC
1. matchsortnoisy2.py - Similar to v1 reordering described in the report. Tries to find matches to the current read (no clean reference). 4 dictionaries. Parameters - matchlen, maxmatch, thresh. Produces outfile with reordered reads. Works with reads containing N.
2. matchsortnoisy8.py - Similar to v2 reordering described in the report. Tries to find matches to the clean reference. 5 dictionaries. Parameters - matchlen, maxmatch, thresh. Produces outfile with reordered reads. Works with reads containing N.
-->
