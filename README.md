# SPRING

[![Build Status](https://travis-ci.org/shubhamchandak94/Spring.svg?branch=master)](https://travis-ci.org/shubhamchandak94/Spring)
### [Bioinformatics publication](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty1015/5232998?guestAccessKey=266a1378-4684-4f04-bb99-6febdf9d1fb9)

SPRING is a compression tool for Fastq files (containing up to 4.29 Billion reads):
- Near-optimal compression ratios for single-end and paired-end datasets
- Fast and memory-efficient decompression
- Supports variable length short reads of length upto 511 bases (without -l flag)
- Supports variable length long reads of arbitrary length (upto 4.29 Billion) (with -l flag). This mode directly applies general purpose compression (BSC) to reads and so compression gains might be lower than those without -l flag.
- Supports lossless compression of reads, quality scores and read identifiers
- Supports reordering of reads (while preserving read pairing information) to boost compression
- Supports quantization of quality values using [QVZ](https://github.com/mikelhernaez/qvz/), [Illumina 8-level binning](https://www.illumina.com/documents/products/whitepapers/whitepaper_datacompression.pdf) and binary thresholding
- Supports decompression of a subset of reads (random access)
- Supports gzipped fastq files as input (output) during (de)compression
- Tested on Linux and macOS

**Note:** If you want to use SPRING only as a ***tool for reordering reads*** (approximately according to genome position), take a look at the [reorder-only branch](https://github.com/shubhamchandak94/Spring/tree/reorder-only).

### Download
```bash
git clone https://github.com/shubhamchandak94/SPRING.git
```

### Install
The instructions below will create the spring executable in the build directory inside SPRING. If you plan to build and run SPRING on separate architectures, then you might need to remove/comment the line ```set(FLAGS "${FLAGS} -march=native")``` in CMakeLists.txt (or use flags based on the target architecture).

On Linux with cmake installed and version at least 3.9 (check using ```cmake --version```):
```bash
cd SPRING
mkdir build
cd build
cmake ..
make
```

On Linux with cmake not installed or with version older than 3.9:
```bash
cd SPRING
mkdir build
cd build
wget https://cmake.org/files/v3.10/cmake-3.10.3.tar.gz
tar -xzf cmake-3.10.3.tar.gz
cd cmake-3.10.3
./configure
make
cd ..
./cmake-3.10.3/bin/cmake ..
make
```

On macOS, install GCC compiler since Clang has issues with OpenMP library:
- Install HomeBrew (https://brew.sh/)
- Install GCC (this step will be faster if Xcode command line tools are already installed using ```xcode-select --install```):
```bash
brew update
brew install gcc@7
```
- Set environment variables:
```bash
export CC=gcc-7
export CXX=g++-7
```
- Delete ```CMakeCache.txt``` (if present) from the build directory
- Follow the steps above for Linux

### Usage
Run the spring executable ```/PATH/TO/spring``` with the options below:
```
Allowed options:
  -h [ --help ]                   produce help message
  -c [ --compress ]               compress
  -d [ --decompress ]             decompress
  --decompress-range arg          --decompress-range start end
                                  (optional) decompress only reads (or read 
                                  pairs for PE datasets) from start to end 
                                  (both inclusive) (1 <= start <= end <= 
                                  num_reads (or num_read_pairs for PE)). If -r 
                                  was specified during compression, the range 
                                  of reads does not correspond to the original 
                                  order of reads in the FASTQ file.
  -i [ --input-file ] arg         input file name (two files for paired end)
  -o [ --output-file ] arg        output file name (for paired end 
                                  decompression, if only one file is specified,
                                  two output files will be created by suffixing
                                  .1 and .2.)
  -w [ --working-dir ] arg (=.)   directory to create temporary files (default 
                                  current directory)
  -t [ --num-threads ] arg (=8)   number of threads (default 8)
  -r [ --allow-read-reordering ]  do not retain read order during compression 
                                  (paired reads still remain paired)
  --no-quality                    do not retain quality values during 
                                  compression
  --no-ids                        do not retain read identifiers during 
                                  compression
  -q [ --quality-opts ] arg       quality mode: possible modes are
                                  1. -q lossless (default)
                                  2. -q qvz qv_ratio (QVZ lossy compression, 
                                  parameter qv_ratio roughly corresponds to 
                                  bits used per quality value)
                                  3. -q ill_bin (Illumina 8-level binning)
                                  4. -q binary thr high low (binary (2-level) 
                                  thresholding, quality binned to high if >= 
                                  thr and to low if < thr)
  -l [ --long ]                   Use for compression of arbitrarily long read 
                                  lengths. Can also provide better compression 
                                  for reads with significant number of indels. 
                                  -r disabled in this mode. For Illumina short 
                                  reads, compression is better without -l flag.
  -g [ --gzipped_fastq ]          enable if compression input is gzipped fastq
                                  or to output gzipped fastq during
                                  decompression
```
Note that the SPRING compressed files are tar archives consisting of the different compressed streams, although we recommend using the `.spring` extension as in the examples shown below.

### Example Usage of SPRING
This section contains several examples for SPRING compression and decompression with various modes and options. The compressed SPRING file uses the `.spring` extension as a convention.

For compressing file_1.fastq and file_2.fastq losslessly using default 8 threads (Lossless).
```bash
./spring -c -i file_1.fastq file_2.fastq -o file.spring
```
For compressing file_1.fastq.gz and file_2.fastq.gz (gzipped fastq files) losslessly using default 8 threads (Lossless).
```bash
./spring -c -i file_1.fastq.gz file_2.fastq.gz -o file.spring -g
```
Using 16 threads (Lossless).
```bash
./spring -c -i file_1.fastq file_2.fastq -o file.spring -t 16
```
Compressing with only paired end info preserved, ids not stored, qualities compressed after Illumina binning (Recommended lossy mode for older Illumina machines. For Novaseq files, lossless quality compression is recommmended).
```bash
./spring -c -i file_1.fastq file_2.fastq -r --no-ids -q ill_bin -o file.spring
```
Compressing with only paired end info preserved, ids not stored, qualities binary thresholded (qv < 20 binned to 6 and qv >= 20 binned to 40).
```bash
./spring -c -i file_1.fastq file_2.fastq -r --no-ids -q binary 20 40 6 -o file.spring
```
Compressing with only paired end info preserved, ids not stored, qualities quantized using qvz with approximately 1 bit used per quality value.
```bash
./spring -c -i file_1.fastq file_2.fastq -r --no-ids -q qvz 1.0 -o file.spring
```
Compressing only reads and ids.
```bash
./spring -c -i file_1.fastq file_2.fastq --no-quality -o file.spring
```
Compressing single-end long read Fastq losslessly.
```bash
./spring -c -l -i file.fastq  -o file.spring
```
For single end file, compressing without order preserved.
```bash
./spring -c -i file.fastq -r -o file.spring
```
For single end file, compressing with order preserved (lossless).
```bash
./spring -c -i file.fastq -o file.spring
```
Decompressing (single end) to file.fastq.
```bash
./spring -d -i file.spring -o file.fastq
```
Decompressing (single end) to file.fastq, only decompress reads from 400 to 10000000.
```bash
./spring -d -i file.spring -o file.fastq --decompress-range 400 1000000
```
Decompressing (paired end) to file.fastq.1 and file.fastq.2.
```bash
./spring -d -i file.spring -o file.fastq
```
Decompressing (paired end) to file_1.fastq and file_2.fastq.
```bash
./spring -d -i file.spring -o file_1.fastq file_2.fastq
```
Decompressing (paired end) to file_1.fastq.gz and file_2.fastq.gz.
```bash
./spring -d -i file.spring -o file_1.fastq.gz file_2.fastq.gz -g
```
Decompressing (paired end) to file_1.fastq and file_2.fastq, only decompress pairs from 4000000 to 8000000.
```bash
./spring -d -i file.spring -o file_1.fastq file_2.fastq --decompress-range 4000000 8000000
```
