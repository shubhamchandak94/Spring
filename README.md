# SPRING

![Build Status](https://travis-ci.org/shubhamchandak94/Spring.svg?branch=develop)

SPRING is a compression tool for Fastq files (containing up to 4.29 Billion reads):
- Near-optimal compression ratios for single-end and paired-end datasets
- Fast and memory-efficient decompression
- Supports variable length short reads of length upto 511 bases (without -l flag)
- Supports variable length long reads of arbitrary length (upto 4.29 Billion) (with -l flag)
- Supports lossless compression of reads, quality scores and read identifiers
- Supports reordering of reads (while preserving read pairing information) to boost compression
- Supports quantization of quality values using Illumina 8-level binning.
- Tested on Linux and macOS

### Download
```bash
git clone -b develop https://github.com/shubhamchandak94/SPRING.git
```
```diff
- TODO: remove -b develop once this moves to master.
```

### Install
The instructions below will create the spring executable in the build directory inside SPRING.

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
wget https://cmake.org/files/v3.12/cmake-3.12.1.tar.gz
tar -xzf cmake-3.12.1.tar.gz
cd cmake-3.12.1
./configure
make
cd ..
./cmake-3.12.1/bin/cmake ..
make
```

On macOS, install GCC compiler since Clang has issues with OpenMP library:
- Install HomeBrew (https://brew.sh/)
- Install GCC (this step will be faster if XCode command line tools are already installed using ```xcode-select --install```):
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
  -i [ --input-file ] arg         input file name (use option twice for paired 
                                  end)
  -o [ --output-file ] arg        output file name (for paired end 
                                  decompression, if only one file is specified,
                                  two output files will be created by suffixing
                                  .1 and .2.)
  -t [ --num-threads ] arg (=8)   number of threads (default 8)
  -w [ --working-dir ] arg (=.)   directory to create temporary files (default 
                                  current directory)
  -r [ --allow_read_reordering ]  do not retain read order during compression 
                                  (paired reads still remain paired).
  --no-quality                    do not retain quality values during 
                                  compression
  --no-ids                        do not retain read identifiers during 
                                  compression
  --ill-bin                       apply Illumina binning to quality scores 
                                  before compression
  -l [ --long ]                   Use for compression of arbitrarily long read 
                                  lengths. Can also provide better compression 
                                  for reads with significant number of indels. 
                                  -r disabled in this mode. For Illumina short
                                  reads, compression is better without -l flag. 
```

### Example Usage of SPRING
```diff
- TODO: update information preserving name.
```

For compressing file_1.fastq and file_2.fastq losslessly using default 8 threads (Perfectly lossless).
```bash
./spring -c -i file_1.fastq -i file_2.fastq -o outputname
```
Using 16 threads (Perfectly lossless).
```bash
./spring -c -i file_1.fastq -i file_2.fastq -o outputname -t 16
```
Compressing with only paired end info preserved, ids not stored, qualities compressed after Illumina binning (Information-preserving mode).
```bash
./spring -c -i file_1.fastq -i file_2.fastq -r --no-ids --ill-bin -o outputname
```
Compressing only reads and ids.
```bash
./spring -c -i file_1.fastq -i file_2.fastq --no-quality -o outputname
```
Compressing single-end long read Fastq losslessly.
```bash
./spring -c -l -i file.fastq  -o outputname
```
For single end file, compressing without order preserved.
```bash
./spring -c -i file.fastq -r -o outputname
```
For single end file, compressing with order preserved (lossless).
```bash
./spring -c -i file.fastq -o outputname
```
Decompressing (single end) to uncompressedfilename.
```bash
./spring -d -i compressedfilename -o uncompressedfilename
```
Decompressing (paired end) to uncompressedfilename.1 and uncompressedfilename.2.
```bash
./spring -d -i compressedfilename -o uncompressedfilename
```
Decompressing (paired end) to file_1.fastq and file_2.fastq.
```bash
./spring -d -i compressedfilename -o file_1.fastq -o file_2.fastq
```
