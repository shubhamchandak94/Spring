# SPRING

[![Build Status](https://travis-ci.org/shubhamchandak94/Spring.svg?branch=reorder-only)](https://travis-ci.org/shubhamchandak94/Spring)
### [Bioinformatics publication](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty1015/5232998?guestAccessKey=266a1378-4684-4f04-bb99-6febdf9d1fb9)

### reorder-only branch
Useful if you are only interested in obtaining the SPRING reordered FASTQ file (instead of the compressed file) when using the -r/--allow-read-reordering mode for SPRING. This can be useful for various applications since SPRING attempts to reorder reads according to their position in the genome. This produces a similar effect as compression (with -r flag) followed by decompression, but is faster. Works with single and paired datasets with at most 4.29 billion short reads of length up to 511 bases. Also supports gzipped input/output. The quality values and read identifiers are preserved and are reordered along with the read sequences. For paired data, the two output files have the paired reads in the same position.

### Download
```bash
git clone -b reorder-only https://github.com/shubhamchandak94/SPRING.git
```

### Install
The instructions below will create the spring-reorder executable in the build directory inside SPRING. If you plan to build and run SPRING on separate architectures, then you might need to remove/comment the line ```set(FLAGS "${FLAGS} -march=native")``` in CMakeLists.txt (or use flags based on the target architecture).

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
brew install gcc@9
```
- Set environment variables:
```bash
export CC=gcc-9
export CXX=g++-9
```
- Delete ```CMakeCache.txt``` (if present) from the build directory
- Follow the steps above for Linux

### Usage
Run the spring-reorder executable ```/PATH/TO/spring-reorder``` with the options below:
```
Allowed options:
  -h [ --help ]                 produce help message
  -i [ --input-file ] arg       input FASTQ file name (two files for paired
                                end)
  -o [ --output-file ] arg      output FASTQ file name (for paired end, if only
                                one file is specified, two output files will be
                                created by suffixing .1 and .2.)
  -w [ --working-dir ] arg (=.) directory to create temporary files (default
                                current directory)
  -t [ --num-threads ] arg (=8) number of threads (default 8)
  --gzipped-input               enable if compression input is gzipped fastq
  --gzipped-output              enable to output gzipped fastq
```
Note that several of the compression related options are not relevant here.

### Resource usage
Note that SPRING uses some temporary disk space, and can fail if the disk space is not sufficient. The additional temporary disk usage in this mode is around 10-20% of the original uncompressed file. The memory usage in this mode is similar to the memory usage for short-read compression in the usual Spring mode (master branch).

### Example Usage
For reordering file_1.fastq and file_2.fastq to file_1.reordered.fastq and file_2.reordered.fastq using default 8 threads.
```bash
./spring-reorder -c -i file_1.fastq file_2.fastq -o file_1.reordered.fastq file_2.reordered.fastq
```

For reordering file_1.fastq to file_1.reordered.fastq using 16 threads.
```bash
./spring-reorder -c -i file_1.fastq -o file_1.reordered.fastq -t 16
```

Gzipped input file:
```bash
./spring-reorder -c -i file_1.fastq.gz -o file_1.reordered.fastq --gzipped-input
```

Produce gzipped output FASTQ files:
```bash
./spring-reorder -c -i file_1.fastq file_2.fastq -o file_1.reordered.fastq.gz file_2.reordered.fastq.gz --gzipped-output
```
