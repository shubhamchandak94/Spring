# qvz

Quality Value Zip (qvz) is a lossy compression algorithm for storing quality values associated
with DNA sequencing. This software implements the qvz algorithm for both encoding and decoding.

## Installing

qvz can be used on Windows, Linux, or Mac. Currently we only provide a source distribution. qvz
has no external dependencies, linking only against libc, libm, and librt.

The distribution is configured out of the box for linux. To build on a mac, copy src/Makefile.apple to
replace src/Makefile. Build with `make` in the toplevel folder. You are responsible for installing the
binary in an appropriate system-wide location (i.e. /usr/bin) if you wish.

There is currently no makefile for Windows. A visual studio project can be made (and testing has been
done on windows to verify compatibility) but the steps are beyond the scope of this guide, because you
must take several additional steps to guarantee a sane environment on windows (such as replacing stdint
and inttypes with correct versions). For this reason we also currently do not distribute a windows build
script.

## Usage

qvz is used from the command line. The general invocation is:

```qvz (options) [input file] [output file]```

Input and output must be files, currently it does not process standard input or output. The input file
must be a file consisting only of quality scores, with one read per line. Thus, the input would consist
of every fourth line in a FASTQ file. The other three lines must be compressed separately.

Available options are:

```
Operating Mode:
-q            Compress the quality score input file (default on)
-x            Extract quality values from input file

Compression Parameters:
-f [ratio]    Compress using a variable allocation of [ratio] bits per bit of input entropy per symbol
-r [rate]     Compress using a fixed allocation of [rate] bits per symbol
-d [M|L|A]    Compress while optimizing for MSE, Log(1+L1), or L1 distortions, respectively (default: MSE)

Clustering Parameters:
-c [#]        Compress using # clusters. Going above 5 is not recommended due to computational complexity (default: 1)
-T [#]        Use # as a threshold for cluster centroid movement distance before declaring an approximate clustering as "good enough"

Extra Options:
-h            Print help summary
-v            Enable verbose progress output
-s            Print summary stats to STDOUT after compression (independent of -v)
-u [file]     Write the quantized but not compressed values of [file] (default: off)
```

## Algorithm

qvz uses the approach described in a paper submitted to Bioinformatics to perform lossy compression. Data is clustered
to reduce global variability, then each cluster is compressed by calculating a set of quantization
matrices that performs optimally under the chosen distortion metric  and the empirical statistics of
the data, using a first order Markov prediction model.

## License
qvz is available under the terms of the GPLv3. See COPYING for more information.

## Bugs and Feedback
Please use GitHub issues to open bug reports or provide feedback about qvz.

## Authors
qvz was created by Greg Malysa, Mikel Hernaez, Idoia Ochoa, Milind Rao, and Karthik Ganesan at
Stanford University.
