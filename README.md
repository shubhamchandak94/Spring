Exactly one of compress or decompress needs to be specified 
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
  -r [ --allow_read_reordering ]  do not retain read order during compression 
                                  (paired reads still remain paired).
  --no-quality                    do not retain quality values during 
                                  compression
  --no-ids                        do not retain read identifiers during 
                                  compression
  -w [ --working-dir ] arg (=.)   directory to create temporary files (default 
                                  current directory)
  --ill-bin                       apply Illumina binning to quality scores 
                                  before compression
  -l [ --long ]                   Use for compression of arbitrarily long read 
                                  lengths. Can also provide better compression 
                                  for reads with significant number of indels. 
                                  Some other options might be disabled in this 
                                  mode.

