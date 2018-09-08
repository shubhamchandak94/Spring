#!/bin/bash

#lossless

file_1=/raid/shubham/HARC/data/SRR554369_1.fastq
file_2=/raid/shubham/HARC/data/SRR554369_2.fastq
logname=logs/8_29_18/SRR554369.log
outputname_lossless=/raid/shubham/HARC/data/SRR554369_lossless.fastore
stdbuf -oL /usr/bin/time -v ./../FaStore/scripts/fastore_compress.sh --lossless --threads 8 --in $file_1 --pair $file_2  --out $outputname_lossless --verbose |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./../FaStore/bin/fastore_pack d -z -i$outputname_lossless -o"$file_1".lossless.d" $file_2".lossless.d"" -v -t8 |& tee -a $logname
rm $file_1".lossless.d" $file_2".lossless.d"

file_1=/raid/shubham/HARC/data/ERR532393_1.fastq
file_2=/raid/shubham/HARC/data/ERR532393_2.fastq
logname=logs/8_29_18/ERR532393.log
outputname_lossless=/raid/shubham/HARC/data/ERR532393_lossless.fastore
stdbuf -oL /usr/bin/time -v ./../FaStore/scripts/fastore_compress.sh --lossless --threads 8 --in $file_1 --pair $file_2  --out $outputname_lossless --verbose |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./../FaStore/bin/fastore_pack d -z -i$outputname_lossless -o"$file_1".lossless.d" $file_2".lossless.d"" -v -t8 |& tee -a $logname
rm $file_1".lossless.d" $file_2".lossless.d"

file_1=/raid/shubham/HARC/data/SRR870667_2.fastq
logname=logs/8_29_18/SRR870667_2.log
outputname_lossless=/raid/shubham/HARC/data/SRR870667_2_lossless.fastore
stdbuf -oL /usr/bin/time -v ./../FaStore/scripts/fastore_compress.sh --lossless --threads 8 --in $file_1  --out $outputname_lossless --verbose |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./../FaStore/bin/fastore_pack d -i$outputname_lossless -o$file_1".lossless.d" -v -t8 |& tee -a $logname
rm $file_1".lossless.d"


file_1=/raid/shubham/HARC/data/NA12878-Rep-1_S1_S2_R1.fastq.trimmed
file_2=/raid/shubham/HARC/data/NA12878-Rep-1_S1_S2_R2.fastq.trimmed
logname=logs/8_29_18/NA12878-Rep-1_S1_S2.trimmed.log
outputname_lossless=/raid/shubham/HARC/data/NA12878-Rep-1_S1_S2.trimmed_lossless.fastore
stdbuf -oL /usr/bin/time -v ./../FaStore/scripts/fastore_compress.sh --lossless --threads 8 --in $file_1 --pair $file_2  --out $outputname_lossless --verbose |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./../FaStore/bin/fastore_pack d -z -i$outputname_lossless -o"$file_1".lossless.d" $file_2".lossless.d"" -v -t8 |& tee -a $logname
rm $file_1".lossless.d" $file_2".lossless.d"

file_1=/raid/shubham/HARC/data/NA12878-Rep-1_S1_L001_R1_001.fastq.trimmed
file_2=/raid/shubham/HARC/data/NA12878-Rep-1_S1_L001_R2_001.fastq.trimmed
logname=logs/8_29_18/NA12878-Rep-1_S1_L001.trimmed.log
outputname_lossless=/raid/shubham/HARC/data/NA12878-Rep-1_S1_L001.trimmed_lossless.fastore
stdbuf -oL /usr/bin/time -v ./../FaStore/scripts/fastore_compress.sh --lossless --threads 8 --in $file_1 --pair $file_2  --out $outputname_lossless --verbose |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./../FaStore/bin/fastore_pack d -z -i$outputname_lossless -o"$file_1".lossless.d" $file_2".lossless.d"" -v -t8 |& tee -a $logname
rm $file_1".lossless.d" $file_2".lossless.d"

file_1=/raid/shubham/HARC/data/HS_28_1.fastq
file_2=/raid/shubham/HARC/data/HS_28_2.fastq
logname=logs/8_29_18/HS_28.log
outputname_lossless=/raid/shubham/HARC/data/HS_28_lossless.fastore
stdbuf -oL /usr/bin/time -v ./../FaStore/scripts/fastore_compress.sh --lossless --threads 8 --in $file_1 --pair $file_2  --out $outputname_lossless --verbose |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./../FaStore/bin/fastore_pack d -z -i$outputname_lossless -o"$file_1".lossless.d" $file_2".lossless.d"" -v -t8 |& tee -a $logname
rm $file_1".lossless.d" $file_2".lossless.d"

file_1=/raid/shubham/HARC/data/SRR062634_1.fastq
file_2=/raid/shubham/HARC/data/SRR062634_2.fastq
logname=logs/8_29_18/SRR062634.log
outputname_lossless=/raid/shubham/HARC/data/SRR062634_lossless.fastore
stdbuf -oL /usr/bin/time -v ./../FaStore/scripts/fastore_compress.sh --lossless --threads 8 --in $file_1 --pair $file_2  --out $outputname_lossless --verbose |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./../FaStore/bin/fastore_pack d -z -i$outputname_lossless -o"$file_1".lossless.d" $file_2".lossless.d"" -v -t8 |& tee -a $logname
rm $file_1".lossless.d" $file_2".lossless.d"

file_1=/raid/shubham/HARC/data/PhiX_100M_1.fastq
file_2=/raid/shubham/HARC/data/PhiX_100M_2.fastq
logname=logs/8_29_18/PhiX_100M.log
outputname_lossless=/raid/shubham/HARC/data/PhiX_100M_lossless.fastore
stdbuf -oL /usr/bin/time -v ./../FaStore/scripts/fastore_compress.sh --lossless --threads 8 --in $file_1 --pair $file_2  --out $outputname_lossless --verbose |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./../FaStore/bin/fastore_pack d -z -i$outputname_lossless -o"$file_1".lossless.d" $file_2".lossless.d"" -v -t8 |& tee -a $logname
rm $file_1".lossless.d" $file_2".lossless.d"
