#!/bin/bash
sync && echo 3 > /proc/sys/vm/drop_caches

logname=logs/8_29_18/additional_experiments.log

file_1=/raid/shubham/HARC/data/NA12878-Rep-1_S1_L001_R1_001.fastq.trimmed
file_2=/raid/shubham/HARC/data/NA12878-Rep-1_S1_L001_R2_001.fastq.trimmed

outputname_lossless=/raid/shubham/HARC/data/NA12878-Rep-1_S1_L001_trimmed.sam
/usr/bin/time -v ./../bwa-0.7.17/bwa mem -t 8 /raid/shubham/human_g1k_v37.fasta $file_1 $file_2 > $outputname_lossless 

outputname_lossless=/raid/shubham/HARC/data/NA12878-Rep-1_S1_L001_trimmed_r_id.spring
stdbuf -oL /usr/bin/time -v ./build/spring -c -i $file_1 $file_2 -r --no-quality -o $outputname_lossless |& tee -a $logname
sync && echo 3 > /proc/sys/vm/drop_caches

outputname_lossless=/raid/shubham/HARC/data/NA12878-Rep-1_S1_L001_trimmed_128.spring
stdbuf -oL /usr/bin/time -v ./build/spring_128 -c -i $file_1 $file_2 -o $outputname_lossless |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./build/spring_128 -d -i $outputname_lossless -o $file_1".lossless.d" $file_2".lossless.d" |& tee -a $logname
rm $file_1".lossless.d" $file_2".lossless.d"
sync && echo 3 > /proc/sys/vm/drop_caches

outputname_lossless=/raid/shubham/HARC/data/NA12878-Rep-1_S1_L001_trimmed_512.spring
stdbuf -oL /usr/bin/time -v ./build/spring_512 -c -i $file_1 $file_2 -o $outputname_lossless |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./build/spring_512 -d -i $outputname_lossless -o $file_1".lossless.d" $file_2".lossless.d" |& tee -a $logname
rm $file_1".lossless.d" $file_2".lossless.d"
sync && echo 3 > /proc/sys/vm/drop_caches
