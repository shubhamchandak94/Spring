#!/bin/bash
sync && echo 3 > /proc/sys/vm/drop_caches

logname=logs/8_29_18/additional_experiments.log

file_1=/raid/shubham/HARC/data/NA12878-Rep-1_S1_L001_R1_001.fastq.trimmed
file_2=/raid/shubham/HARC/data/NA12878-Rep-1_S1_L001_R2_001.fastq.trimmed

outputname_lossless=/raid/shubham/HARC/data/NA12878-Rep-1_S1_L001_trimmed_r_id.spring
stdbuf -oL /usr/bin/time -v ./build/spring -c -i $file_1 $file_2 -r --no-quality -o $outputname_lossless |& tee -a $logname
sync && echo 3 > /proc/sys/vm/drop_caches

outputname_lossless=/raid/shubham/HARC/data/NA12878-Rep-1_S1_L001_trimmed_1d.spring
stdbuf -oL /usr/bin/time -v ./build/spring_1d -c -i $file_1 $file_2 -o $outputname_lossless |& tee -a $logname
sync && echo 3 > /proc/sys/vm/drop_caches

file_1=/raid/shubham/HARC/data/NA12878-Rep-1_S1_S2_R1.fastq.trimmed
file_2=/raid/shubham/HARC/data/NA12878-Rep-1_S1_S2_R2.fastq.trimmed
outputname_lossless=/raid/shubham/HARC/data/NA12878-Rep-1_S1_S2_trimmed_no_early_stop.spring
stdbuf -oL /usr/bin/time -v ./build/spring_no_early_stop -c -i $file_1 $file_2 -o $outputname_lossless |& tee -a $logname
sync && echo 3 > /proc/sys/vm/drop_caches
