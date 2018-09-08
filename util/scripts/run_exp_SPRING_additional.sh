#!/bin/bash
sync && echo 3 > /proc/sys/vm/drop_caches

logname=logs/8_29_18/additional_experiments.log

file_1=/raid/shubham/HARC/data/NA12878-Rep-1_S1_L001_R1_001.fastq.trimmed
file_2=/raid/shubham/HARC/data/NA12878-Rep-1_S1_L001_R2_001.fastq.trimmed

outputname_lossless=/raid/shubham/HARC/data/NA12878-Rep-1_S1_L001_trimmed_lossless.spring
stdbuf -oL /usr/bin/time -v ./build/spring -d -i $outputname_lossless -o $file_1".lossless.d" $file_2".lossless.d" --decompress_range 100000000 101000000 |& tee -a $logname
rm $file_1".lossless.d" $file_2".lossless.d"
sync && echo 3 > /proc/sys/vm/drop_caches
stdbuf -oL /usr/bin/time -v ./build/spring -d -i $outputname_lossless -o $file_1".lossless.d" $file_2".lossless.d" --decompress_range 100000000 110000000 |& tee -a $logname
rm $file_1".lossless.d" $file_2".lossless.d"
sync && echo 3 > /proc/sys/vm/drop_caches
stdbuf -oL /usr/bin/time -v ./build/spring -d -i $outputname_lossless -o $file_1".lossless.d" $file_2".lossless.d" --decompress_range 100000000 200000000 |& tee -a $logname
rm $file_1".lossless.d" $file_2".lossless.d"
sync && echo 3 > /proc/sys/vm/drop_caches

outputname_lossless=/raid/shubham/HARC/data/NA12878-Rep-1_S1_L001_trimmed_lossless_32.spring
stdbuf -oL /usr/bin/time -v ./build/spring -c -i $file_1 $file_2 -o $outputname_lossless -t 32 |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./build/spring -d -i $outputname_lossless -o $file_1".lossless.d" $file_2".lossless.d" -t 32 |& tee -a $logname
rm $file_1".lossless.d" $file_2".lossless.d"
sync && echo 3 > /proc/sys/vm/drop_caches

outputname_lossless=/raid/shubham/HARC/data/NA12878-Rep-1_S1_L001_trimmed_lossless_16.spring
stdbuf -oL /usr/bin/time -v ./build/spring -c -i $file_1 $file_2 -o $outputname_lossless -t 16 |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./build/spring -d -i $outputname_lossless -o $file_1".lossless.d" $file_2".lossless.d" -t 16 |& tee -a $logname
rm $file_1".lossless.d" $file_2".lossless.d"
sync && echo 3 > /proc/sys/vm/drop_caches

outputname_lossless=/raid/shubham/HARC/data/NA12878-Rep-1_S1_L001_trimmed_lossless_4.spring
stdbuf -oL /usr/bin/time -v ./build/spring -c -i $file_1 $file_2 -o $outputname_lossless -t 4 |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./build/spring -d -i $outputname_lossless -o $file_1".lossless.d" $file_2".lossless.d" -t 4 |& tee -a $logname
rm $file_1".lossless.d" $file_2".lossless.d"
sync && echo 3 > /proc/sys/vm/drop_caches

outputname_lossless=/raid/shubham/HARC/data/NA12878-Rep-1_S1_L001_trimmed_r_id.spring
stdbuf -oL /usr/bin/time -v ./build/spring -c -i $file_1 $file_2 -r --no-qualities -o $outputname_lossless |& tee -a $logname
sync && echo 3 > /proc/sys/vm/drop_caches

file_1=/raid/shubham/HARC/data/HS_28_1.fastq
file_2=/raid/shubham/HARC/data/HS_28_2.fastq
outputname_lossless=/raid/shubham/HARC/data/HS_28_qvz.spring
stdbuf -oL /usr/bin/time -v ./build/spring -c -i $file_1 $file_2 -o $outputname_lossless -q qvz 1 -t 32|& tee -a $logname
outputname_lossless=/raid/shubham/HARC/data/HS_28_binary.spring
stdbuf -oL /usr/bin/time -v ./build/spring -c -i $file_1 $file_2 -o $outputname_lossless -q binary 20 40 6 -t 32|& tee -a $logname
