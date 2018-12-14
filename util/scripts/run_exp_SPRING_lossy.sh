#!/bin/bash
file_1=/raid/shubham/HARC/data/NA12878-Rep-1_S1_S2_R1.fastq
file_2=/raid/shubham/HARC/data/NA12878-Rep-1_S1_S2_R2.fastq
logname=logs/8_29_18/NA12878-Rep-1_S1_S2.log
outputname_lossy=/raid/shubham/HARC/data/NA12878-Rep-1_S1_S2_lossy.spring
stdbuf -oL /usr/bin/time -v ./build/spring -c -i $file_1 $file_2 -o $outputname_lossy -r --no-ids |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./build/spring -d -i $outputname_lossy -o $file_1".lossy.d" $file_2".lossy.d"  |& tee -a $logname
rm $file_1".lossy.d" $file_2".lossy.d"

file_1=/raid/shubham/HARC/data/NA12878-Rep-1_S1_L001_R1_001.fastq
file_2=/raid/shubham/HARC/data/NA12878-Rep-1_S1_L001_R2_001.fastq
logname=logs/8_29_18/NA12878-Rep-1_S1_L001.log
outputname_lossy=/raid/shubham/HARC/data/NA12878-Rep-1_S1_L001_lossy.spring
stdbuf -oL /usr/bin/time -v ./build/spring -c -i $file_1 $file_2 -o $outputname_lossy -r --no-ids |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./build/spring -d -i $outputname_lossy -o $file_1".lossy.d" $file_2".lossy.d"  |& tee -a $logname
rm $file_1".lossy.d" $file_2".lossy.d"

file_1=/raid/shubham/HARC/data/SRR554369_1.fastq
file_2=/raid/shubham/HARC/data/SRR554369_2.fastq
logname=logs/8_29_18/SRR554369.log
outputname_lossy=/raid/shubham/HARC/data/SRR554369_lossy.spring
stdbuf -oL /usr/bin/time -v ./build/spring -c -i $file_1 $file_2 -o $outputname_lossy -r --no-ids -q ill_bin |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./build/spring -d -i $outputname_lossy -o $file_1".lossy.d" $file_2".lossy.d"  |& tee -a $logname
rm $file_1".lossy.d" $file_2".lossy.d"

file_1=/raid/shubham/HARC/data/SRR327342_1.fastq
file_2=/raid/shubham/HARC/data/SRR327342_2.fastq
logname=logs/8_29_18/SRR327342.log
outputname_lossy=/raid/shubham/HARC/data/SRR327342_lossy.spring
stdbuf -oL /usr/bin/time -v ./build/spring -c -i $file_1 $file_2 -o $outputname_lossy -r --no-ids -q ill_bin |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./build/spring -d -i $outputname_lossy -o $file_1".lossy.d" $file_2".lossy.d"  |& tee -a $logname
rm $file_1".lossy.d" $file_2".lossy.d"

file_1=/raid/shubham/HARC/data/ERR532393_1.fastq
file_2=/raid/shubham/HARC/data/ERR532393_2.fastq
logname=logs/8_29_18/ERR532393.log
outputname_lossy=/raid/shubham/HARC/data/ERR532393_lossy.spring
stdbuf -oL /usr/bin/time -v ./build/spring -c -i $file_1 $file_2 -o $outputname_lossy -r --no-ids -q ill_bin |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./build/spring -d -i $outputname_lossy -o $file_1".lossy.d" $file_2".lossy.d"  |& tee -a $logname
rm $file_1".lossy.d" $file_2".lossy.d"

file_1=/raid/shubham/HARC/data/SRR870667_2.fastq
logname=logs/8_29_18/SRR870667_2.log
outputname_lossy=/raid/shubham/HARC/data/SRR870667_2_lossy.spring
stdbuf -oL /usr/bin/time -v ./build/spring -c -i $file_1 -o $outputname_lossy -r --no-ids -q ill_bin |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./build/spring -d -i $outputname_lossy -o $file_1".lossy.d" |& tee -a $logname
rm $file_1".lossy.d"

file_1=/raid/shubham/HARC/data/NA12878-Rep-1_S1_S2_R1.fastq.trimmed
file_2=/raid/shubham/HARC/data/NA12878-Rep-1_S1_S2_R2.fastq.trimmed
logname=logs/8_29_18/NA12878-Rep-1_S1_S2.trimmed.log
outputname_lossy=/raid/shubham/HARC/data/NA12878-Rep-1_S1_S2_trimmed_lossy.spring
stdbuf -oL /usr/bin/time -v ./build/spring -c -i $file_1 $file_2 -o $outputname_lossy -r --no-ids |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./build/spring -d -i $outputname_lossy -o $file_1".lossy.d" $file_2".lossy.d"  |& tee -a $logname
rm $file_1".lossy.d" $file_2".lossy.d"

file_1=/raid/shubham/HARC/data/NA12878-Rep-1_S1_L001_R1_001.fastq.trimmed
file_2=/raid/shubham/HARC/data/NA12878-Rep-1_S1_L001_R2_001.fastq.trimmed
logname=logs/8_29_18/NA12878-Rep-1_S1_L001.trimmed.log
outputname_lossy=/raid/shubham/HARC/data/NA12878-Rep-1_S1_L001_trimmed_lossy.spring
stdbuf -oL /usr/bin/time -v ./build/spring -c -i $file_1 $file_2 -o $outputname_lossy -r --no-ids |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./build/spring -d -i $outputname_lossy -o $file_1".lossy.d" $file_2".lossy.d"  |& tee -a $logname
rm $file_1".lossy.d" $file_2".lossy.d"

file_1=/raid/shubham/HARC/data/HS_28_1.fastq
file_2=/raid/shubham/HARC/data/HS_28_2.fastq
logname=logs/8_29_18/HS_28.log
outputname_lossy=/raid/shubham/HARC/data/HS_28_lossy.spring
stdbuf -oL /usr/bin/time -v ./build/spring -c -i $file_1 $file_2 -o $outputname_lossy -r --no-ids -q ill_bin |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./build/spring -d -i $outputname_lossy -o $file_1".lossy.d" $file_2".lossy.d"  |& tee -a $logname
rm $file_1".lossy.d" $file_2".lossy.d"

file_1=/raid/shubham/HARC/data/SRR1770413_1.fastq
file_2=/raid/shubham/HARC/data/SRR1770413_2.fastq
logname=logs/8_29_18/SRR1770413.log
outputname_lossy=/raid/shubham/HARC/data/SRR1770413_lossy.spring
stdbuf -oL /usr/bin/time -v ./build/spring -c -i $file_1 $file_2 -o $outputname_lossy -r --no-ids -q ill_bin |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./build/spring -d -i $outputname_lossy -o $file_1".lossy.d" $file_2".lossy.d"  |& tee -a $logname
rm $file_1".lossy.d" $file_2".lossy.d"

file_1=/raid/shubham/HARC/data/SRR062634_1.fastq
file_2=/raid/shubham/HARC/data/SRR062634_2.fastq
logname=logs/8_29_18/SRR062634.log
outputname_lossy=/raid/shubham/HARC/data/SRR062634_lossy.spring
stdbuf -oL /usr/bin/time -v ./build/spring -c -i $file_1 $file_2 -o $outputname_lossy -r --no-ids -q ill_bin |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./build/spring -d -i $outputname_lossy -o $file_1".lossy.d" $file_2".lossy.d"  |& tee -a $logname
rm $file_1".lossy.d" $file_2".lossy.d"

file_1=/raid/shubham/HARC/data/PhiX_100M_1.fastq
file_2=/raid/shubham/HARC/data/PhiX_100M_2.fastq
logname=logs/8_29_18/PhiX_100M.log
outputname_lossy=/raid/shubham/HARC/data/PhiX_100M_lossy.spring
stdbuf -oL /usr/bin/time -v ./build/spring -c -i $file_1 $file_2 -o $outputname_lossy -r --no-ids |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./build/spring -d -i $outputname_lossy -o $file_1".lossy.d" $file_2".lossy.d"  |& tee -a $logname
rm $file_1".lossy.d" $file_2".lossy.d"
