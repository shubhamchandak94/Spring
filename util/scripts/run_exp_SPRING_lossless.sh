#!/bin/bash
file_1=/raid/shubham/HARC/data/NA12878-Rep-1_S1_S2_R1.fastq
file_2=/raid/shubham/HARC/data/NA12878-Rep-1_S1_S2_R2.fastq
logname=logs/8_29_18/NA12878-Rep-1_S1_S2.log
outputname_lossless=/raid/shubham/HARC/data/NA12878-Rep-1_S1_S2_lossless.spring
stdbuf -oL /usr/bin/time -v ./build/spring -c -i $file_1 $file_2 -o $outputname_lossless |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./build/spring -d -i $outputname_lossless -o $file_1".lossless.d" $file_2".lossless.d"  |& tee -a $logname
rm $file_1".lossless.d" $file_2".lossless.d"

file_1=/raid/shubham/HARC/data/NA12878-Rep-1_S1_L001_R1_001.fastq
file_2=/raid/shubham/HARC/data/NA12878-Rep-1_S1_L001_R2_001.fastq
logname=logs/8_29_18/NA12878-Rep-1_S1_L001.log
outputname_lossless=/raid/shubham/HARC/data/NA12878-Rep-1_S1_L001_lossless.spring
stdbuf -oL /usr/bin/time -v ./build/spring -c -i $file_1 $file_2 -o $outputname_lossless |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./build/spring -d -i $outputname_lossless -o $file_1".lossless.d" $file_2".lossless.d"  |& tee -a $logname
rm $file_1".lossless.d" $file_2".lossless.d"

file_1=/raid/shubham/HARC/data/SRR554369_1.fastq
file_2=/raid/shubham/HARC/data/SRR554369_2.fastq
logname=logs/8_29_18/SRR554369.log
outputname_lossless=/raid/shubham/HARC/data/SRR554369_lossless.spring
stdbuf -oL /usr/bin/time -v ./build/spring -c -i $file_1 $file_2 -o $outputname_lossless |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./build/spring -d -i $outputname_lossless -o $file_1".lossless.d" $file_2".lossless.d"  |& tee -a $logname
rm $file_1".lossless.d" $file_2".lossless.d"

file_1=/raid/shubham/HARC/data/SRR327342_1.fastq
file_2=/raid/shubham/HARC/data/SRR327342_2.fastq
logname=logs/8_29_18/SRR327342.log
outputname_lossless=/raid/shubham/HARC/data/SRR327342_lossless.spring
stdbuf -oL /usr/bin/time -v ./build/spring -c -i $file_1 $file_2 -o $outputname_lossless |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./build/spring -d -i $outputname_lossless -o $file_1".lossless.d" $file_2".lossless.d"  |& tee -a $logname
rm $file_1".lossless.d" $file_2".lossless.d"

file_1=/raid/shubham/HARC/data/ERR532393_1.fastq
file_2=/raid/shubham/HARC/data/ERR532393_2.fastq
logname=logs/8_29_18/ERR532393.log
outputname_lossless=/raid/shubham/HARC/data/ERR532393_lossless.spring
stdbuf -oL /usr/bin/time -v ./build/spring -c -i $file_1 $file_2 -o $outputname_lossless |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./build/spring -d -i $outputname_lossless -o $file_1".lossless.d" $file_2".lossless.d"  |& tee -a $logname
rm $file_1".lossless.d" $file_2".lossless.d"

file_1=/raid/shubham/HARC/data/SRR870667_2.fastq
logname=logs/8_29_18/SRR870667_2.log
outputname_lossless=/raid/shubham/HARC/data/SRR870667_2_lossless.spring
stdbuf -oL /usr/bin/time -v ./build/spring -c -i $file_1 -o $outputname_lossless |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./build/spring -d -i $outputname_lossless -o $file_1".lossless.d" |& tee -a $logname
rm $file_1".lossless.d"

file_1=/raid/shubham/HARC/data/NA12878-Rep-1_S1_S2_R1.fastq.trimmed
file_2=/raid/shubham/HARC/data/NA12878-Rep-1_S1_S2_R2.fastq.trimmed
logname=logs/8_29_18/NA12878-Rep-1_S1_S2.trimmed.log
outputname_lossless=/raid/shubham/HARC/data/NA12878-Rep-1_S1_S2_trimmed_lossless.spring
stdbuf -oL /usr/bin/time -v ./build/spring -c -i $file_1 $file_2 -o $outputname_lossless |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./build/spring -d -i $outputname_lossless -o $file_1".lossless.d" $file_2".lossless.d"  |& tee -a $logname
rm $file_1".lossless.d" $file_2".lossless.d"

file_1=/raid/shubham/HARC/data/NA12878-Rep-1_S1_L001_R1_001.fastq.trimmed
file_2=/raid/shubham/HARC/data/NA12878-Rep-1_S1_L001_R2_001.fastq.trimmed
logname=logs/8_29_18/NA12878-Rep-1_S1_L001.trimmed.log
outputname_lossless=/raid/shubham/HARC/data/NA12878-Rep-1_S1_L001_trimmed_lossless.spring
stdbuf -oL /usr/bin/time -v ./build/spring -c -i $file_1 $file_2 -o $outputname_lossless |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./build/spring -d -i $outputname_lossless -o $file_1".lossless.d" $file_2".lossless.d"  |& tee -a $logname
rm $file_1".lossless.d" $file_2".lossless.d"

file_1=/raid/shubham/HARC/data/HS_28_1.fastq
file_2=/raid/shubham/HARC/data/HS_28_2.fastq
logname=logs/8_29_18/HS_28.log
outputname_lossless=/raid/shubham/HARC/data/HS_28_lossless.spring
stdbuf -oL /usr/bin/time -v ./build/spring -c -i $file_1 $file_2 -o $outputname_lossless |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./build/spring -d -i $outputname_lossless -o $file_1".lossless.d" $file_2".lossless.d"  |& tee -a $logname
rm $file_1".lossless.d" $file_2".lossless.d"

file_1=/raid/shubham/HARC/data/SRR1770413_1.fastq
file_2=/raid/shubham/HARC/data/SRR1770413_2.fastq
logname=logs/8_29_18/SRR1770413.log
outputname_lossless=/raid/shubham/HARC/data/SRR1770413_lossless.spring
stdbuf -oL /usr/bin/time -v ./build/spring -c -i $file_1 $file_2 -o $outputname_lossless |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./build/spring -d -i $outputname_lossless -o $file_1".lossless.d" $file_2".lossless.d"  |& tee -a $logname
rm $file_1".lossless.d" $file_2".lossless.d"

file_1=/raid/shubham/HARC/data/SRR062634_1.fastq
file_2=/raid/shubham/HARC/data/SRR062634_2.fastq
logname=logs/8_29_18/SRR062634.log
outputname_lossless=/raid/shubham/HARC/data/SRR062634_lossless.spring
stdbuf -oL /usr/bin/time -v ./build/spring -c -i $file_1 $file_2 -o $outputname_lossless |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./build/spring -d -i $outputname_lossless -o $file_1".lossless.d" $file_2".lossless.d"  |& tee -a $logname
rm $file_1".lossless.d" $file_2".lossless.d"

file_1=/raid/shubham/HARC/data/SRR1284073.fastq
logname=logs/8_29_18/SRR1284073.log
outputname_lossless=/raid/shubham/HARC/data/SRR1284073_lossless.spring
stdbuf -oL /usr/bin/time -v ./build/spring -c -i $file_1 -l -o $outputname_lossless |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./build/spring -d -i $outputname_lossless -o $file_1".lossless.d" |& tee -a $logname
rm $file_1".lossless.d"

file_1=/raid/shubham/HARC/data/ERR637420.fastq
logname=logs/8_29_18/ERR637420.log
outputname_lossless=/raid/shubham/HARC/data/ERR637420_lossless.spring
stdbuf -oL /usr/bin/time -v ./build/spring -c -i $file_1 -l -o $outputname_lossless |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./build/spring -d -i $outputname_lossless -o $file_1".lossless.d"  |& tee -a $logname
rm $file_1".lossless.d"

file_1=/raid/shubham/HARC/data/PhiX_100M_1.fastq
file_2=/raid/shubham/HARC/data/PhiX_100M_2.fastq
logname=logs/8_29_18/PhiX_100M.log
outputname_lossless=/raid/shubham/HARC/data/PhiX_100M_lossless.spring
stdbuf -oL /usr/bin/time -v ./build/spring -c -i $file_1 $file_2 -o $outputname_lossless |& tee -a $logname
stdbuf -oL /usr/bin/time -v ./build/spring -d -i $outputname_lossless -o $file_1".lossless.d" $file_2".lossless.d"  |& tee -a $logname
rm $file_1".lossless.d" $file_2".lossless.d"
