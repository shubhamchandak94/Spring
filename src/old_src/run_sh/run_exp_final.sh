#!/bin/bash
acc="ERP001775_1"
stdbuf -oL /usr/bin/time -v ./../pigz-2.3.4/pigz -f -p 8 data/fastq/"$acc".dna.d |& tee -a logs/"$acc".log
stdbuf -oL /usr/bin/time -v ./../leon/leon -c -file data/fastq/"$acc".fastq -nb-cores 8 -seq-only -verbose 0 |& tee -a logs/"$acc".log
stdbuf -oL /usr/bin/time -v ./../leon/leon -d -file data/fastq/"$acc".fastq.leon -nb-cores 8 -seq-only -verbose 0 |& tee -a logs/"$acc".log
stdbuf -oL /usr/bin/time -v ./../orcom/bin/orcom_bin e -idata/fastq/"$acc".fastq -odata/fastq/"$acc".bin -t8 |& tee -a logs/"$acc".log
stdbuf -oL /usr/bin/time -v ./../orcom/bin/orcom_pack e -idata/fastq/"$acc".bin -odata/fastq/"$acc".orcom -t8 |& tee -a logs/"$acc".log
stdbuf -oL /usr/bin/time -v ./../orcom/bin/orcom_pack d -idata/fastq/"$acc".orcom -odata/fastq/"$acc"_orcom.dna -t8 |& tee -a logs/"$acc".log



for acc in ERR194146 ERR174310 SRR870667_2 SRR327342_1 SRR554369; do
stdbuf -oL echo $acc |& tee -a logs/"$acc".log
stdbuf -oL /usr/bin/time -v ./run_default.sh -c data/fastq/"$acc".fastq  -p |& tee -a logs/"$acc".log
stdbuf -oL /usr/bin/time -v ./run_default.sh -d data/fastq/"$acc".tar |& tee -a logs/"$acc".log
stdbuf -oL /usr/bin/time -v ./run_default.sh -d data/fastq/"$acc".tar -p |& tee -a logs/"$acc".log
stdbuf -oL /usr/bin/time -v ./../pigz-2.3.4/pigz -f -p 8 data/fastq/"$acc".dna.d |& tee -a logs/"$acc".log
stdbuf -oL /usr/bin/time -v ./../leon/leon -c -file data/fastq/"$acc".fastq -nb-cores 8 -seq-only -verbose 0 |& tee -a logs/"$acc".log
stdbuf -oL /usr/bin/time -v ./../leon/leon -d -file data/fastq/"$acc".fastq.leon -nb-cores 8 -seq-only -verbose 0 |& tee -a logs/"$acc".log
stdbuf -oL /usr/bin/time -v ./../orcom/bin/orcom_bin e -idata/fastq/"$acc".fastq -odata/fastq/"$acc".bin -t8 |& tee -a logs/"$acc".log
stdbuf -oL /usr/bin/time -v ./../orcom/bin/orcom_pack e -idata/fastq/"$acc".bin -odata/fastq/"$acc".orcom -t8 |& tee -a logs/"$acc".log
stdbuf -oL /usr/bin/time -v ./../orcom/bin/orcom_pack d -idata/fastq/"$acc".orcom -odata/fastq/"$acc"_orcom.dna -t8 |& tee -a logs/"$acc".log
done
