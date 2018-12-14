#!/bin/bash

# Download data
# DATA_DIR="data/"
source config.py


usage()
{
cat << EOF
usage: $0 options


OPTIONS:
   -h      Show this message
   -f      download relevant files
   -p      preprocess
   -g      generateConfig
   -c      compress
   -d      Decompress
   -e      Compute Entropy Bound       
 		
EOF
}


# TODO: Make this more general later
download()
{
   
	echo "*** FASTQ Sequences being downloaded ***"
	mkdir -p data/$basename
	mkdir -p data/$basename/output
	if [[ -z "${URL_2}" ]]; then
		wget -O data/$basename/file_1.fastq.gz $URL_1
		gunzip data/$basename/file_1.fastq.gz 
		mv data/$basename/file_1.fastq data/$basename/input.fastq
	else
		echo "Combining 2 FASTQ files"
		wget -O data/$basename/file_1.fastq.gz $URL_1
		wget -O data/$basename/file_2.fastq.gz $URL_2
		gunzip data/$basename/file_1.fastq.gz data/$basename/file_2.fastq.gz
		cat data/$basename/file_1.fastq data/$basename/file_2.fastq > data/$basename/input.fastq
	fi
}

# Removes quality values and N values from the dna and writes reads with N to a separate file
preprocess()
{
	echo "*** Preprocessing ***"
	./src/preprocess_final.out data/$basename/input.fastq data/$basename
	#sort -o data/$basename/input_N.dna data/$basename/input_N.dna -T.
}

generateConfig()
{
	echo "#define maxmatch $maxmatch" > src/cpp/noisy/config.h
	echo "#define thresh $thresh" >> src/cpp/noisy/config.h
	echo "#define numdict $numdict" >> src/cpp/noisy/config.h
	echo "#define maxsearch $maxsearch" >> src/cpp/noisy/config.h
	if [  ${dict1start+x} ]; then echo "#define dict1_start $dict1start" >> src/cpp/noisy/config.h; fi
	if [  ${dict1end+x} ]; then echo "#define dict1_end $dict1end" >> src/cpp/noisy/config.h; fi
	if [  ${dict2start+x} ]; then echo "#define dict2_start $dict2start" >> src/cpp/noisy/config.h; fi
	if [  ${dict2end+x} ]; then echo "#define dict2_end $dict2end" >> src/cpp/noisy/config.h; fi
	if [  ${dict3start+x} ]; then echo "#define dict3_start $dict3start" >> src/cpp/noisy/config.h; fi
	if [  ${dict3end+x} ]; then echo "#define dict3_end $dict3end" >> src/cpp/noisy/config.h; fi
	if [  ${dict4start+x} ]; then echo "#define dict4_start $dict4start" >> src/cpp/noisy/config.h; fi
	if [  ${dict4end+x} ]; then echo "#define dict4_end $dict4end" >> src/cpp/noisy/config.h; fi
	#readlen="$(wc -L < data/$basename/input_clean.dna)"
	readlen="$(head data/$basename/input_clean.dna | wc -L)"
	echo "#define readlen $readlen" >> src/cpp/noisy/config.h
	echo "#define num_thr $num_thr" >> src/cpp/noisy/config.h
}
compress()
{
	g++ src/cpp/noisy/matchsort7_v14.cpp -march=native -O3 -fopenmp -lpthread -std=c++11 -o src/reorder_noisy.out
	mkdir -p data/$basename/output 
	./src/reorder_noisy.out data/$basename
	cp data/$basename/input_N.dna data/$basename/output/input_N.dna
	g++ src/cpp/noisy/encoder.cpp -march=native -O3 -fopenmp -std=c++11 -o src/encoder.out
	./src/encoder.out data/$basename
	
	#remove temporary files
	rm data/$basename/output/temp.dna
	rm data/$basename/output/tempflag.txt
	rm data/$basename/output/temppos.txt
       	
	#compress and create tarball
	7z a data/$basename/output/read_pos.txt.7z data/$basename/output/read_pos.txt -mmt=$num_thr
	7z a data/$basename/output/read_noise.txt.7z data/$basename/output/read_noise.txt -mmt=$num_thr
	7z a data/$basename/output/read_noisepos.txt.7z data/$basename/output/read_noisepos.txt -mmt=$numt_thr
	7z a data/$basename/output/input_N.dna.7z data/$basename/output/input_N.dna -mmt=$numt_thr
	7z a data/$basename/output/read_meta.txt.7z data/$basename/output/read_meta.txt -mmt=$numt_thr
	7z a data/$basename/output/read_rev.txt.7z data/$basename/output/read_rev.txt -mmt=$num_thr
	./src/libbsc/bsc e data/$basename/output/read_seq.txt  data/$basename/output/read_seq.txt.bsc -b512p -tT #-tT for single thread - uses too much memory in multi-threaded
	mv data/$basename/output/read_order.bin data/$basename/
	rm data/$basename/output/*.txt data/$basename/output/*.dna
	tar -cf data/$basename/output.tar data/$basename/output
	rm -r data/$basename/output/
}

decompress()
{
	echo "Decompression ..."
	ls -all data/$basename/output.tar
	tar -xf data/$basename/output.tar
	7z e data/$basename/output/read_pos.txt.7z -odata/$basename/output/
	7z e data/$basename/output/read_noise.txt.7z -odata/$basename/output/
	7z e data/$basename/output/read_noisepos.txt.7z -odata/$basename/output/
	7z e data/$basename/output/input_N.dna.7z -odata/$basename/output/
	7z e data/$basename/output/read_meta.txt.7z -odata/$basename/output/
	7z e data/$basename/output/read_rev.txt.7z -odata/$basename/output/
	./src/libbsc/bsc d data/$basename/output/read_seq.txt.bsc data/$basename/output/read_seq.txt -tT
	./src/decoder.out data/$basename
}

compute_entropy()
{
	# Install MFCompress
	wget  http://sweet.ua.pt/ap/software/mfcompress/MFCompress-linux64-1.01.tgz
	tar -zxvf MFCompress-linux64-1.01.tgz
	mkdir -p util/MFCompress
	mv MFCompress-linux64-1.01/* util/MFCompress
	rm -r MFCompress-linux64-1.01*
	pip install --user tqdm

	# MFCompress permissions
	chmod 741 ./util/MFCompress/MFCompressC
    mkdir -p logs
	echo "Downloading the FASTA File"
	wget -O data/$basename/genome_fasta.fa.gz $URL_genome
	gunzip data/$basename/genome_fasta.fa.gz
	
    echo "Computing FASTA File entropy"
    ./util/MFCompress/MFCompressC -3 data/$basename/genome_fasta.fa

	echo "computing Noise entropy"
        cp config.py logs/"$basename"_entropy_computation.log
	python util/compute_entropy.py data/$basename/input.quality data/$basename/genome_fasta.fa.mfc data/$basename/genome_fasta.fa data/$basename/output.tar | tee -a logs/"$basename"_entropy_computation.log

}

run_others()
{
	/usr/bin/time -v ./../leon/leon -c -file data/$basename/input.fastq -nb-cores $num_thr -seq-only -verbose 0 #produces input.fastq.leon
	ls -all data/$basename/input.fastq.leon
	/usr/bin/time -v ./../leon/leon -d -file data/$basename/input.fastq.leon -nb-cores $num_thr -seq-only -verbose 0 #produces input.fasta.d

	/usr/bin/time -v ./../orcom/bin/orcom_bin e -idata/$basename/input.fastq -odata/$basename/input.bin -t$num_thr
	/usr/bin/time -v ./../orcom/bin/orcom_pack e -idata/$basename/input.bin -odata/$basename/input.orcom -t$num_thr
	ls -all data/$basename/input.orcom.cdna
	/usr/bin/time -v ./../orcom/bin/orcom_pack d -idata/$basename/input.orcom -odata/$basename/input_orcom.dna -t$num_thr

	/usr/bin/time -v ./../pigz-2.3.4/pigz -f -p $num_thr data/$basename/input.dna
	ls -all data/$basename/input.dna.gz
	/usr/bin/time -v ./../pigz-2.3.4/unpigz -f -p $num_thr data/$basename/input.dna.gz
}

#Process the arguments
while getopts hfpcdgeo opt
do
   case "$opt" in
	h) usage; exit 1;;
	f) download;; 
	p) preprocess;;
	g) generateConfig;;
	c) compress;;
	d) decompress;;
	e) compute_entropy;;
	o) run_others;;
	?) usage; exit;;
   esac
done
