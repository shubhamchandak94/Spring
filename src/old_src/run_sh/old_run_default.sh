#!/bin/bash
set -e
# Download data
# DATA_DIR="data/"


usage()
{
cat << EOF
usage: 
Compression - compresses FASTQ reads. Output written to .tar file
./run_default.sh -c PATH_TO_FASTQ [-p] [-t NUM_THREADS]
-p = Preserve order of reads
-t NUM_THREADS - default 8

Decompression - decompresses reads. Output written to .dna.d file
./run_default.sh -d PATH_TO_TAR [-p] [-t NUM_THREADS]
-p = Get reads in original order (slower). Only applicable if -p was used during compression.
-t NUM_THREADS - default 8

Help (this message)
./run_default.sh -h

EOF
exit 0
}

compress()
{
	pathname=$(dirname $filename)
	echo "*** Preprocessing ***"
	echo $filename
	./src/preprocess_final.out $filename $pathname

	readlen="$(head $pathname/input_clean.dna | wc -L)"
	if (($readlen > 256));then
		echo "Maximum read length exceeded" 
		exit 1
	fi
	echo "#define maxmatch $((readlen/2))" > src/cpp/noisy/config.h
	echo "#define thresh 4" >> src/cpp/noisy/config.h
	echo "#define numdict 2" >> src/cpp/noisy/config.h
	echo "#define maxsearch 1000" >> src/cpp/noisy/config.h
	echo "#define dict1_start $(( readlen > 100 ? readlen/2-32 : readlen/2-readlen*32/100 ))" >> src/cpp/noisy/config.h
	echo "#define dict1_end $((readlen/2-1))" >> src/cpp/noisy/config.h
	echo "#define dict2_start $((readlen/2))" >> src/cpp/noisy/config.h
	echo "#define dict2_end $(( readlen > 100 ? readlen/2-1+32 : readlen/2-1+readlen*32/100 ))" >> src/cpp/noisy/config.h

	echo "#define readlen $readlen" >> src/cpp/noisy/config.h
	echo "#define num_thr $num_thr" >> src/cpp/noisy/config.h

	g++ src/cpp/noisy/matchsort7_v14.cpp -march=native -O3 -fopenmp -lpthread -std=c++11 -o src/reorder_noisy.out
	mkdir -p $pathname/output/ 
	./src/reorder_noisy.out $pathname
	mv $pathname/input_N.dna $pathname/output/input_N.dna
	mv $pathname/input_clean.dna $pathname/output/input_clean.dna
	mv $pathname/read_order_N.bin $pathname/output/read_order_N.bin
	g++ src/cpp/noisy/encoder.cpp -march=native -O3 -fopenmp -std=c++11 -o src/encoder.out
	./src/encoder.out $pathname
	
	#remove temporary files
	rm $pathname/output/temp.dna
	rm $pathname/output/tempflag.txt
	rm $pathname/output/temppos.txt
	
	#compress and create tarball
	7z a $pathname/output/read_pos.txt.7z $pathname/output/read_pos.txt -mmt=$num_thr
	7z a $pathname/output/read_noise.txt.7z $pathname/output/read_noise.txt -mmt=$num_thr
	7z a $pathname/output/read_noisepos.txt.7z $pathname/output/read_noisepos.txt -mmt=$numt_thr
	7z a $pathname/output/input_N.dna.7z $pathname/output/input_N.dna -mmt=$numt_thr
	7z a $pathname/output/read_meta.txt.7z $pathname/output/read_meta.txt -mmt=$numt_thr
	7z a $pathname/output/read_rev.txt.7z $pathname/output/read_rev.txt -mmt=$num_thr
	./src/libbsc/bsc e $pathname/output/read_seq.txt $pathname/output/read_seq.txt.bsc -b512p -tT #-tT for single thread - uses too much memory in multi-threaded
	rm $pathname/output/*.txt $pathname/output/*.dna 
	if [[ $preserve_order == "True" ]];then
		7z a $pathname/output/read_order.bin.7z $pathname/output/read_order.bin -mmt=$num_thr
		7z a $pathname/output/read_order_N.bin.7z $pathname/output/read_order_N.bin -mmt=$num_thr
	fi
	rm $pathname/output/*.bin
	tar -cf $pathname/$(basename "$filename" .fastq).tar -C $pathname/output .
	rm -r $pathname/output/
}

decompress()
{
	echo "Decompression ..."
	pathname=$(dirname $filename)
	mkdir -p $pathname/output
	tar -xf $filename -C $pathname/output
	if [[ $preserve_order == "True" ]];then
		if [ ! -f $pathname/output/read_order.bin.7z ];then
			echo "Not compressed using -p flag. Order cannot be restored"
			usage
			exit 1
		fi
	fi
	7z e $pathname/output/read_pos.txt.7z -o$pathname/output/
	7z e $pathname/output/read_noise.txt.7z -o$pathname/output/
	7z e $pathname/output/read_noisepos.txt.7z -o$pathname/output/
	7z e $pathname/output/input_N.dna.7z -o$pathname/output/
	7z e $pathname/output/read_meta.txt.7z -o$pathname/output/
	7z e $pathname/output/read_rev.txt.7z -o$pathname/output/
	./src/libbsc/bsc d $pathname/output/read_seq.txt.bsc $pathname/output/read_seq.txt -tT
	if [[ $preserve_order == "True" ]];then
		7z e $pathname/output/read_order.bin.7z -o$pathname/output/
		7z e $pathname/output/read_order_N.bin.7z -o$pathname/output/
		./src/decoder_preserve.out $pathname
		echo "Restoring order of reads"
		sort -nk1 $pathname/output/output.tmp -o $pathname/output/out_sorted.tmp -T $pathname/output/
		cut -d " " -f 2 $pathname/output/out_sorted.tmp > $pathname/output/output_clean.dna
		./src/merge_N.out $pathname
		echo "Done!"
	else
		./src/decoder.out $pathname
	fi
	rm -r $pathname/output/
	mv $pathname/output.dna $pathname/$(basename "$filename" .tar).dna.d
}

#Initialize variables to default values.
num_thr=8

#Check the number of arguments. If none are passed, print help and exit.
NUMARGS=$#
if [ $NUMARGS -eq 0 ]; then
 usage
fi

mode=''
preserve_order="False"

while getopts ':c:d:t:p' opt; do
  case "$opt" in
    c) [[ -n "$mode" ]] && usage || mode='c' && filename=$OPTARG;;
    d) [[ -n "$mode" ]] && usage || mode='d' && filename=$OPTARG;;
    t) num_thr=$OPTARG;;
    p) preserve_order="True";;
    h) usage ;;
    \?) usage ;;
    *) usage ;;
  esac
done

if [[ $mode == 'c' ]];then
compress
elif [[ $mode == 'd' ]];then
decompress
else
echo "Either -c or -d required"
usage
exit 1
fi;
