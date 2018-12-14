#!/bin/bash

set -e


# path to binaries
#
FASTORE_BIN=/raid/shubham/FaStore/bin/fastore_bin
FASTORE_REBIN=/raid/shubham/FaStore/bin/fastore_rebin
FASTORE_PACK=/raid/shubham/FaStore/bin/fastore_pack
#
#

log()
{
	if [ ! -z ${VERBOSE+x} ]; then
		printf "$1\n"
	fi
}

print_usage()
{
	echo ""
	echo "    FaStore -- a space saving solution for raw sequencing data"
	echo ""
	echo "usage: bash $0 <mode> [--fast] [--threads <th>]"
	echo "          --in <in.fq> [--pair <pair.fq>] --out <archive>"
	echo "          [--signature] [--verbose] [--help]"
	echo ""
	echo "where:"
	echo "    <mode> specifies the compression mode which can be one of:" 
	echo "      --lossless     - lossless mode with only reads reordered" 
	echo "      --lossy_old    - applied Illumina q-scores binning and ids not preserved"
	echo "      --lossy_novaseq- ids not preserved"
	echo "      --max          - applied q-scores binary thresholding and w/o read ids"
	echo "    --fast           - C0 compression mode (C1 by default)"
	echo "    --threads <th>   - the number of processing threads"
	echo "    --in <in.fq>     - the FASTQ file to be processed. Multiple FASTQ files can"
	echo "                       be passed in form: --in \"<in01.fq> <in02.fq> ... \""
	echo "    --pair <pair.fq> - the paired FASTQ file(s). When specified, the compressor"
	echo "                       will assume that the reads are paired with the ones from"
	echo "                       the file specified by '--in'."
	echo "    --out <archive>  - the prefix of the output archive files"
	echo "    --signature <n>  - the length of the signature (default: 8)"
	echo "    --verbose        - print additional information while compressing"
	echo "    --help           - displays this message"
	echo ""
	exit 1
}


# parse input args
#
if [ $# -eq 0 ]
  then
	print_usage
	exit 1
fi

while [[ $# -gt 0 ]]
do
	ARG="$1"

	case $ARG in
		--lossless)
			MODE=0
			shift 1;;
		--lossy_old)
			MODE=1
			shift 1;;
		--lossy_novaseq)
			MODE=2
			shift 1;;
		--max)
			MODE=3
			shift 1;;
		--fast)
			FAST_MODE=1
			shift 1;;
		--in)
			IN="$2"
			shift 2;;
		--pair)
			PAIR="$2"
			PAR_PE="-z"
			shift 2;;
		--out)
			OUT_PFX="$2"
			shift 2;;
		--threads)
			THREADS="$2"
			shift 2;;
		--signature)
			SIG_LEN="$2"
			shift 2;;
		--verbose)
			VERBOSE=1
			shift 1;;
		--help|*) 
			echo "Unkown argument: \"$ARG\""
			print_usage
			exit 1;;
	esac
done


# check params
#
if [ -z ${MODE+x} ]; then
	echo "ERROR: compression mode has not been specified"
	exit 1
fi

if [ -z ${IN+x} ]; then
	echo "ERROR: no input files have been specified"
	exit 1
fi

if [ -z ${THREADS+x} ]; then
	log "WARN: number of threads mode has not been specified --> setting to 1"
	THREADS=1
fi

if [ -z ${SIG_LEN+x} ]; then
	SIG_LEN=8
fi

if [ -z ${OUT_PFX+x} ]; then
	log "WARN: no output files prefix name has been specified --> setting to OUT"
	OUT_PFX="OUT"
fi


# parse params
#
case $MODE in
	0) PAR_ID="-H"; PAR_QUA="-q0";;
	1) PAR_ID=""; PAR_QUA="-q2";;
	2) PAR_ID=""; PAR_QUA="-q0";;
	3) PAR_ID=""; PAR_QUA="-q1";;
esac


# set processing params
#
PAR_BIN_C1="-p$SIG_LEN -s0 -b256"
PAR_REBIN_C1="-r -w1024 -W1024"
PAR_PACK_C1="-r -f256 -c10 -d8 -w1024 -W1024"

PAR_BIN_C0="-p$SIG_LEN -s10 -b256"
PAR_PACK_C0="-f256 -c10 -d8 -w256 -W256"

TH_BIN=$THREADS
TH_REBIN=$THREADS
TH_PACK=$THREADS


# temporary and output files configuration
#
TMP_PFX="__tmp-dna-$RANDOM"
TMP_BIN="$TMP_PFX-bin"
TMP_REBIN="$TMP_PFX-rebin"

OUT_PACK="$OUT_PFX"


if [ -z ${PAIR+x} ]; then
	IN_FQ="$IN"
else
	IN_FQ="$IN $PAIR"
fi

if [ ! -z ${VERBOSE+x} ]; then
	PAR_PACK_VB="-v"
fi

#echo "$PAR_ID, $PAR_QUA, $PAR_PE, $IN, $PAIR, $OUT_PFX, $THREADS"



log "processing files: $IN_FQ"

if [ -z ${FAST_MODE+x} ]; then

	log "\n:: binning ..."
	$FASTORE_BIN e "-i$IN_FQ" "-o$TMP_BIN" "-t$TH_BIN" $PAR_ID $PAR_QUA $PAR_BIN_C1 $PAR_PE
	log "temporary files:"
	log "$(ls -s $TMP_BIN.*)"

	log "\n:: rebinning: 0 -> 2 ..."
	$FASTORE_REBIN e "-i$TMP_BIN" "-o$TMP_REBIN-2" "-t$TH_REBIN" $PAR_REBIN_C1 $PAR_PE -p2
	log "temporary files:"
	log "$(ls -s $TMP_REBIN-2.*)"
	rm $TMP_BIN*

	log "\n:: rebinning: 2 -> 4 ..."
	$FASTORE_REBIN e "-i$TMP_REBIN-2" "-o$TMP_REBIN-4" "-t$TH_REBIN" $PAR_REBIN_C1 $PAR_PE -p4
	log "temporary files:"
	log "$(ls -s $TMP_REBIN-4.*)"
	rm $TMP_REBIN-2*

	log "\n:: rebinning: 4 -> 8 ..."
	$FASTORE_REBIN e "-i$TMP_REBIN-4" "-o$TMP_REBIN-8" "-t$TH_REBIN" $PAR_REBIN_C1 $PAR_PE -p8
	log "temporary files:"
	log "$(ls -s $TMP_REBIN-8.*)"
	rm $TMP_REBIN-4*

	log "\n:: packing ..."
	$FASTORE_PACK e "-i$TMP_REBIN-8" "-o$OUT_PACK" "-t$TH_PACK" $PAR_PACK_C1 $PAR_PE $PAR_PACK_VB
	log "archive files:"
	log "$(ls -s $OUT_PACK.*)"
	rm $TMP_REBIN-8*

else

	log "\n:: binning ..."
	$FASTORE_BIN e "-i$IN_FQ" "-o$TMP_BIN" "-t$TH_BIN" $PAR_ID $PAR_QUA $PAR_BIN_C0 $PAR_PE
	log "temporary files:"
	log "$(ls -s $TMP_BIN.*)"

	log "\n:: packing..."
	$FASTORE_PACK e "-i$TMP_BIN" "-o$OUT_PACK" "-t$TH_PACK" $PAR_PACK_C0 $PAR_PE $PAR_PACK_VB
	log "archive files:"
	log "$(ls -s $OUT_PACK.*)"
	rm $TMP_BIN*

fi
