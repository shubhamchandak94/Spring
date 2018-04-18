#!/bin/bash

# Create dirs
mkdir -p data

#libbsc for read_seq compression
(cd src/libbsc && make)
cp src/libbsc/bsc bin/

#Compilation of files
g++ src/preprocess.cpp -O3 -march=native -fopenmp -std=c++11 -o bin/preprocess.out
g++ src/pack_order.cpp -O3 -march=native -fopenmp -std=c++11 -o bin/pack_order.out
g++ src/pe_encode.cpp -O3 -march=native -fopenmp -std=c++11 -o bin/pe_encode.out
g++ src/unpack_order.cpp -O3 -march=native -fopenmp -std=c++11 -o bin/unpack_order.out
g++ src/reorder_compress_quality_id_pe.cpp src/ID_compression/src/*.c src/qvz/src/*.c -O3 -march=native -fopenmp -Isrc/ID_compression/include -Isrc/qvz/include -DLINUX -std=c++11 -o bin/reorder_compress_quality_id_pe.out -lrt
g++ src/decompress_quality_id_pe.cpp src/ID_compression/src/*.c src/qvz/src/*.c -O3 -march=native -fopenmp -Isrc/ID_compression/include -Isrc/qvz/include -std=c++11 -DLINUX -o bin/decompress_quality_id_pe.out -lrt
g++ src/decompress_quality_id_se.cpp src/ID_compression/src/*.c src/qvz/src/*.c -O3 -march=native -fopenmp -Isrc/ID_compression/include -Isrc/qvz/include -std=c++11 -DLINUX -o bin/decompress_quality_id_se.out -lrt
g++ src/pe_decode.cpp -O3 -march=native -fopenmp -std=c++11 -o bin/pe_decode.out

#The files below use bitset which needs size at compile time, so we compile for multiple sizes here
mkdir -p bin/reorder
for bitset_size in 64 128 192 256 320 384 448 512
do
	max_read_len=$(($bitset_size/2))
	echo "#define MAX_READ_LEN $max_read_len" > src/config.h
	g++ src/reorder.cpp -march=native -O3 -fopenmp -lpthread -std=c++11 -o bin/reorder/reorder_$bitset_size".out"
done

mkdir -p bin/encoder
mkdir -p bin/decoder_se
mkdir -p bin/decoder_pe
for bitset_size in 64 128 192 256 320 384 448 512 576 640 704 768
do
	max_read_len=$(($bitset_size/3))
	echo "#define MAX_READ_LEN $max_read_len" > src/config.h
	g++ src/encoder.cpp -march=native -O3 -fopenmp -lpthread -std=c++11 -o bin/encoder/encoder_$bitset_size".out"
	g++ src/decoder_pe.cpp -march=native -O3 -fopenmp -lpthread -std=c++11 -o bin/decoder_pe/decoder_pe_$bitset_size".out"
	g++ src/decoder_se.cpp -march=native -O3 -fopenmp -lpthread -std=c++11 -o bin/decoder_se/decoder_se_$bitset_size".out"
done
rm src/config.h
