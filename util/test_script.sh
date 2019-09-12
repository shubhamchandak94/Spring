#!/bin/bash

set -e

./spring -c -i ../util/test_1.fastq -o abcd 
./spring -d -i abcd -o tmp
cmp tmp ../util/test_1.fastq


./spring -c -i ../util/test_1.fastq ../util/test_2.fastq -o abcd 
./spring -d -i abcd -o tmp
cmp tmp.1 ../util/test_1.fastq
cmp tmp.2 ../util/test_2.fastq


./spring -c -i ../util/test_1.fastq -o abcd -l 
./spring -d -i abcd -o tmp
cmp tmp ../util/test_1.fastq

./spring -d -i abcd -o tmp.gz -g
gunzip -f tmp.gz
cmp tmp ../util/test_1.fastq

./spring -c -i ../util/test_1.fastq ../util/test_2.fastq -o abcd -l
./spring -d -i abcd -o tmp
cmp tmp.1 ../util/test_1.fastq
cmp tmp.2 ../util/test_2.fastq

./spring -c -i ../util/test_1.fastq.gz ../util/test_2.fastq.gz -o abcd -g
./spring -d -i abcd -o tmp
cmp tmp.1 ../util/test_1.fastq
cmp tmp.2 ../util/test_2.fastq

./spring -c -i ../util/test_1.fastq.gz -o abcd -g
./spring -d -i abcd -o tmp
cmp tmp ../util/test_1.fastq

./spring -d -i abcd -o tmp.gz -g
gunzip -f tmp.gz
cmp tmp ../util/test_1.fastq

./spring -c -i ../util/test_1.fastq.gz ../util/test_2.fastq.gz -o abcd -g
./spring -d -i abcd -o tmp
cmp tmp.1 ../util/test_1.fastq
cmp tmp.2 ../util/test_2.fastq

./spring -d -g -i abcd -o tmp.1.gz tmp.2.gz
gunzip -f tmp.1.gz
gunzip -f tmp.2.gz
cmp tmp.1 ../util/test_1.fastq
cmp tmp.2 ../util/test_2.fastq

./spring -c -i ../util/test_1.fastq -o abcd -t 8 
./spring -d -i abcd -o tmp -t 5
cmp tmp ../util/test_1.fastq

./spring -c -i ../util/test_1.fastq ../util/test_2.fastq -o abcd -t 8
./spring -d -i abcd -o tmp -t 5
cmp tmp.1 ../util/test_1.fastq
cmp tmp.2 ../util/test_2.fastq


./spring -c -i ../util/test_1.fastq -o abcd -r
./spring -d -i abcd -o tmp
sort tmp > tmp.sorted
sort ../util/test_1.fastq > tmp_1.sorted 
cmp tmp.sorted tmp_1.sorted


./spring -c -i ../util/test_1.fastq ../util/test_2.fastq -o abcd -t 8
./spring -d -i abcd -o tmp -t 5
sort tmp.1 > tmp.sorted
sort ../util/test_1.fastq > tmp_1.sorted 
cmp tmp.sorted tmp_1.sorted
sort tmp.2 > tmp.sorted
sort ../util/test_2.fastq > tmp_1.sorted 
cmp tmp.sorted tmp_1.sorted

rm abcd tmp*
