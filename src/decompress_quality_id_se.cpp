#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <omp.h>
#include "sam_block.h"
#include "codebook.h"
#include "qv_compressor.h"
#include "cluster.h"

std::string basedir;
std::string infile_id;
std::string outfile_id;
std::string infile_quality;
std::string outfile_quality;
std::string infile_readlength;
std::string infile_order;
std::string infilenumreads;

int max_readlen, num_thr, num_thr_e;
uint32_t numreads;
std::string preserve_order, quality_mode, preserve_quality, preserve_id;

void decompress_id();
void decompress_quality(uint8_t *readlengths);
void load_readlengths(uint8_t *readlengths);

void decode(char *input_file, char *output_file, struct qv_options_t *opts, uint8_t *read_lengths);

int main(int argc, char** argv)
{
	basedir = std::string(argv[1]);
	infile_id = basedir + "/compressed_id_1.bin";
	outfile_id = basedir + "/id_1.txt";
	infile_quality = basedir + "/compressed_quality_1.bin";
	outfile_quality = basedir + "/quality_1.txt";
	infile_readlength = basedir + "/read_lengths.bin";
	infile_order = basedir + "/read_order.bin";
	infilenumreads = basedir + "/numreads.bin";
	max_readlen = atoi(argv[2]);
	num_thr = atoi(argv[3]);
	num_thr_e = atoi(argv[4]);
	preserve_order = std::string(argv[5]);
	preserve_quality = std::string(argv[6]);
	preserve_id = std::string(argv[7]);
	if(preserve_quality == "True")
		quality_mode = std::string(argv[8]);
	
	if(preserve_quality == "False" && preserve_id == "False")
		return 0;

	std::ifstream f_numreads(infilenumreads, std::ios::binary);
	f_numreads.seekg(4);
	f_numreads.read((char*)&numreads,sizeof(uint32_t));
	f_numreads.close();
 	std::cout << "Decompressing Quality and/or IDs\n";	
	omp_set_num_threads(num_thr);

	if(preserve_quality == "True" && !(quality_mode == "bsc" || quality_mode == "illumina_binning_bsc"))
	{
		uint8_t *readlengths = new uint8_t[numreads];
		load_readlengths(readlengths);
		decompress_quality(readlengths);
	}
	if(preserve_id == "True")
	{
		decompress_id();
	}
}

void decompress_id()
{
	#pragma omp parallel
	{
	int tid = omp_get_thread_num();
	for(int tid_e = tid*num_thr_e/num_thr; tid_e < (tid+1)*num_thr_e/num_thr; tid_e++)
	{
		uint32_t numreads_thr = numreads/num_thr_e;
		if(tid_e == num_thr_e - 1)
			numreads_thr = numreads-numreads_thr*(num_thr_e-1); 
		struct compressor_info_t comp_info;
		comp_info.numreads = numreads_thr;
		comp_info.mode = DECOMPRESSION;
		comp_info.fcomp = fopen((infile_id+"."+std::to_string(tid_e)).c_str(),"r");
		comp_info.f_id = fopen((outfile_id+"."+std::to_string(tid_e)).c_str(),"w");
		decompress((void *)&comp_info);
		fclose(comp_info.fcomp);
		fclose(comp_info.f_id);
	}
	}	
	return;
}

void decompress_quality(uint8_t *readlengths)
{
	#pragma omp parallel
	{
	int tid = omp_get_thread_num();
	for(int tid_e = tid*num_thr_e/num_thr; tid_e < (tid+1)*num_thr_e/num_thr; tid_e++)
	{
		uint64_t start = uint64_t(tid_e)*(numreads/num_thr_e);
		uint8_t* read_lengths;
		read_lengths = readlengths + start;
		struct qv_options_t opts;
		opts.verbose = 0;
		std::string input_file_string = (infile_quality+"."+std::to_string(tid_e));
		char *input_file = new char [input_file_string.length()+1];
		strcpy(input_file,input_file_string.c_str());
		std::string output_file_string = (outfile_quality+"."+std::to_string(tid_e));
		char *output_file = new char [output_file_string.length()+1];
		strcpy(output_file,output_file_string.c_str());
		decode(input_file,output_file,&opts,read_lengths);
	}
	}
	return;
}

void load_readlengths(uint8_t *readlengths)
{
	if(preserve_order == "True")
	{
		std::ifstream in_readlength(infile_readlength,std::ios::binary);
		std::ifstream in_order(infile_order,std::ios::binary);
		uint8_t cur_readlen;
		uint32_t order;
		for(uint32_t i = 0; i < numreads; i++)
		{
			in_readlength.read((char*)&cur_readlen,sizeof(uint8_t));
			in_order.read((char*)&order,sizeof(uint32_t));
			readlengths[order] = cur_readlen;
		}
		in_readlength.close();
		in_order.close();
	}
	else
	{
		std::ifstream in_readlength(infile_readlength,std::ios::binary);
		uint32_t order;
		for(uint32_t i = 0; i < numreads; i++)
		{
			in_readlength.read((char*)&readlengths[i],sizeof(uint8_t));
		}
		in_readlength.close();
	}
	return;
}
