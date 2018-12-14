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
std::string infile_id[2];
std::string outfile_id[2];
std::string infile_quality[2];
std::string outfile_quality[2];
std::string infile_readlength;
std::string infile_order_1;
std::string infile_order_2;
std::string infile_paired_flag_first;

std::string infilenumreads;

int max_readlen, num_thr, num_thr_e;
uint32_t numreads, numreads_by_2;
uint8_t paired_id_code;
std::string preserve_order, quality_mode, preserve_quality, preserve_id;

void decompress_id();
void decompress_quality(uint8_t *readlengths_1, uint8_t *readlengths_2);
void load_readlengths(uint8_t *readlengths_1, uint8_t *readlengths_2);

void decode(char *input_file, char *output_file, struct qv_options_t *opts, uint8_t *read_lengths);

int main(int argc, char** argv)
{
	basedir = std::string(argv[1]);
	infile_id[0] = basedir + "/compressed_id_1.bin";
	infile_id[1] = basedir + "/compressed_id_2.bin";
	outfile_id[0] = basedir + "/id_1.txt";
	outfile_id[1] = basedir + "/id_2.txt";
	infile_quality[0] = basedir + "/compressed_quality_1.bin";
	infile_quality[1] = basedir + "/compressed_quality_2.bin";
	outfile_quality[0] = basedir + "/quality_1.txt";
	outfile_quality[1] = basedir + "/quality_2.txt";
	infile_readlength = basedir + "/read_lengths.bin";
	infile_order_1 = basedir + "/read_order.bin";
	infile_order_2 = basedir + "/read_order_2.bin";
	infile_paired_flag_first = basedir + "/read_paired_flag_first.bin";

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
	f_numreads.read((char*)&paired_id_code,sizeof(uint8_t));
	f_numreads.close();
	numreads_by_2 = numreads/2;

	std::cout << "Decompressing Quality and/or IDs\n";
	omp_set_num_threads(num_thr);
	if(preserve_quality == "True" && !(quality_mode == "bsc" || quality_mode == "illumina_binning_bsc"))
	{
		uint8_t *readlengths_1 = new uint8_t[numreads_by_2];
		uint8_t *readlengths_2 = new uint8_t[numreads_by_2];
		load_readlengths(readlengths_1, readlengths_2);

		decompress_quality(readlengths_1, readlengths_2);
		delete[] readlengths_1;
		delete[] readlengths_2;
	}
	if(preserve_id == "True")
	{
		decompress_id();
	}
}

void decompress_id()
{
	for(int k = 0; k < 2; k++)
	{
		if(paired_id_code != 0 && k==1)
			break;
		#pragma omp parallel
		{
		int tid = omp_get_thread_num();
		for(int tid_e = tid*num_thr_e/num_thr; tid_e < (tid+1)*num_thr_e/num_thr; tid_e++)
		{
			uint32_t numreads_thr = numreads_by_2/num_thr_e;
			if(tid_e == num_thr_e - 1)
				numreads_thr = numreads_by_2-numreads_thr*(num_thr_e-1); 
			struct compressor_info_t comp_info;
			comp_info.numreads = numreads_thr;
			comp_info.mode = DECOMPRESSION;
			comp_info.fcomp = fopen((infile_id[k]+"."+std::to_string(tid_e)).c_str(),"r");
			comp_info.f_id = fopen((outfile_id[k]+"."+std::to_string(tid_e)).c_str(),"w");
			decompress((void *)&comp_info);
			fclose(comp_info.fcomp);
			fclose(comp_info.f_id);
		}
		}	

	}
	return;
}

void decompress_quality(uint8_t *readlengths_1, uint8_t *readlengths_2)
{
	for(int k = 0; k < 2; k++)
	{	
		#pragma omp parallel
		{
		int tid = omp_get_thread_num();
		for(int tid_e = tid*num_thr_e/num_thr; tid_e < (tid+1)*num_thr_e/num_thr; tid_e++)
		{
			uint64_t start = uint64_t(tid_e)*(numreads_by_2/num_thr_e);
			uint8_t* read_lengths;
			if(k == 0)
				read_lengths = readlengths_1 + start;
			else
				read_lengths = readlengths_2 + start;
			struct qv_options_t opts;
			opts.verbose = 0;
			std::string input_file_string = (infile_quality[k]+"."+std::to_string(tid_e));
			char *input_file = new char [input_file_string.length()+1];
			strcpy(input_file,input_file_string.c_str());
			std::string output_file_string = (outfile_quality[k]+"."+std::to_string(tid_e));
			char *output_file = new char [output_file_string.length()+1];
			strcpy(output_file,output_file_string.c_str());
			decode(input_file,output_file,&opts,read_lengths);
		}
		}
	}
	return;
}

void load_readlengths(uint8_t *readlengths_1, uint8_t *readlengths_2)
{
	if(preserve_order == "True")
	{
		std::ifstream in_readlength(infile_readlength,std::ios::binary);
		std::ifstream in_flag_first(infile_paired_flag_first);
		std::ifstream in_order_1(infile_order_1,std::ios::binary);
		std::ifstream in_order_2(infile_order_2,std::ios::binary);
		uint8_t cur_readlen;
		char flag_first;
		uint32_t order;
		for(uint32_t i = 0; i < numreads; i++)
		{
			in_readlength.read((char*)&cur_readlen,sizeof(uint8_t));
			in_flag_first >> flag_first;
			if(flag_first == '1')
			{
				in_order_1.read((char*)&order,sizeof(uint32_t));
				readlengths_1[order] = cur_readlen;
			}
			else
			{
				in_order_2.read((char*)&order,sizeof(uint32_t));
				readlengths_2[order] = cur_readlen;
			}
		}
		in_readlength.close();
		in_flag_first.close();
		in_order_1.close();
		in_order_2.close();
	}
	else
	{
		std::ifstream in_readlength(infile_readlength,std::ios::binary);
		std::ifstream in_flag_first(infile_paired_flag_first);
		std::ifstream in_order_2(infile_order_2,std::ios::binary);
		uint8_t cur_readlen;
		char flag_first;
		uint32_t j = 0;
		uint32_t order;
		for(uint32_t i = 0; i < numreads; i++)
		{
			in_readlength.read((char*)&cur_readlen,sizeof(uint8_t));
			in_flag_first >> flag_first;
			if(flag_first == '1')
			{
				readlengths_1[j] = cur_readlen;
				j++;
			}
			else
			{
				in_order_2.read((char*)&order,sizeof(uint32_t));
				readlengths_2[order] = cur_readlen;
			}
		}
		in_readlength.close();
		in_flag_first.close();
		in_order_2.close();
	}
	return;
}
