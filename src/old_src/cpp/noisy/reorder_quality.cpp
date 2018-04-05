#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>
#include <string>
#include "config.h"

uint32_t numreads;

std::string outfile;
std::string infile;
std::string infile_N;
std::string infile_order;

void reorder_quality();

void getDataParams();//populate numreads and readlen

int main(int argc, char** argv)
{
	std::string basedir = std::string(argv[1]);
	outfile = basedir + "/output/output.quality";
	infile = basedir + "/output/input_clean.quality";
	infile_N = basedir + "/output/input_N.quality";
	infile_order = basedir + "/output/read_order.bin";
	getDataParams();
	reorder_quality();
	return 0;
}

void reorder_quality()
{
	std::ofstream f(outfile);
	std::ifstream f_in(infile);
	std::ifstream f_N(infile_N);
	std::ifstream f_order(infile_order,std::ios::binary);
	uint32_t order;
	uint32_t *reverse_index = new uint32_t [numreads];
	for (uint32_t i = 0; i < numreads; i++)
	{
		f_order.read((char*)&order,sizeof(uint32_t));
		reverse_index[order] = i;
	}
	f_order.close();
	uint32_t max_bin_size = numreads/4;
	char s[readlen+1];
	s[readlen] = '\0';
	for (uint32_t i = 0; i <= numreads/max_bin_size; i++)
	{
		auto numreads_bin = max_bin_size;
		if (i == numreads/max_bin_size)
			numreads_bin = numreads%max_bin_size;
		uint32_t *index_array = new uint32_t [numreads_bin];
		char(*quality_bin)[readlen+1] = new char [numreads_bin][readlen+1];	
		uint32_t pos = 0;
		for(uint32_t j = 0; j < numreads; j++)
		{
			order = reverse_index[j];
			if (order >= i*max_bin_size && order < i*max_bin_size + numreads_bin)
			{
				index_array[order-i*max_bin_size] = pos;
				f_in.seekg(uint64_t(j)*(readlen+1), f_in.beg);
				f_in.getline(quality_bin[pos],readlen+1);
				pos++;
			}
		}
		for(uint32_t j = 0; j < numreads_bin; j++)
		{
			f << quality_bin[index_array[j]] << "\n";
		}
		delete[] index_array;
		delete[] quality_bin;
	}
	f_in.close();
	f << f_N.rdbuf();
	delete[] reverse_index;
	f_N.close();
	f.close();
	return;
}

void getDataParams()
{
	uint32_t number_of_lines = 0;
	std::string line;
	std::ifstream myfile(infile, std::ifstream::in);
	while (std::getline(myfile, line))
		++number_of_lines;
	numreads = number_of_lines;
	myfile.close();
}
