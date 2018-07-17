#include <fstream>
#include <iostream>
#include <string>
#include <cstdio>
#include <cmath>
#include <algorithm>

std::string infile;

void unpack_order();//unpack order into 32 bits per uint

int main(int argc, char** argv)
{
	std::string basedir = std::string(argv[1]);
	infile = basedir + "/read_order.bin";
	unpack_order();
	return 0;
}

void unpack_order()
{
	std::ofstream f_out(infile+".tmp",std::ios::binary);
	std::ifstream f_in(infile,std::ios::binary);
	std::ifstream f_tail(infile+".tail",std::ios::binary);
	
	int numbits;
	uint32_t numreads;
	f_in.read((char*)&numbits,sizeof(int));
	f_in.read((char*)&numreads,sizeof(uint32_t));
	uint32_t order;
	uint32_t order_array[32];
	uint32_t *input_array = new uint32_t [numbits];
	for(uint32_t i = 0; i < numreads/32; i++)
	{
		for(int k = 0; k < numbits; k++)
			f_in.read((char*)&input_array[k], sizeof(uint32_t));
		
		int pos_in_int = 0, pos_in_array = 0;
		
		uint32_t mask = (0xFFFFFFFF)>>(32-numbits);		
		for(int k = 0; k < 32; k++)
		{
			if(32-pos_in_int > numbits)
			{
				order_array[k] = (input_array[pos_in_array]>>pos_in_int) & mask;
				pos_in_int += numbits;
			}
			else if(32-pos_in_int == numbits)
			{
				order_array[k] = (input_array[pos_in_array]>>pos_in_int) & mask;
				pos_in_array++;
				pos_in_int = 0;
			}
			else
			{
				order_array[k] = (input_array[pos_in_array]>>pos_in_int) & mask;
				pos_in_array++;
				order_array[k] += (input_array[pos_in_array]<<(32-pos_in_int)) & mask;
				pos_in_int = numbits - (32-pos_in_int);
			}
		}
		for(int k = 0; k < 32; k++)
			f_out.write((char*)&order_array[k], sizeof(uint32_t));
	}
	f_out << f_tail.rdbuf();
	f_tail.close();
	f_in.close();
	f_out.close();
	remove(infile.c_str());
	rename((infile+".tmp").c_str(),infile.c_str());
	return;
}
