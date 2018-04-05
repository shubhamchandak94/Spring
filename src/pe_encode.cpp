#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdio>

std::string infile_order;
std::string infilenumreads;

std::string outfile_order_paired;
std::string outfile_paired_flag_first;

uint32_t numreads, numreads_by_2; 

void populate_arrays(uint32_t* read_order, uint32_t* read_inverse_order);
//populate arrays:
//read_order = pos in reordered file to pos in original file
//read_inverse_order = pos in original file to pos in reordered file

void write_order_paired(uint32_t* read_order, uint32_t* read_inverse_order);
//write to all output files
//order_paired - store relative position of paired read once per pair (store for the read occuring first in the reordered file)
//For each pair, paired_flag_first stores 1 is 1st read comes first.

void packbits();
//pack flag files into 1 bit per flag

void generate_order_preserve(uint32_t* read_order);
//generate order file for half the reads

int main(int argc, char** argv)
{
	std::string basedir = std::string(argv[1]);
	std::string preserve_order = std::string(argv[2]);
	infilenumreads = basedir + "/numreads.bin";
	infile_order = basedir + "/read_order.bin";
	outfile_order_paired = basedir + "/read_order_paired.bin";
	outfile_paired_flag_first = basedir + "/read_paired_flag_first.bin";
	
	std::ifstream f_numreads(infilenumreads, std::ios::binary);
	f_numreads.seekg(4);
	f_numreads.read((char*)&numreads,sizeof(uint32_t));
	f_numreads.close();
	numreads_by_2 = numreads/2;	
	
	uint32_t *read_order = new uint32_t [numreads];
	uint32_t *read_inverse_order = new uint32_t[numreads];
	populate_arrays(read_order, read_inverse_order);
    	write_order_paired(read_order, read_inverse_order);
	packbits();

	if(preserve_order == "True")
		generate_order_preserve(read_order);
		
	delete[] read_order;
	delete[] read_inverse_order;
	return 0;
}

void populate_arrays(uint32_t* read_order, uint32_t* read_inverse_order)
{
	//read file read_order
	std::ifstream f_order(infile_order, std::ios::binary);
	for(uint32_t i = 0; i < numreads; i++)
	{
		f_order.read((char*)&read_order[i], sizeof(uint32_t));
	}
	f_order.close();
	
	//now fill read_inverse_order
	for(uint32_t i = 0; i < numreads; i++)
	{
		read_inverse_order[read_order[i]] = i;
	}
	return;
}	

void write_order_paired(uint32_t* read_order, uint32_t* read_inverse_order)
{
	std::ofstream f_flag_first(outfile_paired_flag_first);
	std::ofstream f_order_paired(outfile_order_paired,std::ios::binary);
	for(uint32_t i = 0; i < numreads; i++)
	{
		if(read_order[i] < numreads_by_2)//first read of pair
		{
			if(read_inverse_order[read_order[i] + numreads_by_2] > i)//pair not already seen
			{
				uint32_t temp = (read_inverse_order[read_order[i] + numreads_by_2]-i);
				f_order_paired.write((char*)&temp,sizeof(uint32_t));
				f_flag_first << '1';
			}
		}
		else
		{
			if(read_inverse_order[read_order[i] - numreads_by_2] > i)//pair not already seen
			{
				f_flag_first << '0';	
				uint32_t temp = (read_inverse_order[read_order[i] - numreads_by_2]-i);	
				f_order_paired.write((char*)&temp,sizeof(uint32_t));
			}
		}			
	}
	f_flag_first.close();
	f_order_paired.close();
}

void packbits()
{
	//flag_first
	std::ifstream in_flag_first(outfile_paired_flag_first);
	std::ofstream f_flag_first(outfile_paired_flag_first+".tmp",std::ios::binary);
	std::ofstream f_flag_first_tail(outfile_paired_flag_first+".tail");
	
	uint8_t chartoint[128];
	chartoint['0'] = 0;
	chartoint['1'] = 1;
	in_flag_first.close();
	in_flag_first.open(outfile_paired_flag_first);
	char chararray[8];
	uint8_t packedchar;
	for(uint64_t i = 0; i < numreads_by_2/8; i++)
	{
		in_flag_first.read(chararray,8);
		
		packedchar = 128*chartoint[chararray[7]]+ 64*chartoint[chararray[6]]
				+ 32*chartoint[chararray[5]]+ 16*chartoint[chararray[4]]
				+  8*chartoint[chararray[3]]+  4*chartoint[chararray[2]]
				+  2*chartoint[chararray[1]]+  1*chartoint[chararray[0]];
		f_flag_first.write((char*)&packedchar,sizeof(uint8_t));
	}
	f_flag_first.close();
	in_flag_first.read(chararray,numreads_by_2%8);
	for(int i=0; i<numreads_by_2%8;i++)
		f_flag_first_tail << chararray[i];
	f_flag_first_tail.close();
	in_flag_first.close();
	remove((outfile_paired_flag_first).c_str());
	rename((outfile_paired_flag_first+".tmp").c_str(),(outfile_paired_flag_first).c_str());		
	return;
}

void generate_order_preserve(uint32_t* read_order)
{
	std::ofstream fout_order(infile_order,std::ios::binary);
	uint32_t order;
	for(uint32_t i = 0; i < numreads; i++)
	{
		if(read_order[i] < numreads_by_2)
		{
			fout_order.write((char*)&read_order[i],sizeof(uint32_t));
		}
	}
	fout_order.close();
	return;
}
