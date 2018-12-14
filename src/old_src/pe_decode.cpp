#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <string>
#include <cstdio>

std::string outfile_order_2;
std::string infilenumreads;
std::string infile_order_paired;
std::string infile_paired_flag_first;
std::string infile_order;

uint32_t numreads, numreads_by_2;
std::string preserve_order;
void generate_order_from_paired();
void unpackbits();

int main(int argc, char** argv)
{
	std::string basedir = std::string(argv[1]);
	outfile_order_2 = basedir + "/read_order_2.bin";
	infilenumreads = basedir + "/numreads.bin";
	infile_order_paired = basedir + "/read_order_paired.bin";
	infile_paired_flag_first = basedir + "/read_paired_flag_first.bin";
	infile_order = basedir + "/read_order.bin";	
	
	preserve_order = std::string(argv[2]);
	
	std::ifstream f_numreads(infilenumreads, std::ios::binary);
	f_numreads.seekg(4);
	f_numreads.read((char*)&numreads,sizeof(uint32_t));
	f_numreads.close();
	numreads_by_2 = numreads/2;
	unpackbits();
	generate_order_from_paired();
}

void generate_order_from_paired()
{
	bool *flag_first = new bool [numreads];
	uint32_t *read_order = new uint32_t[numreads];
	std::fill(read_order, read_order+numreads, numreads);
	std::ifstream in_order_paired(infile_order_paired,std::ios::binary);
	std::ifstream in_paired_flag_first(infile_paired_flag_first);
	
	char c;
	uint32_t order_paired;
	uint32_t current_pair = 0;
	//decode order_paired and flag_first into read_order 
	for(uint32_t i = 0; i < numreads; i++)
	{
		if(read_order[i] != numreads)//this position already filled
			continue;
		in_paired_flag_first.get(c);
		in_order_paired.read((char*)&order_paired, sizeof(uint32_t));
		if(c == '1')
		{
			read_order[i] = current_pair;
			flag_first[i] = 1;
			read_order[i+order_paired] = current_pair+numreads_by_2;
			flag_first[i+order_paired] = 0;
		}
		else
		{
			if(c!='0') std::cout << c << "\n";
			read_order[i] = current_pair+numreads_by_2;
			flag_first[i] = 0;
			read_order[i+order_paired] = current_pair;
			flag_first[i+order_paired] = 1;
		}
		current_pair++;	
	}
	in_paired_flag_first.close();
	in_order_paired.close();

	//write flag_first to file and replace old file
	std::ofstream out_paired_flag_first(infile_paired_flag_first+".tmp");
	for(uint32_t i = 0; i < numreads; i++)
		out_paired_flag_first << flag_first[i];
	out_paired_flag_first.close();	
	remove(infile_paired_flag_first.c_str());
	rename((infile_paired_flag_first+".tmp").c_str(),infile_paired_flag_first.c_str());

	
	//now we'll need the inverse index, before doing that write to file to save memory
	std::ofstream f_temp(outfile_order_2,std::ios::binary);
	for(uint32_t i = 0; i < numreads; i++)
	{
		f_temp.write((char*)&read_order[i], sizeof(uint32_t));
	}
	f_temp.close();
	
	//build inverse index
	std::ifstream fin_temp(outfile_order_2,std::ios::binary);
	for(uint32_t i = 0; i < numreads; i++)
	{
		fin_temp.read((char*)&order_paired, sizeof(uint32_t));
		read_order[order_paired] = i;
	}
	fin_temp.close();
	
	//now, go through the reads in order, and for each read which is 1st in paired end, 
	//store the location of its pair 
	f_temp.open(outfile_order_2+".tmp",std::ios::binary);
	fin_temp.open(outfile_order_2,std::ios::binary);
	for(uint32_t i = 0; i < numreads; i++)
	{
		fin_temp.read((char*)&order_paired, sizeof(uint32_t));
		if(flag_first[i])
		{
			uint32_t pos_of_pair = read_order[order_paired+numreads_by_2];
			f_temp.write((char*)&pos_of_pair, sizeof(uint32_t));
		}
	}
	fin_temp.close();
	f_temp.close();
	
	//store the inverse index for the 2nd reads of pair
	fin_temp.open(outfile_order_2+".tmp",std::ios::binary);
	for(uint32_t i = 0; i < numreads_by_2; i++)
	{
		fin_temp.read((char*)&order_paired, sizeof(uint32_t));
		read_order[order_paired] = i;
	}
	fin_temp.close();
	remove((outfile_order_2+".tmp").c_str());
	
	//Finally, write the order for second reads of pair 
	f_temp.open(outfile_order_2,std::ios::binary);
	uint32_t j = 0;
	for(uint32_t i = 0; i < numreads; i++)
	{
		if(!flag_first[i])
		{
			f_temp.write((char*)&read_order[i], sizeof(uint32_t));
			read_order[j] = read_order[i];//storing for preserve_order step
			j++;
		}
	}
	f_temp.close();

	if(preserve_order == "True")//in this case, we need to map the order obtained above to the true positions
	{
		//first load from read_order.bin into second half of read_order array
		fin_temp.open(infile_order,std::ios::binary);
		for(uint32_t i = 0; i < numreads_by_2; i++)
			fin_temp.read((char*)&read_order[numreads_by_2+i], sizeof(uint32_t));
		fin_temp.close();
	
		//now compose the two permutations and write result to read_order_2.bin
		f_temp.open(outfile_order_2,std::ios::binary);
		for(uint32_t i = 0; i < numreads_by_2; i++)	
			f_temp.write((char*)&read_order[numreads_by_2+read_order[i]], sizeof(uint32_t));
		f_temp.close();
	}
	delete[] read_order;
	return;
}

void unpackbits()
{
	//flag_first
	std::ifstream in_flag_first(infile_paired_flag_first,std::ios::binary);
	std::ofstream f_flag_first(infile_paired_flag_first+".tmp");
	std::ifstream in_flag_first_tail(infile_paired_flag_first+".tail");
	char inttochar[2];
	inttochar[0] = '0';
	inttochar[1] = '1';
	
	uint8_t packedchar;
	in_flag_first.read((char*)&packedchar,sizeof(uint8_t));
	while(!in_flag_first.eof())
	{	
		f_flag_first << inttochar[packedchar%2];
		packedchar/=2; 	
		f_flag_first << inttochar[packedchar%2];
		packedchar/=2; 	
		f_flag_first << inttochar[packedchar%2];
		packedchar/=2; 	
		f_flag_first << inttochar[packedchar%2];
		packedchar/=2; 	
		f_flag_first << inttochar[packedchar%2];
		packedchar/=2; 	
		f_flag_first << inttochar[packedchar%2];
		packedchar/=2; 	
		f_flag_first << inttochar[packedchar%2];
		packedchar/=2; 	
		f_flag_first << inttochar[packedchar%2];
		packedchar/=2; 	
		in_flag_first.read((char*)&packedchar,sizeof(uint8_t));
	}
	in_flag_first.close();
	f_flag_first << in_flag_first_tail.rdbuf();
	in_flag_first_tail.close();
	f_flag_first.close();
	remove(infile_paired_flag_first.c_str());
	remove((infile_paired_flag_first+".tail").c_str());
	rename((infile_paired_flag_first+".tmp").c_str(),infile_paired_flag_first.c_str());
	return;
}

