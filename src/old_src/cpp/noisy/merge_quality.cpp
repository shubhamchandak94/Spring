#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>
#include <string>


std::string outfile;
std::string infile;
std::string infile_N;
std::string infile_order_N;

void merge_quality();


int main(int argc, char** argv)
{
	std::string basedir = std::string(argv[1]);
	outfile = basedir + "/output/output.quality";
	infile = basedir + "/output/input_clean.quality";
	infile_N = basedir + "/output/input_N.quality";
	infile_order_N = basedir + "/output/read_order_N.bin";
	merge_quality();
	return 0;
}

void merge_quality()
{
	std::ofstream f(outfile);
	std::ifstream f_in(infile);
	std::ifstream f_N(infile_N);
	std::ifstream f_order_N(infile_order_N,std::ios::binary);

	std::string line_N,line;
	uint32_t next,i=0;
	f_order_N.read((char*)&next,sizeof(uint32_t));
	while (std::getline(f_in, line))
	{
		while(i==next && !f_order_N.eof())
		{
			std::getline(f_N, line_N);
			f << line_N << "\n";
			i++;
			f_order_N.read((char*)&next,sizeof(uint32_t));
		}
		f << line << "\n";
		i++;
	}	
	while(i==next && !f_order_N.eof())
	{
		std::getline(f_N, line_N);
		f << line_N << "\n";
		i++;
		f_order_N.read((char*)&next,sizeof(uint32_t));
	}
	f.close();
	f_in.close();
	f_N.close();
	f_order_N.close();
	return;
}

