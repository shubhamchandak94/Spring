#include <iostream>
#include <fstream>
#include <string>

std::string infile;
std::string outfileclean;
std::string outfileN;
std::string outfileorderN;
std::string outfilequality;
std::string outfilequalityN;

void preprocess();

int main(int argc, char** argv)
{
	std::string basedir = std::string(argv[2]);
	
	infile = std::string(argv[1]);
	outfileclean = basedir + "/input_clean.dna";
	outfileN = basedir + "/input_N.dna";
	outfileorderN = basedir + "/read_order_N.bin";
	outfilequality = basedir + "/input_clean.quality";
	outfilequalityN = basedir + "/input_N.quality";
	preprocess();
	std::cout << "Preprocessing Done!\n";
	return 0;
}

void preprocess()
{
	std::string line;
	std::ifstream myfile(infile, std::ifstream::in);
	std::ofstream f_clean(outfileclean);
	std::ofstream f_N(outfileN);
	std::ofstream f_order_N(outfileorderN,std::ios::binary);
	std::ofstream f_quality(outfilequality);
	std::ofstream f_quality_N(outfilequalityN);
	int i = 0;
	uint32_t readnum = 0;
	bool flag_N = false;
	while(std::getline(myfile, line))
	{
		switch(i)
		{
			case 0:	//f_id << line << "\n";
				break;
			case 1: //f << line << "\n";
				if(line.find('N')!=std::string::npos)
				{
					flag_N = true;
					f_N << line << "\n";
					f_order_N.write((char*)&readnum,sizeof(uint32_t));
				}
				else
				{
					flag_N = false;
					f_clean << line << "\n";
				}
				break;
			case 2: break;
			case 3:	//f_quality << line << "\n";
				if(!flag_N)
					f_quality << line << "\n";
				else
					f_quality_N << line << "\n";
				readnum++;
				break;
		}
		i = (i+1)%4;
	}
	return;	
}
