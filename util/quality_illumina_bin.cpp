#include<iostream>
#include<fstream>
#include<string>

char illumina_binning_table[128];
void generate_illumina_binning_table();
void illumina_binning(std::string &quality);

int main(int argc, char** argv)
{
	generate_illumina_binning_table();
	std::string infile = std::string(argv[1]);
	std::string outfile = std::string(argv[2]);
	std::ifstream f_in(infile);
	std::ofstream f_out(outfile);
	std::string quality;
	while(std::getline(f_in,quality))
	{
		illumina_binning(quality);
		f_out << quality << "\n";
	}
	f_in.close();
	f_out.close();
	return 0;
}

void illumina_binning(std::string &quality)
{
	for(uint8_t i = 0; i < quality.length(); i++)
		quality[i] = illumina_binning_table[quality[i]];
	return;
}

void generate_illumina_binning_table()
{
	for(uint8_t i = 0; i <= 33+1; i++)
		illumina_binning_table[i] = 33+0;
	for(uint8_t i = 33+2; i <= 33+9; i++)
		illumina_binning_table[i] = 33+6;
	for(uint8_t i = 33+10; i <= 33+19; i++)
		illumina_binning_table[i] = 33+15;
	for(uint8_t i = 33+20; i <= 33+24; i++)
		illumina_binning_table[i] = 33+22;
	for(uint8_t i = 33+25; i <= 33+29; i++)
		illumina_binning_table[i] = 33+27;
	for(uint8_t i = 33+30; i <= 33+34; i++)
		illumina_binning_table[i] = 33+33;
	for(uint8_t i = 33+35; i <= 33+39; i++)
		illumina_binning_table[i] = 33+37;
	for(uint8_t i = 33+40; i <= 127; i++)
		illumina_binning_table[i] = 33+40;
	return;
}
