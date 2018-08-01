#include<iostream>
#include<fstream>
#include<string>

int main(int argc, char** argv)
{
	std::string infile = std::string(argv[1]);
	std::ifstream f_in(infile);
	std::ofstream f_out_len(infile+".rl.len",std::ios::binary);
	std::ofstream f_out_char(infile+".rl.char");
	char c;
	uint16_t run_length = 0;
	while(f_in >> c)
	{
		if(c == 'F')
			run_length++;
		else
		{
			f_out_len.write((char*)&run_length, sizeof(uint16_t));
			f_out_char << c;
			run_length = 0;
		}
	}
	if(run_length != 0)
		f_out_len.write((char*)&run_length, sizeof(uint16_t));
	return 0;
}
