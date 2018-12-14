#include<iostream>
#include<fstream>
#include<string>

int main(int argc, char** argv)
{
	std::string infile = std::string(argv[1]);
	std::ifstream f_in;
	std::ofstream f_out(infile+".packed",std::ios::binary);
	std::ofstream f_out_tail(infile+".packed.tail");
	char c;
	uint8_t qual_to_int[128];
	qual_to_int['F'] = 0;
	qual_to_int[':'] = 1;
	qual_to_int[','] = 2;
	qual_to_int['#'] = 3;
	
	long file_len = 0;
	f_in.open(infile,std::ios::ate);
	file_len = f_in.tellg();
	std::cout << file_len << "\n";
	f_in.close();
	f_in.open(infile);
	char quality[8];
	uint8_t quality_bin;
	for(long i = 0; i < file_len/4; i++)
	{
		f_in.read(quality,4);
		quality_bin = 64*qual_to_int[quality[3]]+16*qual_to_int[quality[2]]+4*qual_to_int[quality[1]]+qual_to_int[quality[0]];
		f_out.write((char*)&quality_bin, sizeof(uint8_t));
	}
	f_out.close();
	f_in.read(quality,file_len%4);
	for(int i=0; i<file_len%4;i++)
		f_out_tail << quality[i];
	f_out_tail.close();
	f_in.close();
	return 0;
}
