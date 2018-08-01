#include<iostream>
#include<fstream>
#include<string>

int main(int argc, char** argv)
{
	std::string infile = std::string(argv[1]);
	std::ifstream f_in(infile+".packed",std::ios::binary);
	std::ifstream f_in_tail(infile+".packed.tail");
	std::ofstream f_out(infile+".unpacked");
	char c;
	char int_to_qual[4];
	int_to_qual[0] = 'F';
	int_to_qual[1] = ':';
	int_to_qual[2] = ',';
	int_to_qual[3] = '#';
	uint8_t quality_bin;
	f_in.read((char*)&quality_bin,sizeof(uint8_t));
	while(!f_in.eof())
	{
		f_out << int_to_qual[quality_bin%4];
		quality_bin/=4;
		f_out << int_to_qual[quality_bin%4];
		quality_bin/=4;
		f_out << int_to_qual[quality_bin%4];
		quality_bin/=4;
		f_out << int_to_qual[quality_bin%4];
		quality_bin/=4;
		f_in.read((char*)&quality_bin,sizeof(uint8_t));
	}
	f_in.close();
	f_out << f_in_tail.rdbuf();
	f_out.close();
	f_in_tail.close();
	return 0;
}
