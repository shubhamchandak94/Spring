#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

int main(int argc, char **argv)
{
	std::string infile = std::string(argv[1]);
	int max_readlen = atoi(argv[2]);
	std::string s[4];
	std::string outfile = infile + ".split";	
	std::ifstream f_in(infile);
	std::ofstream f_out(outfile);
	while(std::getline(f_in,s[0]))
	{
		std::getline(f_in,s[1]);
		std::getline(f_in,s[2]);
		std::getline(f_in,s[3]);
		if(s[1].length() != s[3].length())
		{
			std::cout << "Quality length does not match read length\n";
			return -1;
		}
		int i;
		for(i = 0; i < (s[1].length()+max_readlen-1)/max_readlen - 1; i++) 
		//rounding up s[1].length()/max_readlen and subtracting 1
		{
			f_out << s[0] << "\n" << s[1].substr(i*max_readlen,max_readlen) << "\n+\n" << s[3].substr(i*max_readlen,max_readlen) <<"\n";
		}
		//last part (only part if readlen <= max_readlen)
		f_out << s[0] << "\n" << s[1].substr(i*max_readlen,s[1].length()-i*max_readlen) << "\n+\n" << s[3].substr(i*max_readlen,s[3].length()-i*max_readlen) <<"\n";
	}
	f_in.close();
	f_out.close();
	return 0;	
}	
