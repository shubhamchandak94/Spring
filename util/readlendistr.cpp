#include <iostream>
#include <fstream>
#include <string>

int main(int argc, char **argv)
{
	std::string file_name = std::string(argv[1]);
	std::ifstream f(file_name);
	long linenum = 0;
	long a[512];
	for(int i = 0; i < 512; i++)
		a[i] = 0;
	std::string line;
	while(std::getline(f,line))
	{
		linenum++;
		if(linenum%4 == 2)
			a[line.length()]++;
	}
	for(int i = 0; i < 512; i++)
		std::cout << i <<"\t" << a[i] <<"\n";
	return 0;
}
