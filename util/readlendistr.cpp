#include <iostream>
#include <fstream>
#include <string>

int main()
{
	std::ifstream f("NA12877-CP09_S9_L001_R2_001.fastq");
	long linenum = 0;
	long a[256];
	for(int i = 0; i < 256; i++)
		a[i] = 0;
	std::string line;
	while(std::getline(f,line))
	{
		linenum++;
		if(linenum%4 == 2)
			a[line.length()]++;
	}
	for(int i = 0; i < 256; i++)
		std::cout << i <<"\t" << a[i] <<"\n";
	return 0;
}
