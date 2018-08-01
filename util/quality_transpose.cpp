#include<iostream>
#include<fstream>
#include<string>
#include<cstdlib>

int main(int argc, char** argv)
{
	std::string infile = std::string(argv[1]);
	std::string outfile = infile+".transposed";
	int readlen = atoi(argv[2]);
	long numreads = atol(argv[3]);
	char *quality_array = new char[numreads*(readlen+1)];
	std::ifstream f_in(infile);
	//read file
	for(long i = 0 ; i < numreads; i++)
		f_in.getline((quality_array+i*(readlen+1)),readlen+1);
	f_in.close();
	//write file
	std::ofstream f_out(outfile);
	for(int j = 0; j < readlen; j++)
	{
		for(long i = 0; i < numreads; i++)
			f_out << quality_array[i*(readlen+1)+j];
		f_out << "\n";
	}
	f_out.close();
	return 0;
}
