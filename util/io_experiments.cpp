#include<iostream>
#include<fstream>
#include<string>

int fun_1(std::string &infile);
int fun_2(std::string &infile, std::string &infile_2);
int fun_3(std::string &infile, std::string &infile_2,int block_size);
int fun_4();

int main(int argc, char** argv)
{
	std::string infile;
	if(argc > 1)
		infile = std::string(argv[1]);
//	fun_1(infile);
	std::string infile_2;
	if(argc > 2)
		infile_2 = std::string(argv[2]);
//	fun_2(infile, infile_2);
//	fun_3(infile,infile_2,1000000);
	fun_4();
}

int fun_1(std::string &infile) {
	std::ifstream f_in(infile);
	std::string a;
	long numlines = 0, num_char = 0;
	while(std::getline(f_in,a)) {
		num_char += a.size();
		numlines+=1;
	}
	std::cout << num_char << "\n";
	std::cout << numlines << "\n";
}

int fun_2(std::string &infile, std::string &infile_2) {
	std::ifstream f_in(infile);
	std::ifstream f_in_2(infile_2);
	std::string a;
	long numlines = 0, num_char = 0;
	long numlines_2 = 0, num_char_2 = 0;
	while(std::getline(f_in,a)) {
		num_char += a.size();
		numlines+=1;
		std::getline(f_in_2,a);
		num_char_2 += a.size();
		numlines_2+=1;
	}
	std::cout << num_char << "\n";
	std::cout << numlines << "\n";
	std::cout << num_char_2 << "\n";
	std::cout << numlines_2 << "\n";
}

int fun_3(std::string &infile, std::string &infile_2, int block_size) {
	std::ifstream f_in(infile);
	std::ifstream f_in_2(infile_2);
	std::string a;
	long numlines = 0, num_char = 0;
	long numlines_2 = 0, num_char_2 = 0;
	bool done_1 = false, done_2 = false;
	while(!done_1 || !done_2) {
		//read block_size lines from each file
		if(!done_1)
		for(int i = 0; i < block_size; i++) {
			if(std::getline(f_in,a)) {
				num_char += a.size();
				numlines+=1;
			}
			else
				done_1 = true;
		}
		if(!done_2)
		for(int i = 0; i < block_size; i++) {
			if(std::getline(f_in_2,a)) {
				num_char_2 += a.size();
				numlines_2+=1;
			}
			else
				done_2 = true;
		}
	}
	std::cout << num_char << "\n";
	std::cout << numlines << "\n";
	std::cout << num_char_2 << "\n";
	std::cout << numlines_2 << "\n";
}

int fun_4() {
	std::ofstream f_out("temptemp");
	for (long i = 0; i < 940000000; i++)
		f_out << std::string (100,'r') << "\n";
}
	
