#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>
#include <string>


long readlen;

std::string outfile;
std::string outfile_clean;
std::string infile_seq;
std::string infile_meta;
std::string infile_pos;
std::string infile_noise;
std::string infile_noisepos;
std::string infile_rev;
std::string infile_order;

char longtochar[] = {'A','C','G','T'};
long chartolong[128];
char dec_noise[128][128];
char chartorevchar[128];

void decode();

void restore_order();

void unpackbits();

void reverse_complement(char* s, char* s1);
std::string buildcontig(std::vector<std::string> reads, std::vector<long> pos);

void writecontig(std::string ref,std::vector<long> pos, std::vector<std::string> reads, std::ofstream& f_seq, std::ofstream& f_pos, std::ofstream& f_noise, std::ofstream& f_noisepos);

void getDataParams();

void setglobalarrays();

uint32_t numreads = 0;

int main(int argc, char** argv)
{
	std::string basedir = std::string(argv[1]);
	outfile = basedir + "/output/output.tmp";
	outfile_clean = basedir + "/output/output_clean.dna";
	infile_seq = basedir + "/output/read_seq.txt";
	infile_meta = basedir + "/output/read_meta.txt";
	infile_pos = basedir + "/output/read_pos.txt";
	infile_noise = basedir + "/output/read_noise.txt";
	infile_noisepos = basedir + "/output/read_noisepos.txt";
	infile_rev = basedir + "/output/read_rev.txt";
	infile_order = basedir + "/output/read_order.bin";
	getDataParams(); //populate readlen
	setglobalarrays();
	decode();
	restore_order();
	return 0;
}

void decode()
{
	std::cout << "Decoding reads\n";
	unpackbits();
	std::ofstream f(outfile);
	std::ifstream f_seq(infile_seq);
	std::ifstream f_pos(infile_pos);
	std::ifstream f_noise(infile_noise);
	std::ifstream f_noisepos(infile_noisepos);
	std::ifstream f_rev(infile_rev);


	char currentread[readlen+1],ref[readlen+1],revread[readlen+1];
	currentread[readlen] = '\0';
	revread[readlen] = '\0';
	ref[readlen] = '\0';
	std::string noise;
	char c;
	long pos; 

	while(f_pos >> std::noskipws >> c)//don't skip whitespaces
	{
		numreads++;
		pos = (unsigned char)(c);
		if(pos!=0)
		{
			for(int i = 0; i <= readlen - 1 - pos;  i++)
				ref[i] = ref[i+pos];
			f_seq.get(ref+readlen-pos,pos+1);
		}
		strcpy(currentread,ref);	
		int prevnoisepos = 0,noisepos;
		std::getline(f_noise,noise);
		for(int i = 0; i < noise.size(); i++)
		{
			c = f_noisepos.get();
			noisepos = (unsigned char)c + prevnoisepos;
			currentread[noisepos] = dec_noise[ref[noisepos]][noise[i]];
			prevnoisepos = noisepos;	
		}
		c = f_rev.get();
		if(c == 'd')
			f << currentread<<"\n";
		else
		{
			reverse_complement(currentread,revread);
			f << revread<<"\n";
		}
	}
	f.close();
	f_seq.close();
	f_pos.close();
	f_noise.close();
	f_noisepos.close();
	f_rev.close();
	std::cout<<"Decoding done\n";
	return;
}

void restore_order()
{
	std::cout << "Restoring order\n";
	uint32_t hundred_mil = 100000000;
	std::string s;
	std::ofstream f_clean(outfile_clean);
	for (uint32_t i = 0; i <= numreads/hundred_mil; i++)
	{
		std::ifstream f_order(infile_order,std::ios::binary);
		std::ifstream f(outfile);
		std::ofstream f_bin(outfile+".bin");
		auto numreads_bin = hundred_mil;
		if (i == numreads/hundred_mil)
			numreads_bin = numreads%hundred_mil;
		uint32_t *index_array = new uint32_t [numreads_bin];	
		uint32_t order,pos = 0;
		for(uint32_t j = 0; j < numreads; j++)
		{
			f_order.read((char*)&order,sizeof(uint32_t));
			if (order >= i*hundred_mil && order < i*hundred_mil + numreads_bin)
			{
				index_array[order-i*hundred_mil] = pos;
				f.seekg(uint64_t(j)*(readlen+1), f.beg);
				std::getline(f,s);
				f_bin << s << "\n";
				pos++;
			}
		}
		f_bin.close();
		std::ifstream f_bin_in(outfile+".bin");
		for(uint32_t j = 0; j < numreads_bin; j++)
		{
			f_bin_in.seekg(uint64_t(index_array[j])*(readlen+1), f_bin_in.beg);
			std::getline(f_bin_in,s);
			f_clean << s << "\n";
		}
		delete[] index_array;	
		f_order.close();
		f.close();
		f_bin_in.close();
	}
	f_clean.close();
	return;
}

void unpackbits()
{
	std::ifstream in_seq(infile_seq,std::ios::binary);
	std::ifstream in_noise(infile_noise,std::ios::binary);
	std::ofstream f_seq(infile_seq+".tmp");
	std::ifstream in_seq_tail(infile_seq+".tail");
	std::ofstream f_noise(infile_noise+".tmp");
	std::ifstream in_noise_tail(infile_noise+".tail");
	char inttobase[4];
	inttobase[0] = 'A';
	inttobase[1] = 'C';
	inttobase[2] = 'G';
	inttobase[3] = 'T';
	
	uint8_t dnabin;
	in_seq.read((char*)&dnabin,sizeof(uint8_t));
	while(!in_seq.eof())
	{	
		f_seq << inttobase[dnabin%4];
		dnabin/=4; 	
		f_seq << inttobase[dnabin%4];
		dnabin/=4; 	
		f_seq << inttobase[dnabin%4];
		dnabin/=4; 	
		f_seq << inttobase[dnabin%4];
		in_seq.read((char*)&dnabin,sizeof(uint8_t));
	}
	in_seq.close();
	f_seq << in_seq_tail.rdbuf();
	in_seq_tail.close();
	f_seq.close();
	remove(infile_seq.c_str());
	remove((infile_seq+".tail").c_str());
	rename((infile_seq+".tmp").c_str(),infile_seq.c_str());		
	
	//noise
	inttobase[0] = '0';
	inttobase[1] = '1';
	inttobase[2] = '2';
	inttobase[3] = '\n';
	
	in_noise.read((char*)&dnabin,sizeof(uint8_t));
	while(!in_noise.eof())
	{	
		f_noise << inttobase[dnabin%4];
		dnabin/=4; 	
		f_noise << inttobase[dnabin%4];
		dnabin/=4; 	
		f_noise << inttobase[dnabin%4];
		dnabin/=4; 	
		f_noise << inttobase[dnabin%4];
		in_noise.read((char*)&dnabin,sizeof(uint8_t));
	}
	in_noise.close();
	f_noise << in_noise_tail.rdbuf();
	in_noise_tail.close();
	f_noise.close();
	remove((infile_noise+".tail").c_str());
	rename((infile_noise+".tmp").c_str(),infile_noise.c_str());		
	
	return;
}

void reverse_complement(char* s, char* s1)
{
	for(int j = 0; j < readlen; j++)
		s1[j] = chartorevchar[s[readlen-j-1]];
	return;
}


void setglobalarrays()
{
	dec_noise['A']['0'] = 'C';
	dec_noise['A']['1'] = 'G';
	dec_noise['A']['2'] = 'T';
	dec_noise['A']['3'] = 'N';
	dec_noise['C']['0'] = 'A';
	dec_noise['C']['1'] = 'G';
	dec_noise['C']['2'] = 'T';
	dec_noise['C']['3'] = 'N';
	dec_noise['G']['0'] = 'T';
	dec_noise['G']['1'] = 'A';
	dec_noise['G']['2'] = 'C';
	dec_noise['G']['3'] = 'N';
	dec_noise['T']['0'] = 'G';
	dec_noise['T']['1'] = 'C';
	dec_noise['T']['2'] = 'A';
	dec_noise['T']['3'] = 'N';
	dec_noise['N']['0'] = 'A';
	dec_noise['N']['1'] = 'G';
	dec_noise['N']['2'] = 'C';
	dec_noise['N']['3'] = 'T';
	chartolong['A'] = 0;
	chartolong['C'] = 1;
	chartolong['G'] = 2;
	chartolong['T'] = 3;
	chartolong['N'] = 4;
	chartorevchar['A'] = 'T';
	chartorevchar['C'] = 'G';
	chartorevchar['G'] = 'C';
	chartorevchar['T'] = 'A';
	return;
}
	

void getDataParams()
{
	std::string line;
	std::ifstream myfile(infile_meta, std::ifstream::in);
	
	myfile >> readlen;
	std::cout << "Read length: " << readlen << std::endl;
	myfile.close();
}


