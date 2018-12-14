#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>
#include <string>
#include <bitset>
#include "config.h"

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
char revinttochar[4] = {'A','G','C','T'};//used in bitsettostring
std::bitset<2*readlen> basemask[readlen][128];//bitset for A,G,C,T at each position 
//used in stringtobitset, chartobitset and bitsettostring
std::bitset<2*readlen> positionmask[readlen];//bitset for each position (1 at two bits and 0 elsewhere)
//used in bitsettostring
std::bitset<2*readlen> mask64;//bitset with 64 bits set to 1 (used in bitsettostring for conversion to ullong)


void decode();

void restore_order();

void unpackbits();

std::bitset<2*readlen> stringtobitset(std::string s);

std::bitset<2*readlen> chartobitset(char *s);

void bitsettostring(std::bitset<2*readlen> b,char *s);

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
//	getDataParams(); //populate readlen
	setglobalarrays();
	decode();
	restore_order();
	return 0;
}

void decode()
{
	std::cout << "Decoding reads\n";
	unpackbits();
	std::ofstream f(outfile,std::ios::binary);
	std::ifstream f_seq(infile_seq);
	std::ifstream f_pos(infile_pos);
	std::ifstream f_noise(infile_noise);
	std::ifstream f_noisepos(infile_noisepos);
	std::ifstream f_rev(infile_rev);


	char currentread[readlen+1],ref[readlen+1],revread[readlen+1];
	std::bitset<2*readlen> b;
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
		{
			b = chartobitset(currentread);
			f.write((char*)&b,sizeof(std::bitset<2*readlen>));
		}
		else
		{
			reverse_complement(currentread,revread);
			b = chartobitset(revread);
			f.write((char*)&b,sizeof(std::bitset<2*readlen>));
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
	uint64_t max_bin_size;
	if(MAX_BIN_SIZE <= 3)
		max_bin_size = uint64_t(3)*200000000/7;
	else 
		max_bin_size = uint64_t(MAX_BIN_SIZE)*200000000/7;
	char s[readlen+1];
	s[readlen] = '\0';
	std::ofstream f_clean(outfile_clean);
	for (uint32_t i = 0; i <= numreads/max_bin_size; i++)
	{
		std::ifstream f_order(infile_order,std::ios::binary);
		std::ifstream f(outfile,std::ios::binary);
		auto numreads_bin = max_bin_size;
		if (i == numreads/max_bin_size)
			numreads_bin = numreads%max_bin_size;
		uint32_t *index_array = new uint32_t [numreads_bin];
		std::bitset<2*readlen> *reads_bin = new std::bitset<2*readlen> [numreads_bin];	
		uint32_t order,pos = 0;
		for(uint32_t j = 0; j < numreads; j++)
		{
			f_order.read((char*)&order,sizeof(uint32_t));
			if (order >= i*max_bin_size && order < i*max_bin_size + numreads_bin)
			{
				index_array[order-i*max_bin_size] = pos;
				f.seekg(uint64_t(j)*sizeof(std::bitset<2*readlen>), f.beg);
				f.read((char*)&reads_bin[pos],sizeof(std::bitset<2*readlen>));
				pos++;
			}
		}
		for(uint32_t j = 0; j < numreads_bin; j++)
		{
			bitsettostring(reads_bin[index_array[j]],s);
			f_clean << s << "\n";
		}
		delete[] index_array;
		delete[] reads_bin;
		f_order.close();
		f.close();
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
	for(int i = 0; i < 64; i++)
		mask64[i] = 1;
	for(int i = 0; i < readlen; i++)
	{
		basemask[i]['A'][2*i] = 0;
		basemask[i]['A'][2*i+1] = 0;
		basemask[i]['C'][2*i] = 0;
		basemask[i]['C'][2*i+1] = 1;
		basemask[i]['G'][2*i] = 1;
		basemask[i]['G'][2*i+1] = 0;
		basemask[i]['T'][2*i] = 1;
		basemask[i]['T'][2*i+1] = 1;
		positionmask[i][2*i] = 1;
		positionmask[i][2*i+1] = 1;
	}		
	return;
}
	
/*
void getDataParams()
{
	std::string line;
	std::ifstream myfile(infile_meta, std::ifstream::in);
	
	myfile >> readlen;
	std::cout << "Read length: " << readlen << std::endl;
	myfile.close();
}
*/

std::bitset<2*readlen> stringtobitset(std::string s)
{
	std::bitset<2*readlen> b;
	for(int i = 0; i < readlen; i++)
		b |= basemask[i][s[i]];
	return b;
}

void bitsettostring(std::bitset<2*readlen> b,char *s)
{
	unsigned long long ull,rem;
	for(int i = 0; i < (2*readlen)/64+1; i++)
	{	
		ull = (b&mask64).to_ullong();
		b>>=64;
		for(int j = 32*i  ; j < 32*i+32 && j < readlen ; j++)
		{
			s[j] = revinttochar[ull%4];	
			ull/=4;
		}
	}
	return;
}

std::bitset<2*readlen> chartobitset(char *s)
{
	std::bitset<2*readlen> b;
	for(int i = 0; i < readlen; i++)
		b |= basemask[i][s[i]];
	return b;
}

