#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <cstring>
#include <string>
#include <numeric>
#include <cstdio>
#include <omp.h>
#include "config.h"

uint32_t numreads;

std::string infile;
std::string infile_flag;
std::string infile_pos;
std::string infile_seq;
std::string outfile_seq;
std::string outfile_meta;
std::string outfile_pos;
std::string outfile_noise;
std::string outfile_noisepos;


char longtochar[] = {'A','C','G','T'};
long chartolong[128];
char enc_noise[128][128];

void encode();

void packbits();

std::string buildcontig(std::vector<std::string> reads, std::vector<long> pos);

void writecontig(std::string ref,std::vector<long> pos, std::vector<std::string> reads, std::ofstream& f_seq, std::ofstream& f_pos, std::ofstream& f_noise, std::ofstream& f_noisepos);

void getDataParams();

void setglobalarrays();

int main(int argc, char** argv)
{
	std::string basedir = std::string(argv[1]);
	infile = basedir + "/output/temp.dna";
	infile_pos = basedir + "/output/temppos.txt";
	infile_flag = basedir + "/output/tempflag.txt";
		
	outfile_seq = basedir + "/output/read_seq.txt";
	outfile_meta = basedir + "/output/read_meta.txt";
	outfile_pos = basedir + "/output/read_pos.txt";
	outfile_noise = basedir + "/output/read_noise.txt";
	outfile_noisepos = basedir + "/output/read_noisepos.txt";
	omp_set_num_threads(num_thr);
	getDataParams(); //populate readlen
	setglobalarrays();
	encode();
	return 0;
}

void encode()
{
	std::cout<<"Encoding reads\n";
	#pragma omp parallel 
	{
	int tid = omp_get_thread_num();
	std::ifstream f(infile);
	std::ifstream in_flag(infile_flag);
	std::ifstream in_pos(infile_pos,std::ios::binary);
	std::ofstream f_seq(outfile_seq+'.'+std::to_string(tid));
	std::ofstream f_pos(outfile_pos+'.'+std::to_string(tid));
	std::ofstream f_noise(outfile_noise+'.'+std::to_string(tid));
	std::ofstream f_noisepos(outfile_noisepos+'.'+std::to_string(tid));
	
	uint64_t i, stop;	
	//doing initial setup and first read
	i = uint64_t(tid)*numreads/omp_get_num_threads();//spread out first read equally
	stop = uint64_t(tid+1)*numreads/omp_get_num_threads();
	if(tid == omp_get_num_threads()-1)
		stop = numreads;
	f.seekg(uint64_t(i)*(readlen+1), f.beg);
	in_flag.seekg(i, in_flag.beg);
	in_pos.seekg(i*sizeof(uint8_t),in_pos.beg);
	std::string current,ref;
	char c;
	std::vector<std::string> reads;
	std::vector<long> pos;
	uint8_t p;
	while(i < stop)
	{
		std::getline(f,current);
		c = in_flag.get();
		in_pos.read((char*)&p,sizeof(uint8_t));
		if(c=='0'||reads.size()>10000000)//so that reads vector doesn't get too large
		{
			if(reads.size()!=0)
			{
				ref = buildcontig(reads,pos);
				writecontig(ref,pos,reads,f_seq,f_pos,f_noise,f_noisepos);
			}
			reads = {current};
			pos = {p};
		}
		else
		{
			reads.push_back(current);
			pos.push_back(p);
		}
		i++;	
					
	}
	ref = buildcontig(reads,pos);
	writecontig(ref,pos,reads,f_seq,f_pos,f_noise,f_noisepos);

	f.close();
	in_flag.close();
	in_pos.close();
	f_seq.close();
	f_pos.close();
	f_noise.close();
	f_noisepos.close();
	}
	
	//Combine files produced by the threads
	std::ofstream f_seq(outfile_seq);
	std::ofstream f_pos(outfile_pos);
	std::ofstream f_noise(outfile_noise);
	std::ofstream f_noisepos(outfile_noisepos);
	std::ofstream f_meta(outfile_meta);
	for(int tid = 0; tid < num_thr; tid++)
	{
		std::ifstream in_seq(outfile_seq+'.'+std::to_string(tid));
		std::ifstream in_pos(outfile_pos+'.'+std::to_string(tid));
		std::ifstream in_noise(outfile_noise+'.'+std::to_string(tid));
		std::ifstream in_noisepos(outfile_noisepos+'.'+std::to_string(tid));
		f_seq << in_seq.rdbuf();
		f_pos << in_pos.rdbuf();
		f_noise << in_noise.rdbuf();
		f_noisepos << in_noisepos.rdbuf();
		remove((outfile_seq+'.'+std::to_string(tid)).c_str());
		remove((outfile_pos+'.'+std::to_string(tid)).c_str());
		remove((outfile_noise+'.'+std::to_string(tid)).c_str());
		remove((outfile_noisepos+'.'+std::to_string(tid)).c_str());
	}
	f_meta << readlen << "\n";
	f_seq.close();
	f_meta.close();
	f_pos.close();
	f_noise.close();
	f_noisepos.close();
	packbits();
	std::cout << "Encoding done\n";
	return;
}

void packbits()
{
	std::ifstream in_seq(outfile_seq);
	std::ifstream in_noise(outfile_noise);
	std::ofstream f_seq(outfile_seq+".tmp",std::ios::binary);
	std::ofstream f_seq_tail(outfile_seq+".tail");
	std::ofstream f_noise(outfile_noise+".tmp",std::ios::binary);
	std::ofstream f_noise_tail(outfile_noise+".tail");
	uint64_t file_len=0;
	char c;
	while(in_seq >> std::noskipws >> c)
		file_len++;
	uint8_t basetoint[128];
	basetoint['A'] = 0;
	basetoint['C'] = 1;
	basetoint['G'] = 2;
	basetoint['T'] = 3;
	
	in_seq.close();
	in_seq.open(outfile_seq);
	char dnabase[4];
	uint8_t dnabin;
	for(uint64_t i = 0; i < file_len/4; i++)
	{
		in_seq.read(dnabase,4);
		
		dnabin = 64*basetoint[dnabase[3]]+16*basetoint[dnabase[2]]+4*
			basetoint[dnabase[1]]+basetoint[dnabase[0]];
		f_seq.write((char*)&dnabin,sizeof(uint8_t));
	}
	f_seq.close();
	in_seq.read(dnabase,file_len%4);
	for(int i=0; i<file_len%4;i++)
		f_seq_tail << dnabase[i];
	f_seq_tail.close();
	in_seq.close();
	remove(outfile_seq.c_str());
	rename((outfile_seq+".tmp").c_str(),outfile_seq.c_str());		
	
	//noise
	file_len=0;
	while(in_noise >> std::noskipws >> c)
		file_len++;
	basetoint['0'] = 0;
	basetoint['1'] = 1;
	basetoint['2'] = 2;
	basetoint['\n'] = 3;
	
	in_noise.close();
	in_noise.open(outfile_noise);
	for(uint64_t i = 0; i < file_len/4; i++)
	{
		in_noise.read(dnabase,4);
		dnabin = 64*basetoint[dnabase[3]]+16*basetoint[dnabase[2]]+4*
			basetoint[dnabase[1]]+basetoint[dnabase[0]];
		f_noise.write((char*)&dnabin,sizeof(uint8_t));
	}
	f_noise.close();
	in_noise.read(dnabase,file_len%4);
	for(int i=0; i<file_len%4;i++)
		f_noise_tail << dnabase[i];
	f_noise_tail.close();
	remove(outfile_noise.c_str());
	rename((outfile_noise+".tmp").c_str(),outfile_noise.c_str());
	
	return;
}


std::string buildcontig(std::vector<std::string> reads, std::vector<long> pos)
{
	if(reads.size() == 1)
		return reads[0];
	std::vector<std::array<long,4>> count(readlen,{0,0,0,0});
	for(long i = 0; i < readlen; i++)
		count[i][chartolong[reads[0][i]]] = 1;
	long prevpos = 0,currentpos;
	for(long j = 1; j < reads.size(); j++)
	{
		count.insert(count.end(),pos[j],{0,0,0,0});
		currentpos = prevpos + pos[j];
		for(long i = 0; i < readlen; i++)
			count[currentpos+i][chartolong[reads[j][i]]] += 1;
		prevpos = currentpos;
	}
	std::string ref(count.size(),'A');
	for(long i = 0; i < count.size(); i++)
	{
		long max = 0,indmax = 0;
		for(long j = 0; j < 4; j++)
			if(count[i][j]>max)
			{
				max = count[i][j];
				indmax = j;
			}
		ref[i] = longtochar[indmax];
	}
	return ref;
}

void writecontig(std::string ref,std::vector<long> pos, std::vector<std::string> reads, std::ofstream& f_seq, std::ofstream& f_pos, std::ofstream& f_noise, std::ofstream& f_noisepos)
{
	f_seq << ref;
	char c;
	if(reads.size() == 1)
	{
		f_noise << "\n";
		c = readlen;//(not pos[0] to handle breaks in read sequence due to limit on reads.size() - can't
				//assume  pos[0] = readlen)
		f_pos << c;
		return;
	}
	long prevj = 0;
	for(long j = 0; j < readlen; j++)
		if(reads[0][j] != ref[j])
		{
			f_noise<<enc_noise[ref[j]][reads[0][j]];
			c = j-prevj;
			f_noisepos<<c;
			prevj = j;
		}
	f_noise << "\n";
	c = readlen;// (to handle breaks in read sequence due to limit on reads.size()
	f_pos << c;
	long prevpos = 0,currentpos;
	for(long i = 1; i < reads.size(); i++)
	{
		currentpos = prevpos + pos[i];
		prevj = 0;
		for(long j = 0; j < readlen; j++)
			if(reads[i][j] != ref[currentpos+j])
			{
				f_noise<<enc_noise[ref[currentpos+j]][reads[i][j]];
				c = j-prevj;
				f_noisepos<<c;
				prevj = j;
			}
		f_noise << "\n";
		c = pos[i];
		f_pos << c;
		prevpos = currentpos;
	}
	return;
}

void setglobalarrays()
{
	//enc_noise uses substitution statistics from Minoche et al. 
	enc_noise['A']['C'] = '0';
	enc_noise['A']['G'] = '1';
	enc_noise['A']['T'] = '2';
	enc_noise['A']['N'] = '3';
	enc_noise['C']['A'] = '0';
	enc_noise['C']['G'] = '1';
	enc_noise['C']['T'] = '2';
	enc_noise['C']['N'] = '3';
	enc_noise['G']['T'] = '0';
	enc_noise['G']['A'] = '1';
	enc_noise['G']['C'] = '2';
	enc_noise['G']['N'] = '3';
	enc_noise['T']['G'] = '0';
	enc_noise['T']['C'] = '1';
	enc_noise['T']['A'] = '2';
	enc_noise['T']['N'] = '3';
	enc_noise['N']['A'] = '0';
	enc_noise['N']['G'] = '1';
	enc_noise['N']['C'] = '2';
	enc_noise['N']['T'] = '3';
	chartolong['A'] = 0;
	chartolong['C'] = 1;
	chartolong['G'] = 2;
	chartolong['T'] = 3;
	chartolong['N'] = 4;
	return;
}
	

void getDataParams()
{
	uint32_t number_of_lines = 0; 
	std::string line;
	std::ifstream myfile(infile, std::ifstream::in);

	while (std::getline(myfile, line))
		++number_of_lines;
	//readlen = read_length;
	numreads = number_of_lines;
	std::cout << "Read length: " << readlen << std::endl;
	std::cout << "Number of reads: " << number_of_lines << std::endl;
	myfile.close();
}
