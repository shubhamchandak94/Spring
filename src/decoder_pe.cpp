#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>
#include <string>
#include <bitset>
#include <cstdio>
#include <omp.h>
#include "config.h"

std::string outfile_1;
std::string outfile_2;
std::string outfile_order_2;
std::string infile_quality[2];
std::string infile_id[2];
std::string infile_seq;
std::string infile_meta;
std::string infile_pos;
std::string infile_noise;
std::string infile_noisepos;
std::string infile_rev;
std::string infile_N;
std::string infile_singleton;
std::string infilenumreads;
std::string infile_order_paired;
std::string infile_paired_flag_first;
std::string infile_order;
std::string infile_readlength;

int max_readlen, num_thr, num_thr_e;
uint32_t numreads, numreads_by_2; 
uint32_t global_counter;
uint8_t paired_id_code;
std::string preserve_order, preserve_quality, preserve_id, fast_mode;
uint64_t max_bin_size;

typedef std::bitset<3*MAX_READ_LEN> bitset;

long chartolong[128];
char dec_noise[128][128];
char chartorevchar[128];
char revinttochar[8] = {'A','N','G','#','C','#','T','#'};//used in bitsettostring
bitset basemask[MAX_READ_LEN][128];//bitset for A,G,C,T at each position 
//used in stringtobitset, chartobitset and bitsettostring
bitset positionmask[MAX_READ_LEN];//bitset for each position (1 at two bits and 0 elsewhere)
//used in bitsettostring
bitset mask63;//bitset with 64 bits set to 1 (used in bitsettostring for conversion to ullong)


void decode();

void restore_order(std::string outfile, std::string infile_order, int filenum, uint32_t *index_array, bitset *reads_bin, uint8_t *read_lengths_bin);
//filenum 1 or 2

void modify_id(std::string &id, uint8_t paired_id_code);

void write_fastq_block(char *read, std::string &quality_file_prefix, std::string &id_file_prefix, std::ifstream &fin_quality, std::ifstream &fin_id, int &current_tid_quality, int &current_tid_id, std::ofstream &fout, bool &paired_id_match, int filenum, uint8_t cur_readlen);

void unpackbits();

bitset chartobitset(char *s, uint8_t readlen);

void bitsettostring(bitset b,char *s, uint8_t readlen);

void reverse_complement(char* s, char* s1, uint8_t readlen);

void setglobalarrays();

int main(int argc, char** argv)
{
	std::string basedir = std::string(argv[1]);
	outfile_1 = basedir + "/output_1.txt";
	outfile_2 = basedir + "/output_2.txt";
	outfile_order_2 = basedir + "/read_order_2.bin";
	infile_quality[0] = basedir + "/quality_1.txt";
	infile_quality[1] = basedir + "/quality_2.txt";
	infile_id[0] = basedir + "/id_1.txt";
	infile_id[1] = basedir + "/id_2.txt";
	infile_seq = basedir + "/read_seq.txt";
	infile_meta = basedir + "/read_meta.txt";
	infile_pos = basedir + "/read_pos.txt";
	infile_noise = basedir + "/read_noise.txt";
	infile_noisepos = basedir + "/read_noisepos.txt";
	infile_rev = basedir + "/read_rev.txt";
	infile_N = basedir + "/unaligned_N.txt";
	infile_singleton = basedir + "/unaligned_singleton.txt";
	infilenumreads = basedir + "/numreads.bin";
	infile_order_paired = basedir + "/read_order_paired.bin";
	infile_paired_flag_first = basedir + "/read_paired_flag_first.bin";
	infile_order = basedir + "/read_order.bin";	
	infile_readlength = basedir + "/read_lengths.bin";
	
	max_readlen = atoi(argv[2]);
	num_thr = atoi(argv[3]);
	num_thr_e = atoi(argv[4]);
	preserve_order = std::string(argv[5]);
	preserve_quality = std::string(argv[6]);
	preserve_id = std::string(argv[7]);
	fast_mode = std::string(argv[8]);
		
	std::ifstream f_numreads(infilenumreads, std::ios::binary);
	f_numreads.seekg(4);
	f_numreads.read((char*)&numreads,sizeof(uint32_t));
	f_numreads.read((char*)&paired_id_code,sizeof(uint8_t));
	f_numreads.close();
	numreads_by_2 = numreads/2;

	omp_set_num_threads(num_thr);
	setglobalarrays();
	
	unpackbits();	
	
	decode();
		
	max_bin_size = 12000000000/(sizeof(bitset)+sizeof(uint32_t)+sizeof(uint8_t));
	if(max_bin_size > numreads_by_2 || fast_mode == "True")
		max_bin_size = numreads_by_2 + 1;//+1 just to make sure there is only one pass over the data
	uint32_t *index_array = new uint32_t [max_bin_size];
	uint8_t *read_lengths_bin = new uint8_t [max_bin_size];
	bitset *reads_bin = new bitset [max_bin_size];
	restore_order(outfile_2, outfile_order_2, 2, index_array, reads_bin, read_lengths_bin);
	if(preserve_order == "True")
		restore_order(outfile_1, infile_order, 1, index_array, reads_bin, read_lengths_bin);
	delete[] index_array;
	delete[] reads_bin;
	delete[] read_lengths_bin;
	return 0;
}

void decode()
{
	std::cout << "Decoding reads\n";
	uint32_t numreads_thr[num_thr];
	//first calculate numreads in each thread, this is needed for accessing flag_first, readlength
	#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		numreads_thr[tid] = 0;
		for(int tid_e = tid*num_thr_e/num_thr; tid_e < (tid+1)*num_thr_e/num_thr; tid_e++)
		{
			std::ifstream f_rev(infile_rev+'.'+std::to_string(tid_e), std::ifstream::ate | std::ifstream::binary);	
			numreads_thr[tid] += f_rev.tellg();//size of f_rev file
			f_rev.close();
		}
	}

	#pragma omp parallel
	{
	int tid = omp_get_thread_num();
	uint32_t pos_in_reordered_reads = 0;
	if(tid != 0)
	{
		for(int i = 0; i < tid; i++)
			pos_in_reordered_reads += numreads_thr[i];
	}
	
	for(int tid_e = tid*num_thr_e/num_thr; tid_e < (tid+1)*num_thr_e/num_thr; tid_e++)
	{
		std::ofstream f_1;
		f_1.open(outfile_1+'.'+std::to_string(tid_e),std::ios::binary);
		std::ofstream f_2(outfile_2+'.'+std::to_string(tid_e),std::ios::binary);
		std::ifstream f_seq(infile_seq+'.'+std::to_string(tid_e));
		std::ifstream f_pos(infile_pos+'.'+std::to_string(tid_e));
		std::ifstream f_noise(infile_noise+'.'+std::to_string(tid_e));
		std::ifstream f_noisepos(infile_noisepos+'.'+std::to_string(tid_e));
		std::ifstream f_rev(infile_rev+'.'+std::to_string(tid_e));
		std::ifstream f_readlength(infile_readlength, std::ios::binary);
		std::ifstream f_flag_first(infile_paired_flag_first);
		f_readlength.seekg((uint64_t)pos_in_reordered_reads*sizeof(uint8_t));
		f_flag_first.seekg(pos_in_reordered_reads);
		char currentread[MAX_READ_LEN+1],ref[MAX_READ_LEN+1],revread[MAX_READ_LEN+1];
		bitset b;
		std::string noise;
		char c;
		long pos;
		char flag_first; 
		uint8_t cur_readlen;
		uint8_t ref_len;
		while(f_pos >> std::noskipws >> c)//don't skip whitespaces
		{
			pos = (unsigned char)(c);
			f_readlength.read((char*)&cur_readlen, sizeof(uint8_t));
			f_flag_first >> flag_first;
			if(pos == max_readlen)
			{
				f_seq.get(ref,cur_readlen+1);
				ref_len = cur_readlen;
				ref[ref_len] = '\0';
			}
			else
			{
				if(pos != 0)
					for(int i = 0; i <= ref_len - 1 - pos;  i++)
						ref[i] = ref[i+pos];
				if(cur_readlen > ref_len - pos)
				{
					f_seq.get(ref+ref_len-pos, cur_readlen-ref_len+pos+1);
					ref_len = cur_readlen;
				}
				else
					ref_len = ref_len - pos;
			}
			strncpy(currentread,ref,cur_readlen);
			currentread[cur_readlen] = '\0';
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
				if(flag_first == '1')
				{
					b = chartobitset(currentread,cur_readlen);
					f_1.write((char*)&cur_readlen,sizeof(uint8_t));
					f_1.write((char*)&b,sizeof(bitset));
				}	
				else
				{
					b = chartobitset(currentread,cur_readlen);
					f_2.write((char*)&cur_readlen,sizeof(uint8_t));
					f_2.write((char*)&b,sizeof(bitset));
				}	
			}
			else
			{
				reverse_complement(currentread,revread,cur_readlen);
				if(flag_first == '1')
				{
					b = chartobitset(revread,cur_readlen);
					f_1.write((char*)&cur_readlen,sizeof(uint8_t));
					f_1.write((char*)&b,sizeof(bitset));
				}	
				else
				{
					b = chartobitset(revread,cur_readlen);
					f_2.write((char*)&cur_readlen,sizeof(uint8_t));
					f_2.write((char*)&b,sizeof(bitset));
				}	
			}
			pos_in_reordered_reads++;
		}
		f_1.close();
		f_2.close();
		f_seq.close();
		f_pos.close();
		f_noise.close();
		f_noisepos.close();
		f_rev.close();
		f_readlength.close();
		f_flag_first.close();
	}
	}
	if(preserve_order == "True")
	{
		std::ofstream f_1;	
		f_1.open(outfile_1,std::ios::binary);
		std::ofstream f_2(outfile_2,std::ios::binary);
		for(int tid_e = 0; tid_e < num_thr_e; tid_e++)
		{
			std::ifstream f_in_1(outfile_1+'.'+std::to_string(tid_e),std::ios::binary);
			std::ifstream f_in_2(outfile_2+'.'+std::to_string(tid_e),std::ios::binary);
			f_1 << f_in_1.rdbuf();
			f_2 << f_in_2.rdbuf();
			f_1.clear();//clear error flags if f_in is empty	
			f_2.clear();//clear error flags if f_in is empty	
			f_in_1.close();
			f_in_2.close();
		}
		uint32_t pos_in_reorder_reads = 0;
		for(int i = 0; i < num_thr; i++)
			pos_in_reorder_reads += numreads_thr[i];
		std::ifstream f_readlength(infile_readlength,std::ios::binary);
		std::ifstream f_flag_first(infile_paired_flag_first);
		f_readlength.seekg((uint64_t)pos_in_reorder_reads*sizeof(uint8_t));
		f_flag_first.seekg(pos_in_reorder_reads);
		std::ifstream f_singleton(infile_singleton);
		char currentread[MAX_READ_LEN+1];
		uint8_t cur_readlen;
		char flag_first;
		bitset b;
		f_readlength.read((char*)&cur_readlen,sizeof(uint8_t));
		f_flag_first >> flag_first;
		if(!f_readlength.eof())
		{
			f_singleton.read(currentread,cur_readlen);
			while(!f_singleton.eof())
			{	
				currentread[cur_readlen] = '\0';			
				if(flag_first=='1')
				{
					b = chartobitset(currentread,cur_readlen);
					f_1.write((char*)&cur_readlen,sizeof(uint8_t));
					f_1.write((char*)&b,sizeof(bitset));
				}
				else
				{			
					b = chartobitset(currentread,cur_readlen);
					f_2.write((char*)&cur_readlen,sizeof(uint8_t));
					f_2.write((char*)&b,sizeof(bitset));
				}
				f_readlength.read((char*)&cur_readlen,sizeof(uint8_t));
				f_flag_first >> flag_first;
				f_singleton.read(currentread,cur_readlen);
				pos_in_reorder_reads++;
			}
		}
		f_singleton.close();
		std::ifstream f_N(infile_N);
		if(!f_readlength.eof())
		{
			f_N.read(currentread,cur_readlen);
			while(!f_N.eof())
			{
				currentread[cur_readlen] = '\0';
				if(flag_first=='1')
				{
					b = chartobitset(currentread,cur_readlen);
					f_1.write((char*)&cur_readlen,sizeof(uint8_t));
					f_1.write((char*)&b,sizeof(bitset));
				}
				else
				{			
					b = chartobitset(currentread,cur_readlen);
					f_2.write((char*)&cur_readlen,sizeof(uint8_t));
					f_2.write((char*)&b,sizeof(bitset));
				}
				f_readlength.read((char*)&cur_readlen,sizeof(uint8_t));
				f_flag_first >> flag_first;
				f_N.read(currentread,cur_readlen);
				pos_in_reorder_reads++;
			}
		}
		f_N.close();
		f_1.close();
		f_2.close();
		f_readlength.close();
		f_flag_first.close();
	}
	else
	{
		std::ifstream fin_quality;
		std::ifstream fin_id;
		bool paired_id_match = false;
		std::string quality, id;
		std::string quality_file_prefix, id_file_prefix;
		int current_tid_quality = 0, current_tid_id = 0;
		quality_file_prefix = infile_quality[0];
		id_file_prefix = infile_id[0];
		if(preserve_quality == "True")
			fin_quality.open(quality_file_prefix+"."+std::to_string(current_tid_quality));
		if(preserve_id == "True")
			fin_id.open(id_file_prefix+"."+std::to_string(current_tid_id));
		global_counter = 1;
		char currentread[MAX_READ_LEN+1];
		std::ofstream f_1;	
		f_1.open(outfile_1);
		std::ofstream f_2(outfile_2,std::ios::binary);
		bitset b;
		uint8_t cur_readlen;
		for(int tid_e = 0; tid_e < num_thr_e; tid_e++)
		{
			std::ifstream f_in_1(outfile_1+'.'+std::to_string(tid_e),std::ios::binary);
			std::ifstream f_in_2(outfile_2+'.'+std::to_string(tid_e),std::ios::binary);
			f_in_1.read((char*)&cur_readlen,sizeof(uint8_t));
			while(!f_in_1.eof())
			{
				f_in_1.read((char*)&b,sizeof(bitset));
				bitsettostring(b,currentread,cur_readlen);
				write_fastq_block(currentread,quality_file_prefix,id_file_prefix,fin_quality,fin_id,current_tid_quality,current_tid_id,f_1,paired_id_match, 1, cur_readlen);
				f_in_1.read((char*)&cur_readlen,sizeof(uint8_t));
			}
	
			f_2 << f_in_2.rdbuf();
			f_1.clear();//clear error flags if f_in is empty	
			f_2.clear();//clear error flags if f_in is empty	
			f_in_1.close();
			f_in_2.close();
		}
		uint32_t pos_in_reorder_reads = 0;
		for(int i = 0; i < num_thr; i++)
			pos_in_reorder_reads += numreads_thr[i];
		std::ifstream f_readlength(infile_readlength,std::ios::binary);
		std::ifstream f_flag_first(infile_paired_flag_first);
		f_readlength.seekg((uint64_t)pos_in_reorder_reads*sizeof(uint8_t));
		f_flag_first.seekg(pos_in_reorder_reads);
		std::ifstream f_singleton(infile_singleton);
		char flag_first;
		f_readlength.read((char*)&cur_readlen,sizeof(uint8_t));
		f_flag_first >> flag_first;
		if(!f_readlength.eof())
		{
			f_singleton.read(currentread,cur_readlen);
			while(!f_singleton.eof())
			{	
				currentread[cur_readlen] = '\0';			
				if(flag_first == '1')
				{
					write_fastq_block(currentread,quality_file_prefix,id_file_prefix,fin_quality,fin_id,current_tid_quality,current_tid_id,f_1,paired_id_match, 1, cur_readlen);
				}
				else
				{			
					b = chartobitset(currentread,cur_readlen);
					f_2.write((char*)&cur_readlen,sizeof(uint8_t));
					f_2.write((char*)&b,sizeof(bitset));
				}
				f_readlength.read((char*)&cur_readlen,sizeof(uint8_t));
				f_flag_first >> flag_first; 
				f_singleton.read(currentread,cur_readlen);
				pos_in_reorder_reads++;
			}
		}
		f_singleton.close();
		std::ifstream f_N(infile_N);
		if(!f_readlength.eof())
		{
			f_N.read(currentread,cur_readlen);
			while(!f_N.eof())
			{
				currentread[cur_readlen] = '\0';
				if(flag_first == '1')
				{
					write_fastq_block(currentread,quality_file_prefix,id_file_prefix,fin_quality,fin_id,current_tid_quality,current_tid_id,f_1,paired_id_match, 1, cur_readlen);
				}
				else
				{			
					b = chartobitset(currentread,cur_readlen);
					f_2.write((char*)&cur_readlen,sizeof(uint8_t));
					f_2.write((char*)&b,sizeof(bitset));
				}
				f_readlength.read((char*)&cur_readlen,sizeof(uint8_t));
				f_flag_first >> flag_first; 
				f_N.read(currentread,cur_readlen);
				pos_in_reorder_reads++;
			}
		}
		f_N.close();
		f_1.close();
		f_2.close();
		f_readlength.close();
		f_flag_first.close();
		if(preserve_quality == "True")
			remove((quality_file_prefix+"."+std::to_string(current_tid_quality)).c_str());
	}
	std::cout<<"Decoding done\n";
	return;
}

void restore_order(std::string outfile, std::string infile_order, int filenum, uint32_t *index_array, bitset *reads_bin, uint8_t* read_lengths_bin)
{
	std::cout << "Restoring order\n";
	char s[MAX_READ_LEN+1];
	
	std::ofstream fout(outfile+".tmp");
	std::ifstream fin_quality;
	std::ifstream fin_id;
	bool paired_id_match = false;
	std::string quality_file_prefix, id_file_prefix;
	int current_tid_quality = 0, current_tid_id = 0;
	if(filenum == 1)
	{
		quality_file_prefix = infile_quality[0];
		id_file_prefix = infile_id[0];
	}	
	else
	{
		quality_file_prefix = infile_quality[1];
		if(paired_id_code != 0)
		{
			id_file_prefix = infile_id[0];
			paired_id_match = true;
		}
		else
			id_file_prefix = infile_id[1];
	}
	global_counter = 1;
	if(preserve_quality == "True")
		fin_quality.open(quality_file_prefix+"."+std::to_string(current_tid_quality));
	if(preserve_id == "True")
		fin_id.open(id_file_prefix+"."+std::to_string(current_tid_id));
	for (uint32_t i = 0; i <= numreads_by_2/max_bin_size; i++)
	{
		std::ifstream f_order(infile_order,std::ios::binary);
		std::ifstream f(outfile,std::ios::binary);
		auto numreads_bin = max_bin_size;
		if (i == numreads_by_2/max_bin_size)
			numreads_bin = numreads_by_2%max_bin_size;

		uint32_t order,pos = 0;
		for(uint32_t j = 0; j < numreads_by_2; j++)
		{
			f_order.read((char*)&order,sizeof(uint32_t));
			if (order >= i*max_bin_size && order < i*max_bin_size + numreads_bin)
			{
				index_array[order-i*max_bin_size] = pos;
				f.seekg(uint64_t(j)*(sizeof(bitset)+sizeof(uint8_t)), f.beg);
				f.read((char*)&read_lengths_bin[pos],sizeof(uint8_t));
				f.read((char*)&reads_bin[pos],sizeof(bitset));
				pos++;
			}
		}
		for(uint32_t j = 0; j < numreads_bin; j++)
		{
			bitsettostring(reads_bin[index_array[j]],s,read_lengths_bin[index_array[j]]);
			write_fastq_block(s,quality_file_prefix,id_file_prefix,fin_quality,fin_id,current_tid_quality,current_tid_id,fout,paired_id_match, filenum, read_lengths_bin[index_array[j]]);
		}

		f_order.close();
		f.close();
	}
	if(preserve_quality == "True")
		remove((quality_file_prefix+"."+std::to_string(current_tid_quality)).c_str());
	fout.close();
	remove(outfile.c_str());
	rename((outfile+".tmp").c_str(),outfile.c_str());
	return;
}

void write_fastq_block(char *read, std::string &quality_file_prefix, std::string &id_file_prefix, std::ifstream &fin_quality, std::ifstream &fin_id, int &current_tid_quality, int &current_tid_id, std::ofstream &fout, bool &paired_id_match, int filenum, uint8_t cur_readlen)
{
	std::string id;
	char quality[max_readlen+1];
	if(preserve_quality == "True")
	{
		fin_quality.getline(quality,cur_readlen+1);	
		if(fin_quality.gcount() == 0)
		{
			fin_quality.close();
			remove((quality_file_prefix+"."+std::to_string(current_tid_quality)).c_str());
			current_tid_quality++;
			fin_quality.open(quality_file_prefix+"."+std::to_string(current_tid_quality));
			fin_quality.getline(quality,cur_readlen+1);
		}
		quality[cur_readlen+1] = '\0';
		fin_quality.clear();
	}
	if(preserve_id == "True")
	{
		if(!std::getline(fin_id,id))
		{
			fin_id.close();
			current_tid_id++;
			fin_id.open(id_file_prefix+"."+std::to_string(current_tid_id));
			std::getline(fin_id,id);
		}
		
		if(paired_id_match)
			modify_id(id,paired_id_code);
	}
	if(preserve_id == "False")
		fout << "@" << global_counter++ << "/" << filenum << "\n";
	else
		fout << id << "\n";
	fout << read << "\n";
	if(preserve_quality == "True")
	{
		fout << "+\n";
		fout << quality << "\n";
	}
}

void modify_id(std::string &id, uint8_t paired_id_code)
{
	if(paired_id_code == 2)
		return;
	else if(paired_id_code == 1)
	{
		id.back() = '2';
		return;
	}
	else if(paired_id_code == 3)	
	{
		int i = 0;
		while(id[i] != ' ')
			i++;
		id[i+1] = '2';
		return;
	}
}

void unpackbits()
{
	#pragma omp parallel
	{
	int tid = omp_get_thread_num();
	for(int tid_e = tid*num_thr_e/num_thr; tid_e < (tid+1)*num_thr_e/num_thr; tid_e++)
	{
		std::ofstream f_seq(infile_seq+'.'+std::to_string(tid_e)+".tmp");
		std::ofstream f_rev(infile_rev+'.'+std::to_string(tid_e)+".tmp");
		std::ifstream in_seq(infile_seq+'.'+std::to_string(tid_e),std::ios::binary);
		std::ifstream in_seq_tail(infile_seq+'.'+std::to_string(tid_e)+".tail");
		std::ifstream in_rev(infile_rev+'.'+std::to_string(tid_e),std::ios::binary);
		std::ifstream in_rev_tail(infile_rev+'.'+std::to_string(tid_e)+".tail");
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
		
		//rev
		inttobase[0] = 'd';
		inttobase[1] = 'r';
		
		in_rev.read((char*)&dnabin,sizeof(uint8_t));
		while(!in_rev.eof())
		{	
			f_rev << inttobase[dnabin%2];
			dnabin/=2; 	
			f_rev << inttobase[dnabin%2];
			dnabin/=2; 	
			f_rev << inttobase[dnabin%2];
			dnabin/=2; 	
			f_rev << inttobase[dnabin%2];
			dnabin/=2; 	
			f_rev << inttobase[dnabin%2];
			dnabin/=2; 	
			f_rev << inttobase[dnabin%2];
			dnabin/=2; 	
			f_rev << inttobase[dnabin%2];
			dnabin/=2; 	
			f_rev << inttobase[dnabin%2];
			in_rev.read((char*)&dnabin,sizeof(uint8_t));
		}
		in_rev.close();
		f_rev << in_rev_tail.rdbuf();
		in_rev_tail.close();
		
		f_seq.close();
		f_rev.close();
		remove((infile_seq+'.'+std::to_string(tid_e)).c_str());
		remove((infile_rev+'.'+std::to_string(tid_e)).c_str());
		remove((infile_seq+'.'+std::to_string(tid_e)+".tail").c_str());
		rename((infile_seq+'.'+std::to_string(tid_e)+".tmp").c_str(),(infile_seq+'.'+std::to_string(tid_e)).c_str());
		remove((infile_rev+'.'+std::to_string(tid_e)+".tail").c_str());
		rename((infile_rev+'.'+std::to_string(tid_e)+".tmp").c_str(),(infile_rev+'.'+std::to_string(tid_e)).c_str());	
	}//for end
	}//parallel end
	
	//singleton
	std::ifstream in_singleton(infile_singleton,std::ios::binary);
	std::ofstream f_singleton(infile_singleton+".tmp");
	std::ifstream in_singleton_tail(infile_singleton+".tail");
	char inttobase[4];
	inttobase[0] = 'A';
	inttobase[1] = 'C';
	inttobase[2] = 'G';
	inttobase[3] = 'T';
	
	uint8_t dnabin;
	in_singleton.read((char*)&dnabin,sizeof(uint8_t));
	while(!in_singleton.eof())
	{	
		f_singleton << inttobase[dnabin%4];
		dnabin/=4; 	
		f_singleton << inttobase[dnabin%4];
		dnabin/=4; 	
		f_singleton << inttobase[dnabin%4];
		dnabin/=4; 	
		f_singleton << inttobase[dnabin%4];
		in_singleton.read((char*)&dnabin,sizeof(uint8_t));
	}
	in_singleton.close();
	f_singleton << in_singleton_tail.rdbuf();
	in_singleton_tail.close();
	f_singleton.close();
	remove(infile_singleton.c_str());
	remove((infile_singleton+".tail").c_str());
	rename((infile_singleton+".tmp").c_str(),infile_singleton.c_str());		
	return;
}

void reverse_complement(char* s, char* s1, uint8_t readlen)
{
	for(int j = 0; j < readlen; j++)
		s1[j] = chartorevchar[s[readlen-j-1]];
	s1[readlen] = '\0';
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
	chartorevchar['N'] = 'N';
	for(int i = 0; i < 63; i++)
		mask63[i] = 1;
	for(int i = 0; i < max_readlen; i++)
	{
		basemask[i]['A'][3*i] = 0;
		basemask[i]['A'][3*i+1] = 0;
		basemask[i]['A'][3*i+2] = 0;
		basemask[i]['C'][3*i] = 0;
		basemask[i]['C'][3*i+1] = 0;
		basemask[i]['C'][3*i+2] = 1;
		basemask[i]['G'][3*i] = 0;
		basemask[i]['G'][3*i+1] = 1;
		basemask[i]['G'][3*i+2] = 0;
		basemask[i]['T'][3*i] = 0;
		basemask[i]['T'][3*i+1] = 1;
		basemask[i]['T'][3*i+2] = 1;
		basemask[i]['N'][3*i] = 1;
		basemask[i]['N'][3*i+1] = 0;
		basemask[i]['N'][3*i+2] = 0;
		positionmask[i][3*i] = 1;
		positionmask[i][3*i+1] = 1;
		positionmask[i][3*i+2] = 1;
	}		
	return;
}
	

void bitsettostring(bitset b,char *s, uint8_t readlen)
{
	unsigned long long ull;
	for(int i = 0; i < (3*readlen)/63+1; i++)
	{	
		ull = (b&mask63).to_ullong();
		b>>=63;
		for(int j = 21*i  ; j < 21*i+21 && j < readlen ; j++)
		{
			s[j] = revinttochar[ull%8];	
			ull/=8;
		}
	}
	s[readlen] = '\0';
	return;
}

bitset chartobitset(char *s, uint8_t readlen)
{
	bitset b;
	for(int i = 0; i < readlen; i++)
		b |= basemask[i][s[i]];
	return b;
}
