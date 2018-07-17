#include <iostream>
#include <fstream>
#include <bitset>
#include <string>
#include <vector>
#include <array>
#include <list>
#include <algorithm>
#include <cstring>
#include <string>
#include <numeric>
#include <cstdio>
#include <cstdlib>
#include <omp.h>
#include "BooPHF.h"
#include "config.h"

uint32_t numreads, numreads_s, numreads_N;
int numdict_s = 2;
uint32_t dict_start[2];
uint32_t dict_end[2];
int thresh_s = 24;
int maxsearch = 1000;
int max_readlen, num_thr;


typedef boomphf::SingleHashFunctor<u_int64_t>  hasher_t;
typedef boomphf::mphf<  u_int64_t, hasher_t  > boophf_t;
typedef std::bitset<3*MAX_READ_LEN> bitset;

class bbhashdict
{
	public:
	boophf_t * bphf;
	uint32_t numkeys;
	uint32_t dict_numreads;//number of reads in this dict (for variable length)
	uint32_t *startpos;
	uint32_t *read_id;
	bool *empty_bin;
	void findpos(int64_t *dictidx, uint32_t &startposidx);
	void remove(int64_t *dictidx, uint32_t &startposidx, uint32_t current);
	bbhashdict()
	{
		bphf = NULL;
		startpos = NULL;
		read_id = NULL;
	}
	~bbhashdict()
	{
		delete[] startpos;
		delete[] read_id;
		delete bphf;
	}	
};

struct contig_reads
{
	std::string read;
	int64_t pos;
	char RC;
	uint32_t order;
	uint8_t read_length;
};

std::string basedir;
std::string infile;
std::string infile_flag;
std::string infile_pos;
std::string infile_seq;
std::string infile_RC;
std::string infile_readlength;
std::string infile_N;
std::string outfile_N;
std::string outfile_seq;
std::string outfile_pos;
std::string outfile_noise;
std::string outfile_noisepos;
std::string outfile_singleton;
std::string infile_order;
std::string infile_order_N;
std::string infilenumreads;

char longtochar[] = {'A','C','G','T','N'};
long chartolong[128];
char enc_noise[128][128];


//Some global arrays (some initialized in setglobalarrays())
char revinttochar[8] = {'A','N','G','#','C','#','T','#'};//used in bitsettostring
char inttochar[] = {'A','C','G','T','N'};
char chartorevchar[128];//A-T etc for reverse complement
bitset basemask[MAX_READ_LEN][128];//bitset for A,G,C,T,N at each position 
//used in stringtobitset, chartobitset and bitsettostring
bitset positionmask[MAX_READ_LEN];//bitset for each position (1 at two bits and 0 elsewhere)
//used in bitsettostring
bitset mask63;//bitset with 63 bits set to 1 (used in bitsettostring for conversion to ullong)

void stringtobitset(std::string s, uint8_t readlen, bitset &b);

void correct_order(uint32_t *order_s);
//correct singleton and clean read order to their actual position in fastq (as opposed to pos in clean reads)

void generatemasks(bitset mask[MAX_READ_LEN][MAX_READ_LEN]);

void readsingletons(bitset *read, uint32_t *order_s, uint8_t *read_lengths_s);

void constructdictionary(bitset *read, bbhashdict *dict, uint8_t *read_lengths_s);

std::string bitsettostring(bitset b, uint8_t readlen);

std::string reverse_complement(std::string s, uint8_t readlen);

void writetofile(bitset *read);

void generateindexmasks(bitset *mask1);

void encode(bitset *read, bbhashdict *dict, uint32_t *order_s, uint8_t *read_lengths_s);

void packbits();

std::string buildcontig(std::list<contig_reads> &current_contig, uint32_t list_size);//using lists to enable faster insertion of singletons

void writecontig(std::string &ref,std::list<contig_reads> &current_contig, std::ofstream& f_seq, std::ofstream& f_pos, std::ofstream& f_noise, std::ofstream& f_noisepos, std::ofstream& f_order, std::ofstream &f_RC, std::ofstream &f_readlength);

void getDataParams();

void setglobalarrays();

int main(int argc, char** argv)
{
	basedir = std::string(argv[1]);
	infile = basedir + "/temp.dna";
	infile_pos = basedir + "/temppos.txt";
	infile_flag = basedir + "/tempflag.txt";
	infile_order = basedir + "/read_order.bin";
	infile_order_N = basedir + "/read_order_N.bin";
	infile_RC = basedir + "/read_rev.txt";
	infile_readlength = basedir + "/read_lengths.bin";
	infile_N = basedir + "/input_N.dna";		
	outfile_N = basedir + "/unaligned_N.txt";
	outfile_seq = basedir + "/read_seq.txt";
	outfile_pos = basedir + "/read_pos.txt";
	outfile_noise = basedir + "/read_noise.txt";
	outfile_noisepos = basedir + "/read_noisepos.txt";
	outfile_singleton = basedir + "/unaligned_singleton.txt";
	infilenumreads = basedir + "/numreads.bin";	

	max_readlen = atoi(argv[2]);
	num_thr = atoi(argv[3]);

	omp_set_num_threads(num_thr);
	getDataParams(); //populate numreads
	setglobalarrays();
	bitset *read = new bitset [numreads_s+numreads_N];
	uint32_t *order_s = new uint32_t [numreads_s+numreads_N];
	uint8_t *read_lengths_s = new uint8_t [numreads_s+numreads_N];
	readsingletons(read,order_s,read_lengths_s);
	correct_order(order_s);

	if(max_readlen > 50)
	{
		dict_start[0] = 0;
		dict_end[0] = 20;
		dict_start[1] = 21;
		dict_end[1] = 41;
	}
	else
	{
		dict_start[0] = 0;
		dict_end[0] = 20*max_readlen/50;
		dict_start[1] = 20*max_readlen/50 + 1;
		dict_end[1] = 41*max_readlen/50;
	}
	bbhashdict *dict = new bbhashdict[numdict_s];
	constructdictionary(read,dict,read_lengths_s);
	encode(read,dict,order_s,read_lengths_s);
	delete[] read;
	delete[] dict;
	return 0;
}

void encode(bitset *read, bbhashdict *dict, uint32_t *order_s, uint8_t *read_lengths_s)
{
	omp_lock_t *read_lock = new omp_lock_t [numreads_s+numreads_N];
	omp_lock_t *dict_lock = new omp_lock_t [numreads_s+numreads_N];
	for(int j = 0; j < numreads_s+numreads_N; j++)
	{
		omp_init_lock(&read_lock[j]);
		omp_init_lock(&dict_lock[j]);
	}
	bool *remainingreads = new bool [numreads_s+numreads_N];
	std::fill(remainingreads, remainingreads+numreads_s+numreads_N,1);

	bitset mask1[numdict_s];
	generateindexmasks(mask1);
	bitset mask[MAX_READ_LEN][MAX_READ_LEN];
	generatemasks(mask);	
	std::cout<<"Encoding reads\n";
	#pragma omp parallel 
	{
	int tid = omp_get_thread_num();
	std::ifstream f(infile+'.'+std::to_string(tid));
	std::ifstream in_flag(infile_flag+'.'+std::to_string(tid));
	std::ifstream in_pos(infile_pos+'.'+std::to_string(tid),std::ios::binary);
	std::ifstream in_order(infile_order+'.'+std::to_string(tid),std::ios::binary);
	std::ifstream in_RC(infile_RC+'.'+std::to_string(tid));
	std::ifstream in_readlength(infile_readlength+'.'+std::to_string(tid),std::ios::binary);
	std::ofstream f_seq(outfile_seq+'.'+std::to_string(tid));
	std::ofstream f_pos(outfile_pos+'.'+std::to_string(tid));
	std::ofstream f_noise(outfile_noise+'.'+std::to_string(tid));
	std::ofstream f_noisepos(outfile_noisepos+'.'+std::to_string(tid));
	std::ofstream f_order(infile_order+'.'+std::to_string(tid)+".tmp",std::ios::binary);
	std::ofstream f_RC(infile_RC+'.'+std::to_string(tid)+".tmp");
	std::ofstream f_readlength(infile_readlength+'.'+std::to_string(tid)+".tmp",std::ios::binary);
	int64_t dictidx[2];//to store the start and end index (end not inclusive) in the dict read_id array
	uint32_t startposidx;//index in startpos
	uint64_t ull;
	bool flag = 0;
	//flag to check if match was found or not
	std::string current,ref;
	bitset forward_bitset,reverse_bitset,b;
	char c,rc;
	std::list<contig_reads> current_contig;
	int64_t p;
	uint8_t rl;
	uint32_t ord, list_size = 0;//list_size variable introduced because list::size() was running very slowly
					// on UIUC machine
	std::list<uint32_t> deleted_rids[numdict_s];
	bool done = false;
	while(!done)
	{
		if(!(in_flag >> c))
			done = true;
		if(!done)
		{
			std::getline(f,current);
			rc = in_RC.get();
			in_pos.read((char*)&p,sizeof(int64_t));
			in_order.read((char*)&ord,sizeof(uint32_t));
			in_readlength.read((char*)&rl,sizeof(uint8_t));
		}
		if(c=='0'|| done || list_size > 10000000)//limit on list size so that memory doesn't get too large
		{
			if(list_size!=0)
			{
				//sort contig according to pos
				current_contig.sort([](const contig_reads &a, const contig_reads &b) 
					{ 
					    return a.pos < b.pos; 
					});				
				//make first pos zero and shift all pos values accordingly
				auto current_contig_it = current_contig.begin();
				int64_t first_pos = (*current_contig_it).pos;
				for(;current_contig_it != current_contig.end(); ++current_contig_it)
					(*current_contig_it).pos -= first_pos;

				ref = buildcontig(current_contig,list_size);
				//try to align the singleton reads to ref
				//first create bitsets from first readlen positions of ref
				forward_bitset.reset();
				reverse_bitset.reset();
				if(ref.size() >= max_readlen)
				{
				stringtobitset(ref.substr(0,max_readlen),max_readlen,forward_bitset);
				stringtobitset(reverse_complement(ref.substr(0,max_readlen),max_readlen),max_readlen,reverse_bitset);
				for(long j = 0; j < ref.size()-max_readlen+1; j++)
				{
					//search for singleton reads
					for(int rev = 0; rev < 2; rev++)
					{
						for(int l = 0; l < numdict_s; l++)
						{
							if(!rev)
								b = forward_bitset&mask1[l];
							else
								b = reverse_bitset&mask1[l];
							ull = (b>>3*dict_start[l]).to_ullong();
							startposidx = dict[l].bphf->lookup(ull);
							if(startposidx >= dict[l].numkeys)//not found
								continue;
							//check if any other thread is modifying same dictpos
							if(!omp_test_lock(&dict_lock[startposidx]))
								continue;
							dict[l].findpos(dictidx,startposidx);
							if(dict[l].empty_bin[startposidx])//bin is empty
							{
								omp_unset_lock(&dict_lock[startposidx]);
								continue;
							}
							uint64_t ull1 = ((read[dict[l].read_id[dictidx[0]]]&mask1[l])>>3*dict_start[l]).to_ullong();
							if(ull == ull1)//checking if ull is actually the key for this bin
							{	
								for (int64_t i = dictidx[1] - 1 ; i >= dictidx[0] && i >= dictidx[1] - maxsearch; i--)
								{
									auto rid = dict[l].read_id[i];
									int hamming;
									if(!rev)
										hamming = ((forward_bitset^read[rid])&mask[0][max_readlen-read_lengths_s[rid]]).count();
									else
										hamming = ((reverse_bitset^read[rid])&mask[0][max_readlen-read_lengths_s[rid]]).count();
									if(hamming<=thresh_s)
									{	
										omp_set_lock(&read_lock[rid]);
										if(remainingreads[rid])
										{
											remainingreads[rid]=0;
											flag = 1;
										}
										omp_unset_lock(&read_lock[rid]);
									}
									if(flag == 1)//match found
									{
										flag = 0;
										list_size++;
										char rc = rev?'r':'d';
										long pos = rev?(j+max_readlen-read_lengths_s[rid]):j;
										std::string read_string = rev?reverse_complement(bitsettostring(read[rid],read_lengths_s[rid]),read_lengths_s[rid]):bitsettostring(read[rid],read_lengths_s[rid]);
										current_contig.push_back({read_string,pos,rc,order_s[rid],read_lengths_s[rid]});								
										for(int l1 = 0;l1 < numdict_s; l1++)
										{	if(read_lengths_s[rid] > dict_end[l1])
												deleted_rids[l1].push_back(rid);
										}
									}
								}
							}
							omp_unset_lock(&dict_lock[startposidx]);
							//delete from dictionaries
							for(int l1= 0; l1 < numdict_s; l1++)
								for(auto it = deleted_rids[l1].begin(); it!=deleted_rids[l1].end();)
								{
									b = read[*it]&mask1[l1];
									ull = (b>>3*dict_start[l1]).to_ullong();
									startposidx = dict[l1].bphf->lookup(ull);
									if(!omp_test_lock(&dict_lock[startposidx]))
									{
										++it;
										continue;
									}
									dict[l1].findpos(dictidx,startposidx);
									dict[l1].remove(dictidx,startposidx,*it);
									it = deleted_rids[l1].erase(it);	
									omp_unset_lock(&dict_lock[startposidx]);
								}
						}
					}
					if(j != ref.size()-max_readlen)//not at last position,shift bitsets
					{
						forward_bitset >>= 3;
						forward_bitset = forward_bitset & mask[0][0];
						forward_bitset |= basemask[max_readlen-1][ref[j+max_readlen]];
						reverse_bitset <<= 3;
						reverse_bitset = reverse_bitset & mask[0][0];
						reverse_bitset |= basemask[0][chartorevchar[ref[j+max_readlen]]];
					}	
							
				}//end for
				}//end if
				//sort contig according to pos
				current_contig.sort([](const contig_reads &a, const contig_reads &b) 
					{ 
					    return a.pos < b.pos; 
					});				
				writecontig(ref,current_contig,f_seq,f_pos,f_noise,f_noisepos,f_order,f_RC,f_readlength);
			}
			if(!done)
			{
				current_contig = {{current,p,rc,ord,rl}};
				list_size = 1;
			}
		}
		else if(c == '1') //read found during rightward search
		{
			current_contig.push_back({current,p,rc,ord,rl});
			list_size++;
		}
	}
	f.close();
	in_flag.close();
	in_pos.close();
	in_order.close();
	in_RC.close();
	in_readlength.close();
	f_seq.close();
	f_pos.close();
	f_noise.close();
	f_noisepos.close();
	f_order.close();
	f_readlength.close();
	f_RC.close();
	}

	//Combine files produced by the threads
	std::ofstream f_order(infile_order);
	std::ofstream f_readlength(infile_readlength);
	for(int tid = 0; tid < num_thr; tid++)
	{
		std::ifstream in_order(infile_order+'.'+std::to_string(tid)+".tmp");
		std::ifstream in_readlength(infile_readlength+'.'+std::to_string(tid)+".tmp");
		f_order << in_order.rdbuf();
		f_order.clear();//clear error flag in case in_order is empty
		f_readlength << in_readlength.rdbuf();
		f_readlength.clear();//clear error flag in case in_readlength is empty
		
		remove((infile_order+'.'+std::to_string(tid)).c_str());
		remove((infile_order+'.'+std::to_string(tid)+".tmp").c_str());
		remove((infile_readlength+'.'+std::to_string(tid)).c_str());
		remove((infile_readlength+'.'+std::to_string(tid)+".tmp").c_str());
		rename((infile_RC+'.'+std::to_string(tid)+".tmp").c_str(),(infile_RC+'.'+std::to_string(tid)).c_str());
		remove((infile_RC+'.'+std::to_string(tid)+".tmp").c_str());
		remove((infile_flag+'.'+std::to_string(tid)).c_str());
		remove((infile_pos+'.'+std::to_string(tid)).c_str());
		remove((infile+'.'+std::to_string(tid)).c_str());
	}
	f_order.close();
	f_readlength.close();
	//write remaining singleton reads now
	std::ofstream f_singleton(outfile_singleton);
	f_order.open(infile_order,std::ios::binary|std::ofstream::app);
	f_readlength.open(infile_readlength,std::ios::binary|std::ofstream::app);
	std::ofstream f_N(outfile_N);
	uint32_t matched_s = numreads_s;
	for (uint32_t i = 0; i < numreads_s; i++)
		if(remainingreads[i] == 1)
		{
			matched_s--;
			f_order.write((char*)&order_s[i],sizeof(uint32_t));
			f_readlength.write((char*)&read_lengths_s[i],sizeof(uint8_t));
			f_singleton << bitsettostring(read[i],read_lengths_s[i]);
		}
	uint32_t matched_N = numreads_N;
	for (uint32_t i = numreads_s; i < numreads_s+numreads_N; i++)
		if(remainingreads[i] == 1)
		{
			matched_N--;
			f_N << bitsettostring(read[i],read_lengths_s[i]);
			f_order.write((char*)&order_s[i],sizeof(uint32_t));
			f_readlength.write((char*)&read_lengths_s[i],sizeof(uint8_t));
		}
	f_order.close();
	f_readlength.close();
	f_N.close();
	f_singleton.close();
	delete[] remainingreads;
	packbits();
	std::cout << "Encoding done:\n"; 
	std::cout << matched_s << " singleton reads were aligned\n";
	std::cout << matched_N << " reads with N were aligned\n";
	return;
}

std::string buildcontig(std::list<contig_reads> &current_contig, uint32_t list_size)
{
	if(list_size == 1)
		return (current_contig.front()).read;
	auto current_contig_it = current_contig.begin();
	int64_t currentpos = 0, currentsize = 0,to_insert;
	std::vector<std::array<long,4>> count;
	for(; current_contig_it != current_contig.end(); ++current_contig_it)
	{
		if(current_contig_it == current_contig.begin())//first read
			to_insert = (*current_contig_it).read_length;
		else
		{
			currentpos = (*current_contig_it).pos;
			if(currentpos + (*current_contig_it).read_length > currentsize)
				to_insert = currentpos + (*current_contig_it).read_length - currentsize;
			else
				to_insert = 0;
		}
		count.insert(count.end(),to_insert,{0,0,0,0});
		currentsize = currentsize + to_insert;
		for(long i = 0; i < (*current_contig_it).read_length; i++)
			count[currentpos+i][chartolong[(*current_contig_it).read[i]]] += 1;
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

void writecontig(std::string &ref,std::list<contig_reads> &current_contig, std::ofstream& f_seq, std::ofstream& f_pos, std::ofstream& f_noise, std::ofstream& f_noisepos, std::ofstream& f_order, std::ofstream &f_RC, std::ofstream &f_readlength)
{
	f_seq << ref;
	char c;
	long prevj = 0;
	auto current_contig_it = current_contig.begin();
	long currentpos,prevpos = 0;
	for(;current_contig_it!=current_contig.end(); ++current_contig_it)
	{
		currentpos = (*current_contig_it).pos;
		prevj = 0;
		for(long j = 0; j < (*current_contig_it).read_length; j++)
			if((*current_contig_it).read[j] != ref[currentpos+j])
			{
				f_noise<<enc_noise[ref[currentpos+j]][(*current_contig_it).read[j]];
				c = j-prevj;
				f_noisepos<<c;
				prevj = j;
			}
		f_noise << "\n";
		if(current_contig_it == current_contig.begin())
			c = max_readlen;
		else
			c = currentpos - prevpos;
		f_pos << c;
		f_order.write((char*)&((*current_contig_it).order),sizeof(uint32_t));
		f_readlength.write((char*)&((*current_contig_it).read_length),sizeof(uint8_t));
		f_RC << (*current_contig_it).RC;
		prevpos = currentpos;
	}
	return;
}

void generatemasks(bitset mask[MAX_READ_LEN][MAX_READ_LEN])
{
	for(int i = 0; i < max_readlen; i++)
	{
		for(int j = 0; j < max_readlen; j++)	
		{
			mask[i][j].reset();
			for(int k = 3*i; k < 3*max_readlen-3*j; k++)
				mask[i][j][k] = 1;
		}
	}
	return;
}

void packbits()
{
	#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		//seq
		std::ifstream in_seq(outfile_seq+'.'+std::to_string(tid));
		std::ofstream f_seq(outfile_seq+'.'+std::to_string(tid)+".tmp",std::ios::binary);
		std::ofstream f_seq_tail(outfile_seq+'.'+std::to_string(tid)+".tail");
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
		in_seq.open(outfile_seq+'.'+std::to_string(tid));
		char dnabase[8];
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
		remove((outfile_seq+'.'+std::to_string(tid)).c_str());
		rename((outfile_seq+'.'+std::to_string(tid)+".tmp").c_str(),(outfile_seq+'.'+std::to_string(tid)).c_str());		
		
		//rev
		std::ifstream in_rev(infile_RC+'.'+std::to_string(tid));
		std::ofstream f_rev(infile_RC+'.'+std::to_string(tid)+".tmp",std::ios::binary);
		std::ofstream f_rev_tail(infile_RC+'.'+std::to_string(tid)+".tail");
		file_len=0;
		while(in_rev >> std::noskipws >> c)
			file_len++;
		basetoint['d'] = 0;
		basetoint['r'] = 1;
		
		in_rev.close();
		in_rev.open(infile_RC+'.'+std::to_string(tid));
		for(uint64_t i = 0; i < file_len/8; i++)
		{
			in_rev.read(dnabase,8);
			dnabin = 128*basetoint[dnabase[7]]+64*basetoint[dnabase[6]]+
				 32*basetoint[dnabase[5]]+16*basetoint[dnabase[4]]+
				 8*basetoint[dnabase[3]]+4*basetoint[dnabase[2]]+
				 2*basetoint[dnabase[1]]+basetoint[dnabase[0]];
			f_rev.write((char*)&dnabin,sizeof(uint8_t));
		}
		f_rev.close();
		in_rev.read(dnabase,file_len%8);
		for(int i=0; i<file_len%8;i++)
			f_rev_tail << dnabase[i];
		f_rev_tail.close();
		remove((infile_RC+'.'+std::to_string(tid)).c_str());
		rename((infile_RC+'.'+std::to_string(tid)+".tmp").c_str(),(infile_RC+'.'+std::to_string(tid)).c_str());
	}
	//singleton
	std::ifstream in_singleton(outfile_singleton);
	std::ofstream f_singleton(outfile_singleton+".tmp",std::ios::binary);
	std::ofstream f_singleton_tail(outfile_singleton+".tail");
	uint64_t file_len=0;
	char c;
	while(in_singleton >> std::noskipws >> c)
		file_len++;
	uint8_t basetoint[128];
	basetoint['A'] = 0;
	basetoint['C'] = 1;
	basetoint['G'] = 2;
	basetoint['T'] = 3;
	in_singleton.close();
	in_singleton.open(outfile_singleton);
	char dnabase[8];
	uint8_t dnabin;
	for(uint64_t i = 0; i < file_len/4; i++)
	{
		in_singleton.read(dnabase,4);
		
		dnabin = 64*basetoint[dnabase[3]]+16*basetoint[dnabase[2]]+4*
			basetoint[dnabase[1]]+basetoint[dnabase[0]];
		f_singleton.write((char*)&dnabin,sizeof(uint8_t));
	}
	f_singleton.close();
	in_singleton.read(dnabase,file_len%4);
	for(int i=0; i<file_len%4;i++)
		f_singleton_tail << dnabase[i];
	f_singleton_tail.close();
	in_singleton.close();
	remove((outfile_singleton).c_str());
	rename((outfile_singleton+".tmp").c_str(),(outfile_singleton).c_str());		
	return;
}

void setglobalarrays()
{
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
	std::ifstream f_numreads(infilenumreads, std::ios::binary);
	uint32_t numreads_clean, numreads_total;
	f_numreads.read((char*)&numreads_clean,sizeof(uint32_t));
	f_numreads.read((char*)&numreads_total,sizeof(uint32_t));
	f_numreads.close();
	
	std::ifstream myfile_s(infile+".singleton", std::ifstream::in);
	numreads_s = 0;
	std::string line; 
	while (std::getline(myfile_s, line))
		++numreads_s;
	myfile_s.close();
	numreads = numreads_clean - numreads_s;
	numreads_N = numreads_total - numreads_clean;

	std::cout << "Maximum Read length: " << max_readlen << std::endl;
	std::cout << "Number of non-singleton reads: " << numreads << std::endl;
	std::cout << "Number of singleton reads: " << numreads_s << std::endl;
	std::cout << "Number of reads with N: " << numreads_N << std::endl;
}

void correct_order(uint32_t *order_s)
{
	uint32_t numreads_total = numreads+numreads_s+numreads_N;
	bool *read_flag_N = new bool [numreads_total]();
	//bool array indicating N reads
	for(uint32_t i = 0; i < numreads_N; i++)
	{
		read_flag_N[order_s[numreads_s+i]] = true;
	}

	uint32_t *cumulative_N_reads = new uint32_t [numreads + numreads_s];
	//number of reads occuring before pos in clean reads
	uint32_t pos_in_clean = 0, num_N_reads_till_now = 0;
	for(uint32_t i = 0; i < numreads_total; i++)
	{
		if(read_flag_N[i] == true)
			num_N_reads_till_now++;
		else
			cumulative_N_reads[pos_in_clean++] = num_N_reads_till_now;
	}
	
	//First correct the order for singletons
	for(uint32_t i = 0; i < numreads_s; i++)
		order_s[i] += cumulative_N_reads[order_s[i]];

	//Now correct for clean reads (this is stored on file)
	for(int tid = 0; tid < num_thr; tid++)
	{
		std::ifstream fin_order(infile_order+'.'+std::to_string(tid), std::ios::binary);
		std::ofstream fout_order(infile_order+'.'+std::to_string(tid)+".tmp", std::ios::binary);
		uint32_t pos;
		fin_order.read((char*)&pos, sizeof(uint32_t));
		while(!fin_order.eof())
		{
			pos += cumulative_N_reads[pos];
			fout_order.write((char*)&pos, sizeof(uint32_t));
			fin_order.read((char*)&pos, sizeof(uint32_t));
		}
		fin_order.close();
		fout_order.close();
		remove((infile_order+'.'+std::to_string(tid)).c_str());
		rename((infile_order+'.'+std::to_string(tid)+".tmp").c_str(),(infile_order+'.'+std::to_string(tid)).c_str());
	}
	remove(infile_order_N.c_str());
	delete[] read_flag_N;
	delete[] cumulative_N_reads;
	return;
}
		
void stringtobitset(std::string s, uint8_t readlen, bitset &b)
{
	for(int i = 0; i < readlen; i++)
		b |= basemask[i][s[i]];
	return;
}

void readsingletons(bitset *read, uint32_t *order_s, uint8_t *read_lengths_s)
{
	//not parallelized right now since these are very small number of reads
	std::ifstream f(infile+".singleton", std::ifstream::in);
	std::string s;
	for(uint32_t i = 0; i < numreads_s; i++)
	{
		std::getline(f,s);
		read_lengths_s[i] = s.length();
		stringtobitset(s,read_lengths_s[i],read[i]);
	}
	f.close();
	f.open(infile_N);
	for(uint32_t i = numreads_s; i < numreads_s+numreads_N; i++)
	{
		std::getline(f,s);
		read_lengths_s[i] = s.length();
		stringtobitset(s,read_lengths_s[i],read[i]);
	}
	std::ifstream f_order_s(infile_order+".singleton", std::ios::binary);
	for(uint32_t i = 0; i < numreads_s; i++)
		f_order_s.read((char*)&order_s[i],sizeof(uint32_t));
	f_order_s.close();
	std::ifstream f_order_N(infile_order_N, std::ios::binary);
	for(uint32_t i = numreads_s; i < numreads_s+numreads_N; i++)
		f_order_N.read((char*)&order_s[i],sizeof(uint32_t));	
	f_order_N.close();
}
/*
void readsingletons(bitset *read, uint32_t *order_s)
{
	#pragma omp parallel
	{
	int tid = omp_get_thread_num();
	uint32_t i, stop;	
	//doing initial setup and first read
	i = uint64_t(tid)*numreads_s/omp_get_num_threads();//spread out first read equally
	stop = uint64_t(tid+1)*numreads_s/omp_get_num_threads();
	if(tid == omp_get_num_threads()-1)
		stop = numreads_s;
	std::ifstream f(infile+".singleton", std::ifstream::in);
	f.seekg(uint64_t(i)*(readlen+1), f.beg);
	std::string s;
	while(i < stop)
	{
		std::getline(f,s);
		read[i] = stringtobitset(s);
		i++;
	}
	f.close();
	}
	#pragma omp parallel
	{
	int tid = omp_get_thread_num();
	uint32_t i, stop;	
	//doing initial setup and first read
	i = uint64_t(tid)*numreads_N/omp_get_num_threads();//spread out first read equally
	stop = uint64_t(tid+1)*numreads_N/omp_get_num_threads();
	if(tid == omp_get_num_threads()-1)
		stop = numreads_N;
	std::ifstream f(infile_N, std::ifstream::in);
	f.seekg(uint64_t(i)*(readlen+1), f.beg);
	std::string s;
	while(i < stop)
	{
		std::getline(f,s);
		read[numreads_s+i] = stringtobitset(s);
		i++;
	}
	f.close();
	}
	std::ifstream f_order_s(infile_order+".singleton", std::ios::binary);
	for(uint32_t i = 0; i < numreads_s; i++)
		f_order_s.read((char*)&order_s[i],sizeof(uint32_t));
	f_order_s.close();
	std::ifstream f_order_N(infile_order_N, std::ios::binary);
	for(uint32_t i = numreads_s; i < numreads_s+numreads_N; i++)
		f_order_N.read((char*)&order_s[i],sizeof(uint32_t));	
	f_order_N.close();
	return;
}
*/

void generateindexmasks(bitset *mask1)
//masks for dictionary positions
{
	for(int j = 0; j < numdict_s; j++)
		mask1[j].reset();
	for(int j = 0; j < numdict_s; j++)
		for(int i = 3*dict_start[j]; i < 3*(dict_end[j]+1); i++)
			mask1[j][i] = 1;
	return;
}


void constructdictionary(bitset *read, bbhashdict *dict, uint8_t *read_lengths_s)
{
	bitset mask[numdict_s];
	generateindexmasks(mask);
	for(int j = 0; j < numdict_s; j++)
	{
		uint64_t *ull = new uint64_t[numreads_s+numreads_N];
		#pragma omp parallel
		{
		bitset b;
		int tid = omp_get_thread_num();
		uint32_t i, stop;
		i = uint64_t(tid)*(numreads_s+numreads_N)/omp_get_num_threads();
		stop = uint64_t(tid+1)*(numreads_s+numreads_N)/omp_get_num_threads();
		if(tid == omp_get_num_threads()-1)
			stop = numreads_s+numreads_N;
		//compute keys and and store in ull
		for(; i < stop; i++)
		{
			b = read[i]&mask[j];
			ull[i] = (b>>3*dict_start[j]).to_ullong();
		}
		}//parallel end
		
		//remove keys corresponding to reads shorter than dict_end[j]
		dict[j].dict_numreads = 0;
		for(uint32_t i = 0; i < numreads_s+numreads_N; i++)
		{
			if(read_lengths_s[i] > dict_end[j])
			{
				ull[dict[j].dict_numreads] = ull[i];
				dict[j].dict_numreads++;
			}
		}
		
		//write ull to file
		#pragma omp parallel
		{
		int tid = omp_get_thread_num();
		std::ofstream foutkey(basedir+std::string("/keys.bin.")+std::to_string(tid),std::ios::binary);
		uint32_t i, stop;
		i = uint64_t(tid)*dict[j].dict_numreads/omp_get_num_threads();
		stop = uint64_t(tid+1)*dict[j].dict_numreads/omp_get_num_threads();
		if(tid == omp_get_num_threads()-1)
			stop = dict[j].dict_numreads;
		for(; i < stop; i++)
			foutkey.write((char*)&ull[i], sizeof(uint64_t));
		foutkey.close();
		}//parallel end

		//deduplicating ull
		std::sort(ull,ull+dict[j].dict_numreads);
		uint32_t k = 0;
		for (uint32_t i = 1; i < dict[j].dict_numreads; i++) 
		        if (ull[i] != ull[k])         
				ull[++k] = ull[i];
		dict[j].numkeys = k+1;
		//construct mphf
		auto data_iterator = boomphf::range(static_cast<const u_int64_t*>(ull), static_cast<const u_int64_t*>(ull+dict[j].numkeys));
		double gammaFactor = 5.0;//balance between speed and memory
		dict[j].bphf = new boomphf::mphf<u_int64_t,hasher_t>(dict[j].numkeys,data_iterator,num_thr,gammaFactor,true,false);
	
		delete[] ull;

		//compute hashes for all reads
		#pragma omp parallel
		{
		int tid = omp_get_thread_num();	
		std::ifstream finkey(basedir+std::string("/keys.bin.")+std::to_string(tid),std::ios::binary);
		std::ofstream fouthash(basedir+std::string("/hash.bin.")+std::to_string(tid)+'.'+std::to_string(j),std::ios::binary);
		uint64_t currentkey,currenthash;
		uint32_t i, stop;
		i = uint64_t(tid)*dict[j].dict_numreads/omp_get_num_threads();
		stop = uint64_t(tid+1)*dict[j].dict_numreads/omp_get_num_threads();
		if(tid == omp_get_num_threads()-1)
			stop = dict[j].dict_numreads;
		for(; i < stop; i++)
		{
			finkey.read((char*)&currentkey, sizeof(uint64_t));
			currenthash = (dict[j].bphf)->lookup(currentkey);
			fouthash.write((char*)&currenthash, sizeof(uint64_t));
		}
		finkey.close();
		remove((basedir+std::string("/keys.bin.")+std::to_string(tid)).c_str());
		fouthash.close();
		}//parallel end
		
	}
		
	//for rest of the function, use numdict threads to parallelize
	omp_set_num_threads(std::min(numdict_s,num_thr));
	#pragma omp parallel
	{
	#pragma omp for	
	for(int j = 0; j < numdict_s; j++)
	{
		//fill startpos by first storing numbers and then doing cumulative sum
		dict[j].startpos = new uint32_t[dict[j].numkeys+1]();//1 extra to store end pos of last key
		uint64_t currenthash;
		for(int tid = 0; tid < num_thr; tid++)
		{
			std::ifstream finhash(basedir+std::string("/hash.bin.")+std::to_string(tid)+'.'+std::to_string(j),std::ios::binary);
			finhash.read((char*)&currenthash,sizeof(uint64_t));
			while(!finhash.eof())
			{
				dict[j].startpos[currenthash+1]++;
				finhash.read((char*)&currenthash,sizeof(uint64_t));
			}
			finhash.close();
		}
	
		dict[j].empty_bin = new bool[dict[j].numkeys]();
		for(uint32_t i = 1; i < dict[j].numkeys; i++)
			dict[j].startpos[i] =  dict[j].startpos[i] +  dict[j].startpos[i-1];

		//insert elements in the dict array
		dict[j].read_id = new uint32_t[dict[j].dict_numreads];
		uint32_t i = 0;
		for(int tid = 0; tid < num_thr; tid++)
		{
			std::ifstream finhash(basedir+std::string("/hash.bin.")+std::to_string(tid)+'.'+std::to_string(j),std::ios::binary);
			finhash.read((char*)&currenthash,sizeof(uint64_t));
			while(!finhash.eof())
			{
				while(read_lengths_s[i] <= dict_end[j])
					i++;
				dict[j].read_id[dict[j].startpos[currenthash]++] = i;
				i++;
				finhash.read((char*)&currenthash,sizeof(uint64_t));
			}
			finhash.close();
			remove((basedir+std::string("/hash.bin.")+std::to_string(tid)+'.'+std::to_string(j)).c_str());
		}
		
		//correcting startpos array modified during insertion
		for(int64_t i = dict[j].numkeys; i >= 1 ; i--)	
			dict[j].startpos[i] = dict[j].startpos[i-1];
		dict[j].startpos[0] = 0;
	}//for end
	}//parallel end
	omp_set_num_threads(num_thr);
	return;
}

void bbhashdict::findpos(int64_t *dictidx, uint32_t &startposidx)
{
	dictidx[0] = startpos[startposidx];
	auto endidx = startpos[startposidx+1];
	if(read_id[endidx-1] == numreads)//means exactly one read has been removed
		dictidx[1] = endidx-1;
	else if(read_id[endidx-1] == numreads+1)//means two or more reads have been removed (in this case second last entry stores the number of reads left)
		dictidx[1] = dictidx[0] + read_id[endidx-2];
	else
		dictidx[1] = endidx;//no read deleted
	return;
}

void bbhashdict::remove(int64_t *dictidx, uint32_t &startposidx, uint32_t current)
{
	auto size = dictidx[1] - dictidx[0];
	if(size == 1)//just one read left in bin
	{
		empty_bin[startposidx] = 1;
		return; //need to keep one read to check during matching
	}
	uint32_t pos = std::lower_bound(read_id+dictidx[0],read_id+dictidx[1],current)-(read_id+dictidx[0]);
	for(uint32_t i = dictidx[0]+pos; i < dictidx[1]-1; i++)
		read_id[i] = read_id[i+1];
	auto endidx = startpos[startposidx+1];
	if(dictidx[1] == endidx)//this is first read to be deleted
		read_id[endidx-1] = numreads;
	else if(read_id[endidx-1] == numreads)//exactly one read has been deleted till now
	{
		read_id[endidx-1] = numreads + 1;
		read_id[endidx-2] = size - 1;//number of reads left in bin
	}
	else//more than two reads have been deleted
		read_id[endidx-2]--;

	return;	
}

std::string bitsettostring(bitset b, uint8_t readlen)
{
	char s[MAX_READ_LEN+1];
	s[readlen] = '\0';
	unsigned long long ull,rem;
	for(int i = 0; i < 3*readlen/63+1; i++)
	{	
		ull = (b&mask63).to_ullong();
		b>>=63;
		for(int j = 21*i  ; j < 21*i+21 && j < readlen ; j++)
		{
			s[j] = revinttochar[ull%8];	
			ull/=8;
		}
	}
	std::string s1 = s;
	return s1;
}


std::string reverse_complement(std::string s, uint8_t readlen)
{
	std::string s1(s);
	for(int j = 0; j < readlen; j++)
		s1[j] = chartorevchar[s[readlen-j-1]];
	return s1;
}
