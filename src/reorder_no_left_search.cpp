#include <iostream>
#include <fstream>
#include <bitset>
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>
#include <string>
#include "config.h"
#include "BooPHF.h"
#include <omp.h>
#include <atomic>
#include <cstdio>
#include <cstdlib>

uint32_t num_locks = 0x1000000; //limits on number of locks (power of 2 for fast mod)

typedef boomphf::SingleHashFunctor<u_int64_t>  hasher_t;
typedef boomphf::mphf<  u_int64_t, hasher_t  > boophf_t;
typedef std::bitset<2*MAX_READ_LEN> bitset;

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
		delete[] empty_bin;
		delete bphf;
	}	
};

uint32_t numreads = 0;

int maxshift, thresh = 4, numdict = 2, maxsearch = 1000, num_thr, max_readlen;
int dict_start[2];
int dict_end[2]; 

std::string basedir;
std::string infile;
std::string infilenumreads;
std::string outfile;
std::string outfileRC;
std::string outfileflag;
std::string outfilepos;
std::string outfileorder;
std::string outfilereadlength;

//Some global arrays (some initialized in setglobalarrays())
char revinttochar[4] = {'A','G','C','T'};//used in bitsettostring
char inttochar[] = {'A','C','G','T'};
char chartorevchar[128];//A-T etc for reverse complement
int chartoint[128];//A-0,C-1 etc. used in updaterefcount

bitset basemask[MAX_READ_LEN][128];//bitset for A,G,C,T at each position 
//used in stringtobitset, chartobitset and bitsettostring
bitset positionmask[MAX_READ_LEN];//bitset for each position (1 at two bits and 0 elsewhere)
//used in bitsettostring
bitset mask64;//bitset with 64 bits set to 1 (used in bitsettostring for conversion to ullong)

void stringtobitset(std::string &s, uint8_t readlen, bitset &b);

void bitsettostring(bitset b,char *s, uint8_t readlen);

void readDnaFile(bitset *read, uint8_t *read_lengths);

void constructdictionary(bitset *read, bbhashdict *dict, uint8_t *read_lengths);

//void constructdictionary(bitset *read, spp::sparse_hash_map<uint64_t,uint32_t*> *dict);

void generatemasks(bitset mask[MAX_READ_LEN][MAX_READ_LEN]);
//mask for zeroing the end bits (needed while reordering to compute Hamming distance between shifted reads)

void reorder(bitset *read, bbhashdict *dict, uint8_t *read_lengths);

bool search_match(bitset &ref, bitset *mask1, omp_lock_t *dict_lock, omp_lock_t *read_lock, bitset mask[MAX_READ_LEN][MAX_READ_LEN], uint8_t *read_lengths, bool *remainingreads, bitset *read, bbhashdict *dict, uint32_t &k, bool rev, int shift, int &ref_len);

void writetofile(bitset *read, uint8_t *read_lengths);

void updaterefcount(bitset &cur, bitset &ref, bitset &revref, int count[][MAX_READ_LEN], bool resetcount, bool rev, int shift, uint8_t cur_readlen, int &ref_len);
//update the clean reference and count matrix

void reverse_complement(char *s, char *s1, uint8_t readlen);

void chartobitset(char *s, uint8_t readlen, bitset &b);

void setglobalarrays();
//setting arrays like chartoint etc.

int main(int argc, char** argv)
{
	basedir = std::string(argv[1]);
	infile = basedir + "/input_clean.dna";
	outfile = basedir + "/temp.dna";
	outfileRC = basedir + "/read_rev.txt";
	outfileflag = basedir + "/tempflag.txt";
	outfilepos = basedir + "/temppos.txt";
	outfileorder = basedir + "/read_order.bin";
	outfilereadlength = basedir + "/read_lengths.bin";
	infilenumreads = basedir + "/numreads.bin";
	
	max_readlen = atoi(argv[2]);
	num_thr = atoi(argv[3]);
	maxshift = max_readlen/2;
	dict_start[0] =  max_readlen > 100 ? max_readlen/2-32 : max_readlen/2-max_readlen*32/100;
	dict_end[0] = max_readlen/2-1;
	dict_start[1] = max_readlen/2;
	dict_end[1] = max_readlen > 100 ? max_readlen/2-1+32 : max_readlen/2-1+max_readlen*32/100;
	
	std::ifstream f_numreads(infilenumreads, std::ios::binary);
	f_numreads.read((char*)&numreads,sizeof(uint32_t));	
	omp_set_num_threads(num_thr);	
	setglobalarrays();
	bitset *read = new bitset [numreads];
	uint8_t *read_lengths = new uint8_t[numreads];
	std::cout << "Reading file: " << infile << std::endl;
	readDnaFile(read,read_lengths);
	bbhashdict *dict = new bbhashdict [numdict];
	std::cout << "Constructing dictionaries\n";
	constructdictionary(read,dict,read_lengths);
	std::cout << "Reordering reads\n";
	reorder(read,dict,read_lengths);
	std::cout << "Writing to file\n";
	writetofile(read,read_lengths);	
	delete[] read;
	delete[] dict;
	std::cout << "Done!\n";
	return 0;
}

void setglobalarrays()
{
	chartorevchar['A'] = 'T';
	chartorevchar['C'] = 'G';
	chartorevchar['G'] = 'C';
	chartorevchar['T'] = 'A';
	chartoint['A'] = 0;
	chartoint['C'] = 1;
	chartoint['G'] = 2;
	chartoint['T'] = 3;

	for(int i = 0; i < 64; i++)
		mask64[i] = 1;
	for(int i = 0; i < max_readlen; i++)
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
	

void stringtobitset(std::string &s, uint8_t readlen, bitset &b)
{
	for(int i = 0; i < readlen; i++)
		b |= basemask[i][s[i]];
	return;
}



//TODO: Check parallelization of this function
void readDnaFile(bitset *read, uint8_t *read_lengths)
{
	#pragma omp parallel
	{
	int tid = omp_get_thread_num();
	std::ifstream f(infile, std::ifstream::in);
	std::string s;
	uint32_t i = 0;
	while(std::getline(f,s))
	{	
		if(i%num_thr == tid)
		{
			read_lengths[i] = s.length();
			stringtobitset(s,read_lengths[i],read[i]);			
			i++;
		}	
		else
		{
			i++;
			continue;
		}
	}	
	f.close();
	}
	remove(infile.c_str());
	return;
}

void generateindexmasks(bitset *mask1)
//masks for dictionary positions
{
	for(int j = 0; j < numdict; j++)
		mask1[j].reset();
	for(int j = 0; j < numdict; j++)
		for(int i = 2*dict_start[j]; i < 2*(dict_end[j]+1); i++)
			mask1[j][i] = 1;
	return;
}


void constructdictionary(bitset *read, bbhashdict *dict, uint8_t *read_lengths)
{
	bitset mask[numdict];
	generateindexmasks(mask);
	for(int j = 0; j < numdict; j++)
	{
		uint64_t *ull = new uint64_t[numreads];
		#pragma omp parallel
		{
		bitset b;
		int tid = omp_get_thread_num();
		uint32_t i, stop;
		i = uint64_t(tid)*numreads/omp_get_num_threads();
		stop = uint64_t(tid+1)*numreads/omp_get_num_threads();
		if(tid == omp_get_num_threads()-1)
			stop = numreads;
		//compute keys and and store in ull
		for(; i < stop; i++)
		{
			b = read[i]&mask[j];
			ull[i] = (b>>2*dict_start[j]).to_ullong();
		}
		}//parallel end
		
		//remove keys corresponding to reads shorter than dict_end[j]
		dict[j].dict_numreads = 0;
		for(uint32_t i = 0; i < numreads; i++)
		{
			if(read_lengths[i] > dict_end[j])
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
	omp_set_num_threads(std::min(numdict,num_thr));
	#pragma omp parallel
	{
	#pragma omp for	
	for(int j = 0; j < numdict; j++)
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
				while(read_lengths[i] <= dict_end[j])
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
//	std::move(read_id+dictidx[0]+pos+1,read_id+dictidx[1],read_id+dictidx[0]+pos);
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

void reorder(bitset *read, bbhashdict *dict, uint8_t *read_lengths)
{
	omp_lock_t *dict_lock = new omp_lock_t [num_locks];
	omp_lock_t *read_lock = new omp_lock_t [num_locks];
	for(int j = 0; j < num_locks; j++)
	{
		omp_init_lock(&dict_lock[j]);
		omp_init_lock(&read_lock[j]);
	}
	bitset mask[MAX_READ_LEN][MAX_READ_LEN];
	generatemasks(mask);
	bitset mask1[numdict];
	generateindexmasks(mask1);
	bool *remainingreads = new bool [numreads];
	std::fill(remainingreads, remainingreads+numreads,1);

	//we go through remainingreads array from behind as that speeds up deletion from bin arrays
	
	uint32_t firstread = 0, unmatched[num_thr];
	#pragma omp parallel
	{
	int tid = omp_get_thread_num();
	std::string tid_str = std::to_string(tid);	
	std::ofstream foutRC(outfileRC + '.' + tid_str,std::ofstream::out);
	std::ofstream foutflag(outfileflag + '.' + tid_str,std::ofstream::out);
	std::ofstream foutpos(outfilepos + '.' + tid_str,std::ofstream::out|std::ios::binary);
	std::ofstream foutorder(outfileorder + '.' + tid_str,std::ofstream::out|std::ios::binary);
	std::ofstream foutorder_s(outfileorder + ".singleton." + tid_str,std::ofstream::out|std::ios::binary);
	std::ofstream foutlength(outfilereadlength+'.'+tid_str,std::ofstream::out|std::ios::binary);
		
	unmatched[tid] = 0;
	bitset ref,revref, b;

	uint32_t first_rid;
	//first_rid represents first read of contig, used for left searching
	int count[4][MAX_READ_LEN];
	int64_t dictidx[2];//to store the start and end index (end not inclusive) in the dict read_id array
	uint32_t startposidx;//index in startpos
	bool flag = 0, done = 0, prev_unmatched = false, left_search_start = false, left_search = false;
	//left_search_start - true when left search is starting (to avoid double deletion of first read)
	//left_search - true during left search - needed so that we know when to pick arbitrary read
	uint32_t current, prev;
	uint64_t ull;
	int ref_len;
	int64_t ref_pos, cur_read_pos;
	//starting positions of the ref and the current read in the contig, can be negative during left search or due to RC
	//useful for sorting according to starting position in the encoding stage.
	
	//flag to check if match was found or not
	
	int64_t remainingpos = numreads-1;//used for searching next unmatched read when no match is found
	#pragma omp critical
	{//doing initial setup and first read
		current = firstread;	
		firstread += numreads/omp_get_num_threads();//spread out first read equally
		remainingreads[current] = 0;
		unmatched[tid]++;
	}
	#pragma omp barrier
	updaterefcount(read[current],ref,revref,count,true,false,0,read_lengths[current],ref_len);
	cur_read_pos = 0;
	ref_pos = 0;
	first_rid = current;
	prev_unmatched = true;	
	prev = current;
	while(!done)
	{
		//delete the read from the corresponding dictionary bins (unless we are starting left search)
		if(!left_search_start)
		{
			for(int l = 0; l < numdict; l++)
			{
				if(read_lengths[current] <= dict_end[l])
					continue;
				b = read[current]&mask1[l];
				ull = (b>>2*dict_start[l]).to_ullong();
				startposidx = dict[l].bphf->lookup(ull);
				//check if any other thread is modifying same dictpos
				omp_set_lock(&dict_lock[startposidx & 0xFFFFFF]);
				dict[l].findpos(dictidx,startposidx);
				dict[l].remove(dictidx,startposidx,current);
				omp_unset_lock(&dict_lock[startposidx & 0xFFFFFF]);
			}
		}
		else
		{
			left_search_start = false;
		}
		flag = 0;
		uint32_t k;
		for(int shift = 0; shift < maxshift; shift++)
		{
			//find forward match
			flag = search_match(ref,mask1,dict_lock,read_lock,mask,read_lengths,remainingreads,read,dict,k,false, shift,ref_len);
			if(flag == 1)
			{
				current = k;
				int ref_len_old = ref_len;
				updaterefcount(read[current],ref,revref,count,false,false,shift,read_lengths[current],ref_len);
				if(!left_search)
				{
					cur_read_pos = ref_pos + shift;
					ref_pos = cur_read_pos;
				}	
				else
				{
					cur_read_pos = ref_pos + ref_len_old - shift - read_lengths[current];
					ref_pos = ref_pos + ref_len_old - shift - ref_len;
					
				}
				if(prev_unmatched == true)//prev read not singleton, write it now
				{
					foutRC << 'd';
					foutorder.write((char*)&prev,sizeof(uint32_t));
					foutflag << 0;//for unmatched
					int64_t zero = 0;
					foutpos.write((char*)&zero,sizeof(int64_t));
					foutlength.write((char*)&read_lengths[prev],sizeof(uint8_t));
				}	
				foutRC << (left_search?'r':'d');
				foutorder.write((char*)&current,sizeof(uint32_t));
				foutflag << 1;//for matched
				foutpos.write((char*)&cur_read_pos,sizeof(int64_t));
				foutlength.write((char*)&read_lengths[current],sizeof(uint8_t));
				
				prev_unmatched = false;
				break;
			}

			//find reverse match
			flag = search_match(revref,mask1,dict_lock,read_lock,mask,read_lengths,remainingreads,read,dict,k,true,shift,ref_len);
			if(flag == 1)
			{
				current = k;
				int ref_len_old = ref_len;
				updaterefcount(read[current],ref,revref,count,false,true,shift,read_lengths[current],ref_len);
				if(!left_search)
				{
					cur_read_pos = ref_pos + ref_len_old + shift - read_lengths[current];
					ref_pos = ref_pos + ref_len_old + shift - ref_len;
				}
				else
				{
					cur_read_pos = ref_pos - shift;
					ref_pos = cur_read_pos;
				}
				if(prev_unmatched == true)//prev read not singleton, write it now
				{
					foutRC << 'd';
					foutorder.write((char*)&prev,sizeof(uint32_t));
					foutflag << 0;//for unmatched
					int64_t zero = 0;
					foutpos.write((char*)&zero,sizeof(int64_t));
					foutlength.write((char*)&read_lengths[prev],sizeof(uint8_t));
				}
				foutRC << (left_search?'d':'r');
				foutorder.write((char*)&current,sizeof(uint32_t));
				foutflag << 1;//for matched
				foutpos.write((char*)&cur_read_pos,sizeof(int64_t));
				foutlength.write((char*)&read_lengths[current],sizeof(uint8_t));

				prev_unmatched = false;
				break;
			}
		
			revref<<=2;
			ref>>=2;
		}
		if(flag == 0)//no match found
		{
/*
			if(!left_search) //start left search
			{
				left_search = true;
				left_search_start = true;
				//update ref and count with RC of first_read
				updaterefcount(read[first_rid],ref,revref,count,true,true,0,read_lengths[first_rid],ref_len);
				ref_pos = 0;
				cur_read_pos = 0;
			}
			else //left search done, now pick arbitrary read and start new contig 
*/			{
				left_search = false;
				for(int64_t j = remainingpos; j>=0; j--)
				{
					if(remainingreads[j] == 1)
					{
						omp_set_lock(&read_lock[j & 0xffffff]);
						if(remainingreads[j])//checking again inside critical block
						{
							current = j;
							remainingpos = j-1;
							remainingreads[j] = 0;
							flag = 1;
							unmatched[tid]++;
						}
						omp_unset_lock(&read_lock[j & 0xffffff]);
						if(flag == 1)
							break;
					}
				}		
				if(flag == 0)
				{
					if(prev_unmatched == true)//last read was singleton, write it now
					{
						foutorder_s.write((char*)&prev,sizeof(uint32_t));
					}
					done = 1;//no reads left
				}
				else
				{
					updaterefcount(read[current],ref,revref,count,true,false,0,read_lengths[current],ref_len);
					ref_pos = 0;
					cur_read_pos = 0;
					if(prev_unmatched == true)//prev read singleton, write it now
					{
						foutorder_s.write((char*)&prev,sizeof(uint32_t));
					}
					prev_unmatched = true;
					first_rid = current;
					prev = current;
				}
			}
			
		}
	}//while(!done) end
	
	foutRC.close();
	foutorder.close();
	foutflag.close();
	foutpos.close();
	foutorder_s.close();
	foutlength.close();
//	std::cout << tid << ":Done"<<"\n";
	}//parallel end
	
	delete[] remainingreads;
	std::cout << "Reordering done, "<< std::accumulate(unmatched,unmatched+num_thr,0) <<" were unmatched\n";
	return;
}


void generatemasks(bitset mask[MAX_READ_LEN][MAX_READ_LEN])
{
	for(int i = 0; i < max_readlen; i++)
	{
		for(int j = 0; j < max_readlen; j++)	
		{
			mask[i][j].reset();
			for(int k = 2*i; k < 2*max_readlen-2*j; k++)
				mask[i][j][k] = 1;
		}
	}
	return;
}


bool search_match(bitset &ref, bitset *mask1, omp_lock_t *dict_lock, omp_lock_t *read_lock, bitset mask[MAX_READ_LEN][MAX_READ_LEN], uint8_t *read_lengths, bool *remainingreads, bitset *read, bbhashdict *dict, uint32_t &k, bool rev, int shift, int &ref_len)
//right_aligned means that ref is right aligned. In this case reads need to be shifted before hamming.
//right_alignment is useful for rev and for left search to make sure that the starting positions of reads is monotone
{
	bitset b;
	uint64_t ull;
	int64_t dictidx[2];//to store the start and end index (end not inclusive) in the dict read_id array
	uint32_t startposidx;//index in startpos
	bool flag = 0;
	for(int l = 0; l < numdict; l++)
	{
		//check if dictionary index is within bound
		if(!rev)
		{
			if(dict_end[l]+shift >= ref_len)
				continue;
		}
		else
		{
			if(dict_end[l] >= ref_len+shift || dict_start[l] <= shift)
				continue;
		}
		b = ref&mask1[l];
		ull = (b>>2*dict_start[l]).to_ullong();
		startposidx = dict[l].bphf->lookup(ull);
		if(startposidx >= dict[l].numkeys)//not found
			continue;
		//check if any other thread is modifying same dictpos
		omp_set_lock(&dict_lock[startposidx & 0xFFFFFF]);
		dict[l].findpos(dictidx,startposidx);
		if(dict[l].empty_bin[startposidx])//bin is empty
		{
			omp_unset_lock(&dict_lock[startposidx & 0xFFFFFF]);
			continue;
		}
		uint64_t ull1 = ((read[dict[l].read_id[dictidx[0]]]&mask1[l])>>2*dict_start[l]).to_ullong();
		if(ull == ull1)//checking if ull is actually the key for this bin
		{	
			for (int64_t i = dictidx[1] - 1 ; i >= dictidx[0] && i >= dictidx[1] - maxsearch; i--)
			{
				auto rid = dict[l].read_id[i];
				int hamming;
				if(!rev)
					hamming = ((ref^read[rid])&mask[0][max_readlen-std::min<int>(ref_len-shift,read_lengths[rid])]).count();
				else
					hamming = ((ref^read[rid])&mask[shift][max_readlen-std::min<int>(ref_len+shift,read_lengths[rid])]).count();	
				if(hamming<=thresh)
				{	
					omp_set_lock(&read_lock[rid & 0xFFFFFF]);
					if(remainingreads[rid])
					{
						remainingreads[rid]=0;
						k = rid;
						flag = 1;
					}
					omp_unset_lock(&read_lock[rid & 0xFFFFFF]);
					if(flag == 1)
						break;
				}
			}
		}
		omp_unset_lock(&dict_lock[startposidx & 0xFFFFFF]);
		if(flag == 1)
			break;
	}
	return flag;
}

void writetofile(bitset *read, uint8_t *read_lengths)
{

	//convert bitset to string for all num_thr files in parallel
	#pragma omp parallel
	{
	int tid = omp_get_thread_num();
	std::string tid_str = std::to_string(tid);
	std::ofstream fout(outfile + '.' + tid_str,std::ofstream::out);
	std::ofstream fout_s(outfile + ".singleton." + tid_str,std::ofstream::out);
	std::ifstream finRC(outfileRC + '.' + tid_str,std::ifstream::in);
	std::ifstream finorder(outfileorder + '.' + tid_str,std::ifstream::in|std::ios::binary);
	std::ifstream finorder_s(outfileorder + ".singleton." + tid_str,std::ifstream::in|std::ios::binary);
	char s[MAX_READ_LEN+1],s1[MAX_READ_LEN+1];
	uint32_t current;
	char c;
	while(finRC >> std::noskipws >> c)//read character by character
	{
		finorder.read((char*)&current, sizeof (uint32_t));
		bitsettostring(read[current],s,read_lengths[current]);
		if(c == 'd')
		{
			fout<<s<<"\n";
		}
		else
		{
			reverse_complement(s,s1,read_lengths[current]);
			fout<<s1<<"\n";
		}
	}
	finorder_s.read((char*)&current,sizeof(uint32_t));
	while(!finorder_s.eof())
	{
		bitsettostring(read[current],s,read_lengths[current]);
		fout_s << s << "\n";
		finorder_s.read((char*)&current,sizeof(uint32_t));
	}
	
	fout.close();
	fout_s.close();
	finRC.close();
	finorder.close();
	finorder_s.close();	

	}

	//Now combine the num_thr order files
	std::ofstream fout_s(outfile+".singleton",std::ofstream::out);
	std::ofstream foutorder_s(outfileorder+".singleton",std::ofstream::out|std::ios::binary);
	for(int tid = 0; tid < num_thr; tid++)
	{
		std::string tid_str = std::to_string(tid);
		std::ifstream fin_s(outfile + ".singleton." + tid_str,std::ifstream::in);
		std::ifstream finorder_s(outfileorder + ".singleton." + tid_str,std::ifstream::in|std::ios::binary);
			
		fout_s << fin_s.rdbuf();//write entire file
		foutorder_s << finorder_s.rdbuf();

		//clear error flags which occur on rdbuffing empty file
		fout_s.clear();
		foutorder_s.clear();

		fin_s.close();
		finorder_s.close();
	
		remove((outfile + ".singleton." + tid_str).c_str());
		remove((outfileorder + ".singleton." + tid_str).c_str());
	}
	fout_s.close();
	foutorder_s.close();
	return;
}

void bitsettostring(bitset b, char *s, uint8_t readlen)
{
	unsigned long long ull;
	for(int i = 0; i < 2*readlen/64+1; i++)
	{	
		ull = (b&mask64).to_ullong();
		b>>=64;
		for(int j = 32*i  ; j < 32*i+32 && j < readlen ; j++)
		{
			s[j] = revinttochar[ull%4];	
			ull/=4;
		}
	}
	s[readlen] = '\0';
	return;
}



void reverse_complement(char* s, char* s1, uint8_t readlen)
{
	for(int j = 0; j < readlen; j++)
		s1[j] = chartorevchar[s[readlen-j-1]];
	s1[readlen] = '\0';
	return;
}

void chartobitset(char *s, uint8_t readlen, bitset &b)
{
	b.reset();
	for(int i = 0; i < readlen; i++)
		b |= basemask[i][s[i]];
	return;
}

void updaterefcount(bitset &cur, bitset &ref, bitset &revref, int count[][MAX_READ_LEN], bool resetcount, bool rev, int shift, uint8_t cur_readlen, int &ref_len)
//for var length, shift represents shift of start positions, if read length is small, may not need to shift actually
{
	char s[MAX_READ_LEN+1],s1[MAX_READ_LEN+1],*current;
	bitsettostring(cur,s,cur_readlen);
	if(rev == false)
		current = s;
	else
	{
		reverse_complement(s,s1,cur_readlen);
		current = s1;
	}

	if(resetcount == true)//resetcount - unmatched read so start over
	{
		std::memset(count, 0, sizeof(count[0][0]) * 4 * MAX_READ_LEN);	
		for(int i = 0; i < cur_readlen; i++)
		{	
			count[chartoint[current[i]]][i] = 1;
		}
		ref_len = cur_readlen;
	}
	else
	{
		if(!rev)
		{
			//shift count
			for(int i = 0; i < ref_len-shift; i++)
			{	
				for(int j = 0; j < 4; j++)
					count[j][i] = count[j][i+shift];
				if(i < cur_readlen)
					count[chartoint[current[i]]][i] += 1;
			}
			
			//for the new positions set count to 1
			for(int i = ref_len-shift; i < cur_readlen; i++)
			{	
				for(int j = 0; j < 4; j++)
					count[j][i] = 0;
				count[chartoint[current[i]]][i] = 1;
			}
			ref_len = std::max<int>(ref_len-shift,cur_readlen);
		}
		else//reverse case is quite complicated in the variable length case
		{
			if(cur_readlen-shift >= ref_len)
			{
				for(int i = cur_readlen-shift-ref_len; i < cur_readlen-shift; i++)
				{
					for(int j = 0; j < 4; j++)
						count[j][i] = count[j][i-(cur_readlen-shift-ref_len)];
					count[chartoint[current[i]]][i] += 1;
				}
				for(int i = 0; i < cur_readlen-shift-ref_len; i++)
				{	
					for(int j = 0; j < 4; j++)
						count[j][i] = 0;
					count[chartoint[current[i]]][i] = 1;
				}
				for(int i = cur_readlen-shift; i < cur_readlen; i++)
				{	
					for(int j = 0; j < 4; j++)
						count[j][i] = 0;
					count[chartoint[current[i]]][i] = 1;
				}
				ref_len = cur_readlen;
			}
			else if(ref_len + shift <= max_readlen)
			{
				for(int i = ref_len-cur_readlen+shift; i < ref_len; i++)		
					count[chartoint[current[i-(ref_len-cur_readlen+shift)]]][i] += 1;
				for(int i = ref_len; i < ref_len+shift; i++)
				{
					for(int j = 0; j < 4; j++)
						count[j][i] = 0;
					count[chartoint[current[i-(ref_len-cur_readlen+shift)]]][i] = 1;
				}
				ref_len = ref_len+shift;
			}
			else
			{
				for(int i = 0; i < max_readlen - shift; i++)
				{
					for(int j = 0; j < 4; j++)
						count[j][i] = count[j][i+(ref_len+shift-max_readlen)];
				}
				for(int i = max_readlen - cur_readlen; i < max_readlen - shift; i++)
					count[chartoint[current[i-(max_readlen-cur_readlen)]]][i] += 1;
				for(int i = max_readlen - shift; i < max_readlen; i++)
				{
					for(int j = 0; j < 4; j++)
						count[j][i] = 0;
					count[chartoint[current[i-(max_readlen-cur_readlen)]]][i] = 1;
				}
				ref_len = max_readlen;	
			}
		}	
		//find max of each position to get ref

		for(int i = 0; i < ref_len; i++)
		{
			int max = 0,indmax = 0;
			for(int j = 0; j < 4; j++)
				if(count[j][i]>max)
				{
					max = count[j][i];
					indmax = j;
				}
			current[i] = inttochar[indmax];
		}
	}
	chartobitset(current,ref_len,ref);
	char revcurrent[MAX_READ_LEN];
	reverse_complement(current,revcurrent,ref_len);
	chartobitset(revcurrent,ref_len,revref);
	
	return;
}
