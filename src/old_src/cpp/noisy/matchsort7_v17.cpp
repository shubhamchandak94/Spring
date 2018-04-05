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

//#define readlen 100
//#define num_thr 8

typedef boomphf::SingleHashFunctor<u_int64_t>  hasher_t;
typedef boomphf::mphf<  u_int64_t, hasher_t  > boophf_t;

class bbhashdict
{
	public:
	boophf_t * bphf;
	uint32_t numkeys;
	uint32_t *startpos;
	uint32_t *read_id;
	bool *empty_bin;
	void findpos(uint32_t *dictidx, uint32_t &startposidx);
	void remove(uint32_t *dictidx, uint32_t &startposidx, uint32_t current);
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

uint32_t numreads = 0;

std::string infile;
std::string outfile;
std::string outfileRC;
std::string outfileflag;
std::string outfilepos;
std::string outfileorder;
std::string outdir;

//Some global arrays (some initialized in setglobalarrays())
char revinttochar[4] = {'A','G','C','T'};//used in bitsettostring
char inttochar[] = {'A','C','G','T'};
char chartorevchar[128];//A-T etc for reverse complement
int chartoint[128];//A-0,C-1 etc. used in updaterefcount
int *dict_start;
int *dict_end; 
uint8_t *hamming32; 
std::bitset<2*readlen> basemask[readlen][128];//bitset for A,G,C,T at each position 
//used in stringtobitset, chartobitset and bitsettostring
std::bitset<2*readlen> positionmask[readlen];//bitset for each position (1 at two bits and 0 elsewhere)
//used in bitsettostring
std::bitset<2*readlen> mask64;//bitset with 64 bits set to 1 (used in bitsettostring for conversion to ullong)
std::bitset<2*readlen> first32,last32;

std::bitset<2*readlen> stringtobitset(std::string s);

void bitsettostring(std::bitset<2*readlen> b,char *s);

void readDnaFile(std::bitset<2*readlen> *read, uint32_t *read1, uint32_t* read2);

void getDataParams();

void constructdictionary(std::bitset<2*readlen> *read, bbhashdict *dict);

//void constructdictionary(std::bitset<2*readlen> *read, spp::sparse_hash_map<uint64_t,uint32_t*> *dict);

void generatemasks(std::bitset<2*readlen> *mask,std::bitset<2*readlen> *revmask);
//mask for zeroing the end bits (needed while reordering to compute Hamming distance between shifted reads)

void reorder(std::bitset<2*readlen> *read, uint32_t *read1, uint32_t* read2, bbhashdict *dict);

void writetofile(std::bitset<2*readlen> *read);

void updaterefcount(std::bitset<2*readlen> cur, std::bitset<2*readlen> &ref, std::bitset<2*readlen> &revref, int count[][readlen], bool resetcount, bool rev, int shift);
//update the clean reference and count matrix

void reverse_complement(char *s, char *s1);

std::bitset<2*readlen> chartobitset(char &s);

void setglobalarrays();
//setting arrays like chartoint etc.

int main(int argc, char** argv)
{
	std::string basedir = std::string(argv[1]);
	outdir = basedir + "/output/";
	infile = basedir + "/input_clean.dna";
	outfile = basedir + "/output/temp.dna";
	outfileRC = basedir + "/output/read_rev.txt";
	outfileflag = basedir + "/output/tempflag.txt";
	outfilepos = basedir + "/output/temppos.txt";
	outfileorder = basedir + "/output/read_order.bin";
	getDataParams(); //populate numreads, readlen
	omp_set_num_threads(num_thr);	
	setglobalarrays();
	std::bitset<2*readlen> *read = new std::bitset<2*readlen> [numreads];
	uint32_t *read1 = new uint32_t [numreads];
	uint32_t *read2 = new uint32_t [numreads];
	std::cout << "Reading file: " << infile << std::endl;
	readDnaFile(read, read1, read2);
	std::cout << read[0] <<"\n";
	std::cout << read1[0] <<"\n";
	std::cout << read2[0] <<"\n";
	std::cout << hamming32[0] << "\t" << hamming32[65535] <<"\t" <<  hamming32[10000] << "\n";
	std::cout << first32 << "\n";
	std::cout << last32 << "\n";
	bbhashdict dict[numdict];
	std::cout << "Constructing dictionaries\n";
	constructdictionary(read,dict);
	std::cout << "Reordering reads\n";
	reorder(read,read1, read2, dict);
	std::cout << "Writing to file\n";
	writetofile(read);	
	delete[] read,read1,read2;
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
	#if numdict == 1
	{
		dict_start = new int[1];
		dict_end = new int[1];
		dict_start[0] = dict1_start;
		dict_end[0] = dict1_end;
	}
	#elif numdict == 2
	{
		dict_start = new int[2];
		dict_end = new int[2];
		dict_start[0] = dict1_start;
		dict_end[0] = dict1_end;
		dict_start[1] = dict2_start;
		dict_end[1] = dict2_end;
	}
	#elif numdict == 3
	{
		dict_start = new int[3];
		dict_end = new int[3];
		dict_start[0] = dict1_start;
		dict_end[0] = dict1_end;
		dict_start[1] = dict2_start;
		dict_end[1] = dict2_end;
		dict_start[2] = dict3_start;
		dict_end[2] = dict3_end;
	}
	#elif numdict == 4
	{
		dict_start = new int[4];
		dict_end = new int[4];
		dict_start[0] = dict1_start;
		dict_end[0] = dict1_end;
		dict_start[1] = dict2_start;
		dict_end[1] = dict2_end;
		dict_start[2] = dict3_start;
		dict_end[2] = dict3_end;
		dict_start[3] = dict4_start;
		dict_end[3] = dict4_end;
	}
	#endif
	hamming32 = new uint8_t[4294967296];
	if(hamming32 == NULL) 
	{
		std::cout << "asfdasfasdfasdfafds\n\n\n";
		return;
	}
	for(uint64_t i = 0; i < 4294967296; i++)
	{
		std::bitset<32> b(i);
		hamming32[i] = b.count();
	}
	for(int i = 0; i < 32; i++)
		first32[i] = 1;
	for(int i = 2*readlen - 1; i >= 2*readlen -32; i--)
		last32[i] = 1; 	
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
	

std::bitset<2*readlen> stringtobitset(std::string s)
{
	std::bitset<2*readlen> b;
	for(int i = 0; i < readlen; i++)
		b |= basemask[i][s[i]];
	return b;
}

void getDataParams()
{
	uint32_t number_of_lines = 0;
	std::string line;
	std::ifstream myfile(infile, std::ifstream::in);

	std::getline(myfile, line);
	int read_length = line.length();
	myfile.close();
	myfile.open(infile);

	while (std::getline(myfile, line))
	{
	++number_of_lines;
	int _len = line.length(); 
	if( _len != read_length) 
	{
		std::cout << "Read lengths are not the same!: " << read_length << " , " << _len << std::endl;
	}
	}
	//readlen = read_length;
	numreads = number_of_lines;
	std::cout << "Read length: " << read_length << std::endl;
	std::cout << "Number of reads: " << number_of_lines << std::endl;
	myfile.close();
}

void readDnaFile(std::bitset<2*readlen> *read, uint32_t *read1, uint32_t* read2)
{
	#pragma omp parallel
	{
	int tid = omp_get_thread_num();
	uint32_t i, stop;	
	//doing initial setup and first read
	i = uint64_t(tid)*numreads/omp_get_num_threads();//spread out first read equally
	stop = uint64_t(tid+1)*numreads/omp_get_num_threads();
	if(tid == omp_get_num_threads()-1)
		stop = numreads;
	std::ifstream f(infile, std::ifstream::in);
	f.seekg(uint64_t(i)*(readlen+1), f.beg);
	std::string s;
	while(i < stop)
	{
		std::getline(f,s);
		read[i] = stringtobitset(s);
		read1[i] = (read[i]&first32).to_ullong();
		read2[i] = ((read[i]&last32)>>(2*readlen-32)).to_ullong();
		i++;
	}
	f.close();
	}
	return;
}
	
void generateindexmasks(std::bitset<2*readlen> *mask1)
//masks for dictionary positions
{
	for(int j = 0; j < numdict; j++)
		mask1[j].reset();
	for(int j = 0; j < numdict; j++)
		for(int i = 2*dict_start[j]; i < 2*(dict_end[j]+1); i++)
			mask1[j][i] = 1;
	return;
}


void constructdictionary(std::bitset<2*readlen> *read, bbhashdict *dict)
{
	std::bitset<2*readlen> mask[numdict];
	generateindexmasks(mask);
	omp_set_num_threads(std::min(numdict,num_thr));
	//Parallelizing construction of the multiple dictionaries
	#pragma omp parallel 
	{
	#pragma omp for
	for(int j = 0; j < numdict; j++)
	{ 	
		std::bitset<2*readlen> b;
		uint64_t *ull = new uint64_t[numreads];
		std::ofstream foutkey(outdir+std::string("keys.bin.")+std::to_string(j),std::ofstream::out|std::ios::binary);
		//compute keys and write to file and store in ull
		for(uint32_t i = 0; i < numreads; i++)
		{
			b = read[i]&mask[j];
			ull[i] = (b>>2*dict_start[j]).to_ullong();
			foutkey.write((char*)&ull[i], sizeof(uint64_t));
		}
		foutkey.close();
		//deduplicating ull
		std::sort(ull,ull+numreads);
		uint32_t k = 0;
		for (uint32_t i = 1; i < numreads; i++) 
		        if (ull[i] != ull[k])         
				ull[++k] = ull[i];
		dict[j].numkeys = k+1;
		//construct mphf
		auto data_iterator = boomphf::range(static_cast<const u_int64_t*>(ull), static_cast<const u_int64_t*>(ull+dict[j].numkeys));
		double gammaFactor = 5.0;//balance between speed and memory
		dict[j].bphf = new boomphf::mphf<u_int64_t,hasher_t>(dict[j].numkeys,data_iterator,std::max(1,num_thr/numdict),gammaFactor,true,false);
	
		delete[] ull;

		//fill startpos by first storing numbers and then doing cumulative sum
		dict[j].startpos = new uint32_t[dict[j].numkeys+1];//1 extra to store end pos of last key
		std::fill(dict[j].startpos,dict[j].startpos+dict[j].numkeys+1,0);
		std::ifstream finkey(outdir+std::string("keys.bin.")+std::to_string(j),std::ifstream::in|std::ios::binary);
		uint64_t currentkey;
		for(uint32_t i = 0; i < numreads; i++)
		{
			finkey.read((char*)&currentkey, sizeof(uint64_t));
			dict[j].startpos[((dict[j].bphf)->lookup(currentkey))+1]++;	
		}
		dict[j].empty_bin = new bool[dict[j].numkeys]();
		for(uint32_t i = 1; i < dict[j].numkeys; i++)
			dict[j].startpos[i] =  dict[j].startpos[i] +  dict[j].startpos[i-1];

		//insert elements in the dict array
		dict[j].read_id = new uint32_t[numreads];
		finkey.seekg(0);
		for(uint32_t i = 0; i < numreads; i++)
		{
			finkey.read((char*)&currentkey, sizeof(uint64_t));
			auto idx = (dict[j].bphf)->lookup(currentkey);
			dict[j].read_id[dict[j].startpos[idx]++] = i;
		}
		finkey.close();
		remove((outdir+std::string("keys.bin.")+std::to_string(j)).c_str());
		
		//correcting startpos array modified during insertion
		for(int64_t i = dict[j].numkeys; i >= 1 ; i--)	
			dict[j].startpos[i] = dict[j].startpos[i-1];
		dict[j].startpos[0] = 0;
	}//for end
		
	}//parallel end
	
	omp_set_num_threads(num_thr);
	return;
}

void bbhashdict::findpos(uint32_t *dictidx, uint32_t &startposidx)
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

void bbhashdict::remove(uint32_t *dictidx, uint32_t &startposidx, uint32_t current)
{
	auto size = dictidx[1] - dictidx[0];
	if(size == 1)//just one read left in bin
	{
		empty_bin[startposidx] = 1;
		return; //need to keep one read to check during matching
	}
	uint32_t pos = std::lower_bound(read_id+dictidx[0],read_id+dictidx[1],current)-(read_id+dictidx[0]);
	
	std::move(read_id+dictidx[0]+pos+1,read_id+dictidx[1],read_id+dictidx[0]+pos);
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

void reorder(std::bitset<2*readlen> *read, uint32_t *read1, uint32_t* read2, bbhashdict *dict)
{
	omp_lock_t dict_lock[numdict];//one lock for each dict
	for(int j = 0; j < numdict; j++)
		omp_init_lock(&dict_lock[j]);
	uint8_t _readlen = readlen;//used for writing to binary file
	std::bitset<2*readlen> mask[maxmatch];
	std::bitset<2*readlen> revmask[maxmatch];
	generatemasks(mask,revmask);
	uint32_t maskfirst[maxmatch], masklast[maxmatch], revmaskfirst[maxmatch], revmasklast[maxmatch];
	for (int i = 0; i < maxmatch; i++)
	{
		maskfirst[i] = (mask[i]&first32).to_ullong();
		revmaskfirst[i] = (revmask[i]&first32).to_ullong();
		masklast[i] = ((mask[i]&last32)>>(2*readlen-32)).to_ullong();
		revmasklast[i] = ((revmask[i]&last32)>>(2*readlen-32)).to_ullong();
	}	
	std::bitset<2*readlen> mask1[numdict];
	generateindexmasks(mask1);
	std::atomic<bool> *remainingreads = new std::atomic<bool>[numreads];
	std::fill(remainingreads, remainingreads+numreads,1);
	std::atomic<int64_t> remainingpos(numreads-1);//used for searching next unmatched read when no match is found
	//we go through remainingreads array from behind as that speeds up deletion from bin arrays
	uint32_t beingwritten[numdict][num_thr], beingread[numdict][num_thr];
	for(int j = 0; j <numdict; j++)
		for(int i = 0; i <num_thr; i++)
		{
			beingwritten[j][i] = numreads;//numreads means not writing anywhere right now
			beingread[j][i] = numreads;
		}
	uint32_t firstread = 0,unmatched = 0;
	#pragma omp parallel
	{
	int tid = omp_get_thread_num();
	std::string tid_str = std::to_string(tid);	
	std::ofstream foutRC(outfileRC + '.' + tid_str,std::ofstream::out);
	std::ofstream foutflag(outfileflag + '.' + tid_str,std::ofstream::out);
	std::ofstream foutpos(outfilepos + '.' + tid_str,std::ofstream::out|std::ios::binary);
	std::ofstream foutorder(outfileorder + '.' + tid_str,std::ofstream::out|std::ios::binary);
	std::bitset<2*readlen> ref,revref,b;
	uint32_t ref1,ref2,revref1,revref2;
	int count[4][readlen];
	uint32_t dictidx[2];//to store the start and end index (end not inclusive) in the dict read_id array
	uint32_t startposidx;//index in startpos
	bool flag = 0, done = 0;
	uint32_t current;
	uint64_t ull,foo1=0,bar1=0,foo2=0,bar2=0;
	//flag to check if match was found or not
	
	#pragma omp critical
	{//doing initial setup and first read
		current = firstread;	
		firstread += numreads/omp_get_num_threads();//spread out first read equally
		remainingreads[current] = 0;
		unmatched++;
	}
	#pragma omp barrier
	updaterefcount(read[current],ref,revref,count,true,false,0);
	foutRC << 'd';
	foutorder.write((char*)&current,sizeof(uint32_t));
	foutflag << 0;//for unmatched
	foutpos.write((char*)&_readlen,sizeof(uint8_t));
	uint32_t numdone= 0;
	while(!done)
	{
		numdone++;
		if(numdone%1000000==0)
			std::cout<<tid<<":"<<numdone<<"\n";
		//delete the read from the corresponding dictionary bins
		for(int l = 0; l < numdict; l++)
		{	b = read[current]&mask1[l];
			ull = (b>>2*dict_start[l]).to_ullong();
			startposidx = dict[l].bphf->lookup(ull);
			//check if any other thread is modifying same dictpos
			int i;
			bool go_on = 0; 
			while(!go_on)//make sure no other thread is reading or writing to dictbin
			{
				omp_set_lock(&dict_lock[l]);
				for(i=0;i<num_thr;i++)
					if(beingwritten[l][i]==startposidx || beingread[l][i]==startposidx)
						break;
				if(i==num_thr)
				{
					beingwritten[l][tid] = startposidx;
					go_on = 1;
				}
				omp_unset_lock(&dict_lock[l]);
			}
			dict[l].findpos(dictidx,startposidx);
			dict[l].remove(dictidx,startposidx,current);

			#pragma omp atomic write
			beingwritten[l][tid] = numreads;
		}
		flag = 0;
		uint32_t k;
		for(int j = 0; j < maxmatch; j++)
		{
			//find forward match
			for(int l = 0; l < numdict; l++)
			{
				if(dict_end[l]+j >= readlen)
					continue;
				b = ref&mask1[l];
				ull = (b>>2*dict_start[l]).to_ullong();
				startposidx = dict[l].bphf->lookup(ull);
				if(startposidx >= dict[l].numkeys)//not found
					continue;
				//check if any other thread is modifying same dictpos
				int i;
				bool go_on = 0; 
				while(!go_on)//make sure no other thread is reading or writing to dictbin
				{
					omp_set_lock(&dict_lock[l]);
					for(i=0;i<num_thr;i++)
						if(beingwritten[l][i]==startposidx)
							break;
					if(i==num_thr)
					{
						beingread[l][tid] = startposidx;
						go_on = 1;
					}
					omp_unset_lock(&dict_lock[l]);
				}
				dict[l].findpos(dictidx,startposidx);
				if(dict[l].empty_bin[startposidx])//bin is empty
				{
					#pragma omp atomic write
					beingread[l][tid] = numreads;
					continue;
				}
				uint64_t ull1 = ((read[dict[l].read_id[dictidx[0]]]&mask1[l])>>2*dict_start[l]).to_ullong();
				if(ull == ull1)//checking if ull is actually the key for this bin
				{	
					ref1 = (ref&first32).to_ullong();
					ref2 = ((ref&last32)>>(2*readlen-32)).to_ullong();
					for (int64_t i = dictidx[1] - 1 ; i >= dictidx[0] ; i--)
					{
						auto rid = dict[l].read_id[i];
						if(hamming32[ref1^(read1[rid]&maskfirst[j])]+hamming32[ref2^(read2[rid]&masklast[j])] <= thresh)
						{foo1++;if((ref^(read[rid]&mask[j])).count()<=thresh)
						{	
							#pragma omp critical (readfound)
							if(remainingreads[rid])
							{
								remainingreads[rid]=0;
								k = rid;
								flag = 1;
							}
							if(flag == 1)
								break;
						}}else bar1++;
					}
				}
				#pragma omp atomic write
				beingread[l][tid] = numreads;
				
				if(flag == 1)
				{
					current = k;
					updaterefcount(read[current],ref,revref,count,false,false,j);
					foutRC << 'd';
					foutorder.write((char*)&current,sizeof(uint32_t));
					foutflag << 1;//for matched
					foutpos.write((char*)&j,sizeof(uint8_t));
					break;
				}
				
			}
			if(flag==1)
				break;	

			//find reverse match
			for(int l = 0; l < numdict; l++)
			{
				if(dict_start[l] <= j)
					continue;
				b = revref&mask1[l];
				ull = (b>>2*dict_start[l]).to_ullong();
				startposidx = dict[l].bphf->lookup(ull);
				if(startposidx >= dict[l].numkeys)//not found
					continue;
				//check if any other thread is modifying same dictpos
				int i;
				bool go_on = 0; 
				while(!go_on)//make sure no other thread is reading or writing to dictbin
				{
					omp_set_lock(&dict_lock[l]);
					for(i=0;i<num_thr;i++)
						if(beingwritten[l][i]==startposidx)
							break;
					if(i==num_thr)
					{
						beingread[l][tid] = startposidx;
						go_on = 1;
					}
					omp_unset_lock(&dict_lock[l]);
				}
				dict[l].findpos(dictidx,startposidx);
				if(dict[l].empty_bin[startposidx])//bin is empty
				{	
					#pragma omp atomic write
					beingread[l][tid] = numreads;
					continue;
				}
				uint64_t ull1 = ((read[dict[l].read_id[dictidx[0]]]&mask1[l])>>2*dict_start[l]).to_ullong();
				if(ull == ull1)//checking if ull is actually the key for this bin
				{
					revref1 = (revref&first32).to_ullong();
					revref2 = ((revref&last32)>>(2*readlen-32)).to_ullong();
					for (int64_t i = dictidx[1] - 1 ; i >= dictidx[0] ; i--)
					{
						auto rid = dict[l].read_id[i];
						if(hamming32[revref1^(read1[rid]&revmaskfirst[j])]+hamming32[revref2^(read2[rid]&revmasklast[j])] <= thresh)
						{foo2++;if((revref^(read[rid]&revmask[j])).count()<=thresh)
						{	
							#pragma omp critical (readfound)
							if(remainingreads[rid])
							{
								remainingreads[rid]=0;
								k = rid;
								flag = 1;
							}
							if(flag == 1)
								break;
						}}else bar2++;
					}
				}
				#pragma omp atomic write
				beingread[l][tid] = numreads;
				if(flag == 1)
				{
					current = k;
					updaterefcount(read[current],ref,revref,count,false,true,j);
					foutRC << 'r';
					foutorder.write((char*)&current,sizeof(uint32_t));
					foutflag << 1;//for matched
					foutpos.write((char*)&j,sizeof(uint8_t));
					break;
				}
			}	
			if(flag==1)
				break;
		
			revref<<=2;
			ref>>=2;
		}
		if(flag == 0)//no match found
		{
			for(int64_t j = remainingpos; j>=0; j--)
			
				if(remainingreads[j] == 1)
				{
					#pragma omp critical (readfound)	
					if(remainingreads[j])//checking again inside critical block
					{
						current = j;
						remainingpos = j-1;
						remainingreads[j] = 0;
						flag = 1;
						unmatched++;
					}
					if(flag == 1)
						break;
				}
					
			if(flag == 0)
				done = 1;//no reads left
			else
			{
				updaterefcount(read[current],ref,revref,count,true,false,0);
				foutRC << 'd';
				foutorder.write((char*)&current,sizeof(uint32_t));
				foutflag << 0;//for unmatched
				foutpos.write((char*)&_readlen,sizeof(uint8_t));
			}
		}
	}//while(!done) end
	
	foutRC.close();
	foutorder.close();
	foutflag.close();
	foutpos.close();
	#pragma omp critical
	{std::cout << tid << ":Done"<<"\n";
	std::cout << "foo1 " << foo1 <<" bar1 " << bar1 << "\n";
	std::cout << "foo2 " << foo2 <<" bar2 " << bar2 << "\n";}
	}//parallel end
		
	delete[] remainingreads;
		
	std::cout << "Reordering done, "<<unmatched<<" were unmatched\n";

	return;
}


void generatemasks(std::bitset<2*readlen> *mask,std::bitset<2*readlen> *revmask)
{
	for(int i = 0; i < maxmatch; i++)
	{	
		mask[i].reset();
		revmask[i].reset();
		for(int j = 0; j < 2*readlen - 2*i; j++)
			mask[i][j] = 1;
		for(int j = 2*i; j < 2*readlen; j++)
			revmask[i][j] = 1; 	
	}
	return;
}



void writetofile(std::bitset<2*readlen> *read)
{

	//convert bitset to string for all num_thr files in parallel
	#pragma omp parallel
	{
	int tid = omp_get_thread_num();	
	std::string tid_str = std::to_string(tid);
	std::ofstream fout(outfile + '.' + tid_str,std::ofstream::out);
	std::ifstream finRC(outfileRC + '.' + tid_str,std::ifstream::in);
	std::ifstream finorder(outfileorder + '.' + tid_str,std::ifstream::in|std::ios::binary);
	char s[readlen+1],s1[readlen+1];
	s[readlen] = '\0';
	s1[readlen] = '\0';
	uint32_t current;
	char c;
	while(finRC >> std::noskipws >> c)//read character by character
	{
		finorder.read((char*)&current, sizeof (uint32_t));
		bitsettostring(read[current],s);
		if(c == 'd')
		{
			fout<<s<<"\n";
		}
		else
		{
			reverse_complement(s,s1);
			fout<<s1<<"\n";
		}
	}
	
	fout.close();
	finRC.close();
	finorder.close();	

	}

	//Now combine the num_thr files
	std::ofstream fout(outfile,std::ofstream::out);
	std::ofstream foutRC(outfileRC,std::ofstream::out);
	std::ofstream foutflag(outfileflag,std::ofstream::out);
	std::ofstream foutpos(outfilepos,std::ofstream::out|std::ios::binary);
	std::ofstream foutorder(outfileorder,std::ofstream::out|std::ios::binary);
	for(int tid = 0; tid < num_thr; tid++)
	{
		std::string tid_str = std::to_string(tid);
		std::ifstream fin(outfile + '.' + tid_str,std::ifstream::in);
		std::ifstream finRC(outfileRC + '.' + tid_str,std::ifstream::in);
		std::ifstream finflag(outfileflag + '.' + tid_str,std::ifstream::in);
		std::ifstream finpos(outfilepos + '.' + tid_str,std::ifstream::in|std::ios::binary);
		std::ifstream finorder(outfileorder + '.' + tid_str,std::ifstream::in|std::ios::binary);
		
		fout << fin.rdbuf();//write entire file
		foutRC << finRC.rdbuf();
		foutflag << finflag.rdbuf();
		foutpos << finpos.rdbuf();
		foutorder << finorder.rdbuf();
		
		fin.close();
		finRC.close();
		finflag.close();
		finorder.close();
		finpos.close();
		
		remove((outfile + '.' + tid_str).c_str());
		remove((outfileRC + '.' + tid_str).c_str());
		remove((outfileflag + '.' + tid_str).c_str());
		remove((outfilepos + '.' + tid_str).c_str());
		remove((outfileorder + '.' + tid_str).c_str());
	}
	fout.close();
	foutRC.close();
	foutflag.close();
	foutorder.close();
	foutpos.close();
	return;
}

void bitsettostring(std::bitset<2*readlen> b,char *s)
{
	unsigned long long ull,rem;
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
	return;
}

std::bitset<2*readlen> chartobitset(char *s)
{
	std::bitset<2*readlen> b;
	for(int i = 0; i < readlen; i++)
		b |= basemask[i][s[i]];
	return b;
}

void reverse_complement(char* s, char* s1)
{
	for(int j = 0; j < readlen; j++)
		s1[j] = chartorevchar[s[readlen-j-1]];
	return;
}

void updaterefcount(std::bitset<2*readlen> cur, std::bitset<2*readlen> &ref, std::bitset<2*readlen> &revref, int count[][readlen], bool resetcount, bool rev, int shift)
{
	char s[readlen+1],s1[readlen+1],*current;
	bitsettostring(cur,s);
	if(rev == false)
		current = s;
	else
	{
		reverse_complement(s,s1);
		current = s1;
	}

	if(resetcount == true)//resetcount - unmatched read so start over
	{
		std::memset(count, 0, sizeof(count[0][0]) * 4 * readlen);	
		for(int i = 0; i < readlen; i++)
		{	
			count[chartoint[current[i]]][i] = 1;
		}

	}
	else
	{
		//shift count and find current by max
		for(int i = 0; i < readlen-shift; i++)
		{	
			for(int j = 0; j < 4; j++)
				count[j][i] = count[j][i+shift];
			count[chartoint[current[i]]][i] += 1;
			
			int max = 0,indmax = 0;
			for(int j = 0; j < 4; j++)
				if(count[j][i]>max)
				{
					max = count[j][i];
					indmax = j;
				}
			current[i] = inttochar[indmax];
		}
		//for the new positions make current same as the current read and count to 1
		for(int i = readlen-shift; i < readlen; i++)
		{	
			for(int j = 0; j < 4; j++)
				count[j][i] = 0;
			count[chartoint[current[i]]][i] = 1;
		}
	}
	ref = chartobitset(current);
	char revcurrent[readlen];
	reverse_complement(current,revcurrent);
	revref = chartobitset(revcurrent);
	return;
}
