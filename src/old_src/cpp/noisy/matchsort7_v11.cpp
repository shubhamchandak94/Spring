#include <iostream>
#include <fstream>
#include <bitset>
#include <sparsepp/spp.h>
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>
#include <string>
#include "config.h"
#include <omp.h>
#include <atomic>
//#define readlen 100
//#define num_thr 8

uint32_t numreads = 0;

std::string infile;
std::string outfile;
std::string outfileRC;
std::string outfileflag;
std::string outfilepos;
std::string outfileorder;

//Some global arrays (some initialized in setglobalarrays())
char revinttochar[4] = {'A','G','C','T'};//used in bitsettostring
char inttochar[] = {'A','C','G','T'};
char chartorevchar[128];//A-T etc for reverse complement
int chartoint[128];//A-0,C-1 etc. used in updaterefcount 
std::bitset<2*readlen> basemask[readlen][128];//bitset for A,G,C,T at each position 
//used in stringtobitset, chartobitset and bitsettostring
std::bitset<2*readlen> positionmask[readlen];//bitset for each position (1 at two bits and 0 elsewhere)
//used in bitsettostring
std::bitset<2*readlen> mask64;//bitset with 64 bits set to 1 (used in bitsettostring for conversion to ullong)

std::bitset<2*readlen> stringtobitset(std::string s);

void bitsettostring(std::bitset<2*readlen> b,char *s);

void readDnaFile(std::bitset<2*readlen> *read);

void getDataParams();

void constructdictionary(std::bitset<2*readlen> *read, spp::sparse_hash_map<uint64_t,uint32_t*> *dict);

void generatemasks(std::bitset<2*readlen> *mask,std::bitset<2*readlen> *revmask);
//mask for zeroing the end bits (needed while reordering to compute Hamming distance between shifted reads)

void reorder(std::bitset<2*readlen> *read, spp::sparse_hash_map<uint64_t,uint32_t*> *dict,std::vector<uint32_t> &sortedorder,std::vector<bool> &revcomp,std::vector<bool>& flagvec, std::vector<uint32_t>& readpos);

void writetofile(std::bitset<2*readlen> *read, std::vector<uint32_t> &sortedorder,std::vector<bool> &revcomp,std::vector<bool> &flagvec, std::vector<uint32_t>& readpos);

void updaterefcount(std::bitset<2*readlen> cur, std::bitset<2*readlen> &ref, std::bitset<2*readlen> &revref, int count[][readlen], bool resetcount, bool rev, int shift);
//update the clean reference and count matrix

void reverse_complement(char *s, char *s1);

std::bitset<2*readlen> chartobitset(char &s);

void setglobalarrays();
//setting arrays like chartoint etc.

int main(int argc, char** argv)
{
	std::string basedir = std::string(argv[1]);
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
	std::cout << "Reading file: " << infile << std::endl;
	readDnaFile(read);
	std::cout << "Constructing dictionaries\n";
	spp::sparse_hash_map<uint64_t,uint32_t*> dict[numdict];
	constructdictionary(read,dict);
	std::vector<uint32_t> sortedorder,readpos;
	std::vector<bool> revcomp,flagvec;
	std::cout << "Reordering reads\n";
	reorder(read,dict,sortedorder,revcomp,flagvec,readpos);
	std::cout << "Writing to file\n";
	writetofile(read,sortedorder,revcomp,flagvec,readpos);	
	delete[] read;
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
	myfile.seekg(0, myfile.beg);

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


void readDnaFile(std::bitset<2*readlen> *read)
{
	std::ifstream f(infile, std::ifstream::in);
	f.seekg(0, f.beg);
	std::string s;
	for(uint32_t i = 0; i < numreads; i++)
	{
		f >> s;
		read[i] = stringtobitset(s);
	}
	f.close();
	return;
}
		
void generateindexmasks(std::bitset<2*readlen> *mask1)
//masks for dictionary positions
{
	for(int i = 0; i < numdict; i++)
		mask1[i].reset();
	for(int i = 2*dict1_start; i < 2*(dict1_end+1); i++)
		mask1[0][i] = 1;
	for(int i = 2*dict2_start; i < 2*(dict2_end+1); i++)
		mask1[1][i] = 1;
	return;
}


void constructdictionary(std::bitset<2*readlen> *read, spp::sparse_hash_map<uint64_t,uint32_t*> *dict)
{

	std::bitset<2*readlen> mask[numdict];
	generateindexmasks(mask);
	int dict_start[2] = {dict1_start,dict2_start};
	//Parallelizing construction of the multiple dictionaries
	#pragma omp parallel num_threads(numdict)
	{
	#pragma omp for
	for(int j = 0; j < numdict; j++)
	{	
		std::bitset<2*readlen> b;
		uint64_t ull;
		//find number of times each key occurs
		for(uint32_t i = 0; i < numreads; i++)
		{
			b = read[i]&mask[j];
			ull = (b>>2*dict_start[j]).to_ullong();
			if(dict[j].count(ull) == 1)
				(*dict[j][ull])++;
			else
			{
				dict[j][ull]=new uint32_t;
				(*dict[j][ull]) = 1;
			}
		}
		//allocate memory for each bin (number of reads with given key + 1) 1 for storing the length
		for(auto it = dict[j].begin(); it !=  dict[j].end(); ++it)
		{
			uint32_t binsize = *(it->second);
			delete it->second;
			dict[j][it->first] =  new uint32_t[binsize+1];
			dict[j][it->first][0] = 1;
		}
		//fill in the read ids in each bin, dict[j][ull][0] stores the position where next id is put - at the
		//end it stores the size of the array
		for(uint32_t i = 0; i < numreads; i++)
		{
			b = read[i]&mask[j];
			ull = (b>>2*dict_start[j]).to_ullong();
			dict[j][ull][dict[j][ull][0]++] = i;
		}

	}
	}
	return;
}

void reorder(std::bitset<2*readlen> *read, spp::sparse_hash_map<uint64_t,uint32_t*> *dict, std::vector<uint32_t> &sortedorder,std::vector<bool> &revcomp,std::vector<bool> &flagvec, std::vector<uint32_t>& readpos)
{
	std::bitset<2*readlen> mask[maxmatch];
	std::bitset<2*readlen> revmask[maxmatch];
	generatemasks(mask,revmask);
	std::bitset<2*readlen> mask1[numdict];
	generateindexmasks(mask1);
	std::atomic<bool> *remainingreads = new std::atomic<bool>[numreads];
	std::fill(remainingreads, remainingreads+numreads,1);
	std::atomic<int64_t> remainingpos(numreads-1);//used for searching next unmatched read when no match is found
	//we go through remainingreads array from behind as that speeds up deletion from bin arrays
	int dict_start[2] = {dict1_start,dict2_start};
	std::vector<uint32_t> sortedorder_thr[num_thr], readpos_thr[num_thr];
	std::vector<bool> revcomp_thr[num_thr], flagvec_thr[num_thr];
	uint32_t *beingwritten[num_thr] = {}, *beingread[num_thr] = {};
	uint32_t firstread = 0;
	#pragma omp parallel
	{	
	std::bitset<2*readlen> ref,revref,b;
	int count[4][readlen];
	bool flag = 0, done = 0;
	int tid = omp_get_thread_num();
	uint32_t current;
	uint64_t ull;
	//flag to check if match was found or not
	
	#pragma omp critical
	{//doing initial setup and first read
		current = firstread;	
		firstread += numreads/omp_get_num_threads();//spread out first read equally
		remainingreads[current] = 0;
	}
	#pragma omp barrier	
	updaterefcount(read[current],ref,revref,count,true,false,0);
	revcomp_thr[tid].push_back(0);
	sortedorder_thr[tid].push_back(current);
	flagvec_thr[tid].push_back(0);//for unmatched
	readpos_thr[tid].push_back(readlen);
	uint32_t numdone= 0;
	while(!done)
	{	numdone++;
		if(numdone%1000000==0)
			std::cout<<tid<<":"<<numdone<<"\n";
		//delete the read from the corresponding dictionary bins
		for(int l = 0; l < numdict; l++)
		{	b = read[current]&mask1[l];
			ull = (b>>2*dict_start[l]).to_ullong();
			auto dictbin = dict[l][ull];
			//check if any other thread is modifying dictbin
			int i;
			bool go_on = 0; 
			while(!go_on)//make sure no other thread is reading or writing to dictbin
			{
				#pragma omp critical (dictbin)
				{
					for(i=0;i<num_thr;i++)
						if(beingwritten[i]==dictbin || beingread[i]==dictbin)
							break;
					if(i==num_thr)
					{
						beingwritten[tid] = dictbin;
						go_on = 1;
					}
				}
			}						
			uint32_t pos = std::lower_bound(dictbin+1,dictbin+dictbin[0],current)-dictbin; 
			//binary search since dict[l][ull] is sorted array
			std::move(dictbin+pos+1,dictbin+dictbin[0],dictbin+pos);
			dictbin[0]--;//decrement number of elements
			#pragma omp atomic write
			beingwritten[tid] = NULL;
		}
		flag = 0;
		uint32_t k;
		for(int j = 0; j < maxmatch; j++)
		{
			//find forward match
			for(int l = 0; l < numdict; l++)
			{
				b = ref&mask1[l];
				ull = (b>>2*dict_start[l]).to_ullong();
				if(dict[l].count(ull) == 1)
				{
					auto s = dict[l][ull];
					int i;
					bool go_on = 0;
					while(!go_on)//make sure no other thread is writing to s
					{
						#pragma omp critical (dictbin)
						{
							for(i=0;i<num_thr;i++)
								if(beingwritten[i]==s)
									break;
							if(i==num_thr)
							{
								beingread[tid] = s;
								go_on = 1;
							}
						}
					}						
					long N = s[0]-1;
					for (long i = N ; i >= 1; i--)
						if((ref^(read[s[i]]&mask[j])).count()<=thresh)
						{	
							#pragma omp critical (readfound)
							if(remainingreads[s[i]])
							{
								remainingreads[s[i]]=0;
								k = s[i];
								flag = 1;
							}
							if(flag == 1)
								break;
						}
					#pragma omp atomic write
					beingread[tid] = NULL;
				}
				
				if(flag == 1)
				{
					current = k;
					updaterefcount(read[current],ref,revref,count,false,false,j);
					revcomp_thr[tid].push_back(0);
					sortedorder_thr[tid].push_back(current);
					flagvec_thr[tid].push_back(1);//for matched
					readpos_thr[tid].push_back(j);
					break;
				}
				
			}
			if(flag==1)
				break;	

			//find reverse match
			for(int l = 0; l < numdict; l++)
			{
				b = revref&mask1[l];
				ull = (b>>2*dict_start[l]).to_ullong();
				if(dict[l].count(ull) == 1)
				{
					auto s = dict[l][ull];
					int i;
					bool go_on = 0;
					while(!go_on)//make sure no other thread is writing to s
					{
						#pragma omp critical (dictbin)
						{
							for(i=0;i<num_thr;i++)
								if(beingwritten[i]==s)
									break;
							if(i==num_thr)
							{
								beingread[tid] = s;
								go_on = 1;
							}
						}
					}						
					long N = s[0]-1;
					for (long i = N ; i >= 1; i--)
						if((revref^(read[s[i]]&revmask[j])).count()<=thresh)
						{	
							#pragma omp critical (readfound)
							if(remainingreads[s[i]])
							{
								remainingreads[s[i]]=0;
								k = s[i];
								flag = 1;
							}
							if(flag == 1)
								break;
						}
					#pragma omp atomic write
					beingread[tid] = NULL;
				}
				if(flag == 1)
				{
					current = k;
					updaterefcount(read[current],ref,revref,count,false,true,j);
					revcomp_thr[tid].push_back(1);
					sortedorder_thr[tid].push_back(current);
					flagvec_thr[tid].push_back(1);//for matched
					readpos_thr[tid].push_back(j);
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
					}
					if(flag == 1)
						break;
				}
					
			if(flag == 0)
				done = 1;//no reads left
			else
			{
				updaterefcount(read[current],ref,revref,count,true,false,0);
				revcomp_thr[tid].push_back(0);
				sortedorder_thr[tid].push_back(current);
				flagvec_thr[tid].push_back(0);//for unmatched
				readpos_thr[tid].push_back(readlen);
			}
		}
	}//while(!done) end

	std::cout << tid << ":Done"<<"\n";
	}//parallel end

	//delete stuff
	for(int j = 0; j < numdict; j++)
		for(auto it = dict[j].begin(); it !=  dict[j].end(); ++it)
		{
			delete[] it->second;
		}
		
	delete[] remainingreads;
		
	//Combine sortedorder etc. from different threads
	for(int i = 0; i < num_thr; i++)
	{
		sortedorder.insert(sortedorder.end(),sortedorder_thr[i].begin(),sortedorder_thr[i].end());
		flagvec.insert(flagvec.end(),flagvec_thr[i].begin(),flagvec_thr[i].end());
		revcomp.insert(revcomp.end(),revcomp_thr[i].begin(),revcomp_thr[i].end());
		readpos.insert(readpos.end(),readpos_thr[i].begin(),readpos_thr[i].end());
		//free memory
		std::vector<uint32_t>().swap(sortedorder_thr[i]);
		std::vector<uint32_t>().swap(readpos_thr[i]);
		std::vector<bool>().swap(flagvec_thr[i]);
		std::vector<bool>().swap(revcomp_thr[i]);
	}
	uint32_t unmatched = numreads;
	for(uint32_t i = 0; i < numreads; i++)
		unmatched -= flagvec[i];

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



void writetofile(std::bitset<2*readlen> *read, std::vector<uint32_t> &sortedorder,std::vector<bool> &revcomp,std::vector<bool> &flagvec, std::vector<uint32_t>& readpos)
{
	std::ofstream fout(outfile,std::ofstream::out);
	std::ofstream foutRC(outfileRC,std::ofstream::out);
	std::ofstream foutflag(outfileflag,std::ofstream::out);
	std::ofstream foutpos(outfilepos,std::ofstream::out);
	std::ofstream foutorder(outfileorder,std::ofstream::out|std::ios::binary);
	std::vector<uint32_t>::iterator it1,it4;
	std::vector<bool>::iterator it2,it3;
	char s[readlen+1],s1[readlen+1];
	s[readlen] = '\0';
	s1[readlen] = '\0';
	for (it1 = sortedorder.begin(),it2 = revcomp.begin(),it3 = flagvec.begin(),it4 = readpos.begin(); it1 != sortedorder.end(); ++it1,++it2,++it3,++it4)
	{
		
		foutorder.write ((char*)&(*it1), sizeof (uint32_t));
		foutflag << *it3;
		foutpos << *it4 << "\n";
		bitsettostring(read[*it1],s);
		if(*it2 == 0)
		{
			fout<<s<<"\n";
			foutRC << 'd';
		}
		else
		{
			reverse_complement(s,s1);
			fout<<s1<<"\n";
			foutRC << 'r';
		}
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
