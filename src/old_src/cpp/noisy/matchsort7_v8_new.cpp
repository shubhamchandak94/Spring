#include <iostream>
#include <fstream>
#include <bitset>
#include <unordered_set>
#include <sparsepp/spp.h>
#include <string>
#include <vector>
#include <algorithm>
#include <set>
#include <cstring>
#include <string>
#include "config.h"

//#define readlen 100


int numreads = 0;

std::string infile;
std::string outfile;
std::string outfileRC;
std::string outfileflag;
std::string outfilepos;


void generateindexmasks(std::bitset<2*readlen> *mask1)
{
	for(int i = 0; i < numdict; i++)
		mask1[i].reset();
	for(int i = 2*dict1_start; i < 2*(dict1_end+1); i++)
		mask1[0][i] = 1;
	for(int i = 2*dict2_start; i < 2*(dict2_end+1); i++)
		mask1[1][i] = 1;
	return;
}

char inttochar[] = {'A','C','G','T'};
char chartorevchar[128];
int chartoint[128];
char chartobit[128][2];
int charinttoint[128];

std::bitset<2*readlen> stringtobitset(std::string s);

void bitsettostring(std::bitset<2*readlen> b,char *s);

void readDnaFile(std::bitset<2*readlen> *read);

void getDataParams();

void constructdictionary(std::bitset<2*readlen> *read, spp::sparse_hash_map<std::bitset<2*readlen>,std::vector<int>> *dict);

void generatemasks(std::bitset<2*readlen> *mask,std::bitset<2*readlen> *revmask);

void reorder(std::bitset<2*readlen> *read, spp::sparse_hash_map<std::bitset<2*readlen>,std::vector<int>> *dict,std::vector<int> &sortedorder,std::vector<int> &revcomp,std::vector<int>& flagvec, std::vector<int>& readpos);

void writetofile(std::bitset<2*readlen> *read, std::vector<int> &sortedorder,std::vector<int> &revcomp,std::vector<int> &flagvec, std::vector<int>& readpos);

void updaterefcount(std::bitset<2*readlen> cur, std::bitset<2*readlen> &ref, std::bitset<2*readlen> &revref, int count[][readlen], bool resetcount, bool rev, int shift);

void reverse_complement(char *s, char *s1);

std::bitset<2*readlen> chartobitset(char &s);

void setglobalarrays();

int main(int argc, char** argv)
{
	std::string basedir = std::string(argv[1]);
	infile = basedir + "/input_clean.dna";
	outfile = basedir + "/output/temp.dna";
	outfileRC = basedir + "/output/read_rev.txt";
	outfileflag = basedir + "/output/tempflag.txt";
	outfilepos = basedir + "/output/temppos.txt";
	getDataParams(); //populate numreads, readlen
	
	setglobalarrays();
	std::bitset<2*readlen> *read = new std::bitset<2*readlen> [numreads];
	std::cout << "Reading file: " << infile << std::endl;
	readDnaFile(read);
	std::cout << "Constructing dictionaries\n";
	//using vector instead of list to save some space (at the cost of linear time required to delete elements)
	spp::sparse_hash_map<std::bitset<2*readlen>,std::vector<int>> *dict = new spp::sparse_hash_map<std::bitset<2*readlen>,std::vector<int>> [numdict];
	constructdictionary(read,dict);
	std::vector<int> sortedorder,revcomp,flagvec,readpos;
	std::cout << "Reordering reads\n";
	reorder(read,dict,sortedorder,revcomp,flagvec,readpos);
	std::cout << "Writing to file\n";
	writetofile(read,sortedorder,revcomp,flagvec,readpos);	
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
	charinttoint['0'] = 0;
	charinttoint['1'] = 1;
	chartobit['A'][0] = '0';
	chartobit['A'][1] = '0';
	chartobit['C'][0] = '0';
	chartobit['C'][1] = '1';
	chartobit['G'][0] = '1';
	chartobit['G'][1] = '0';
	chartobit['T'][0] = '1';
	chartobit['T'][1] = '1';
	return;
}
	

std::bitset<2*readlen> stringtobitset(std::string s)
{
	char *s2 = &(s[0]);
	char s1[2*readlen];
	for(int i = 0; i < readlen; i++)
	{
		s1[2*(readlen-i-1)+1] = chartobit[s2[i]][0];
		s1[2*(readlen-i-1)] = chartobit[s2[i]][1];
	}
	std::bitset<2*readlen> b(s1);
	return b;
}

void getDataParams()
{
	int number_of_lines = 0;
    std::string line;
    std::ifstream myfile(infile, std::ifstream::in);
    
    std::getline(myfile, line);
    int read_length = line.length();
    std::cout << "Length of file: " << read_length << std::endl;
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
	for(int i = 0; i < numreads; i++)
	{
		if(i%1000000 == 0)
			std::cout << i/1000000 << "M reads read"<<"\n";
		f >> s;
		read[i] = stringtobitset(s);
	}
	f.close();
	return;
}
		

void constructdictionary(std::bitset<2*readlen> *read, spp::sparse_hash_map<std::bitset<2*readlen>,std::vector<int>> *dict)
{
	std::bitset<2*readlen> b;
	std::bitset<2*readlen> *mask = new std::bitset<2*readlen> [numdict];
	generateindexmasks(mask);
	for(int i = 0; i < numreads; i++)
	{	
		if(i%1000000 == 0)
			std::cout << i/1000000 << "M reads done"<<"\n";
		for(int j = 0; j < numdict; j++)
		{
			b = read[i]&mask[j];
			if(dict[j].count(b) == 1)
				dict[j][b].push_back(i);
			else
				dict[j][b] = {i};
		}
	}
	return;
}

void reorder(std::bitset<2*readlen> *read, spp::sparse_hash_map<std::bitset<2*readlen>,std::vector<int>> *dict, std::vector<int> &sortedorder,std::vector<int> &revcomp,std::vector<int> &flagvec, std::vector<int>& readpos)
{	
	std::bitset<2*readlen> *mask = new std::bitset<2*readlen> [maxmatch];
	std::bitset<2*readlen> *revmask = new std::bitset<2*readlen> [maxmatch];
	std::bitset<2*readlen> ref,revref;
	int count[4][readlen];
	generatemasks(mask,revmask);
	std::bitset<2*readlen> *mask1 = new std::bitset<2*readlen> [numdict];
	generateindexmasks(mask1);
	std::unordered_set<int> remainingreads;

	for(int i = 0; i < numreads; i++)
		remainingreads.insert(i);
	int unmatched = 0;
	int current = 0;
	int flag;
	//flag to check if match was found or not, revflag to check if the current read is forward or reverse 
	sortedorder.push_back(current);
	revcomp.push_back(0);//for direct
	flagvec.push_back(0);//for unmatched
	readpos.push_back(readlen);
	updaterefcount(read[current],ref,revref,count,true,false,0);
	std::bitset<2*readlen> *b = new std::bitset<2*readlen> [numdict];
	std::bitset<2*readlen> b1,b2;
	long long hammingcount = 0;
	for(int i = 1; i < numreads; i++)
	{
		if(i%1000000 == 0)
			std::cout<< i/1000000 << "M reads done. Hamming count = "<< hammingcount <<"\n";
		remainingreads.erase(current);
		for(int l = 0; l < numdict; l++)
		{	b[l] = read[current]&mask1[l];
			int pos = std::lower_bound(dict[l][b[l]].begin(),dict[l][b[l]].end(),current)-dict[l][b[l]].begin();
			//binary search since dict[b] is sorted
			dict[l][b[l]].erase(dict[l][b[l]].begin()+pos);
			if(dict[l][b[l]].size() == 0)
				dict[l].erase(b[l]);
		}
			
		b1 = ref;
		b2 = revref;
		
		flag = 0;
		for(int j = 0; j < maxmatch; j++)
		{
			std::vector<int> s;
			for(int l = 0; l < numdict; l++)
			{
				b[l] = b1&mask1[l];
				if(dict[l].count(b[l]) == 1)
				{
					if(s.empty())
						s = dict[l][b[l]];
					else
					{
						std::vector<int> temp;
						std::set_union(s.begin(),s.end(),dict[l][b[l]].begin(),dict[l][b[l]].end(),std::back_inserter(temp));
						s = temp;
					}
				}
			}
			if(!s.empty())
			{
				for (std::vector<int>::iterator it = s.begin() ; it != s.end(); ++it)
				{
					int k = *it;
					hammingcount++;
					if((b1^(read[k]&mask[j])).count()<=thresh)
					{
						current = k;
						flag = 1;
						updaterefcount(read[current],ref,revref,count,false,false,j);
						revcomp.push_back(0);
						sortedorder.push_back(current);
						flagvec.push_back(1);//for matched
						readpos.push_back(j);
						break;
					}			
					
				}
				if(flag == 1)
					break;
			}
			b1>>=2;

			std::vector<int> s1;
			for(int l = 0; l < numdict; l++)
			{
				b[l] = b2&mask1[l];
				if(dict[l].count(b[l]) == 1)
				{
					if(s1.empty())
						s1 = dict[l][b[l]];
					else
					{
						std::vector<int> temp;
						std::set_union(s1.begin(),s1.end(),dict[l][b[l]].begin(),dict[l][b[l]].end(),std::back_inserter(temp));
						s1 = temp;
					}
				}
			}
			if(!s1.empty())
			{
				for (std::vector<int>::iterator it = s1.begin() ; it != s1.end(); ++it)
				{
					int k = *it;
					hammingcount++;
					if((b2^(read[k]&revmask[j])).count()<=thresh)
					{
						current = k;
						flag = 1;
						updaterefcount(read[current],ref,revref,count,false,true,j);
						revcomp.push_back(1);
						sortedorder.push_back(current);
						flagvec.push_back(1);
						readpos.push_back(j);
						break;
					}			
					
				}
				if(flag == 1)
					break;
			}
			b2<<=2;
		}	
		if(flag == 0)
		{
			unmatched += 1;
			current = *remainingreads.begin();
			revcomp.push_back(0);
			updaterefcount(read[current],ref,revref,count,true,false,0);
			sortedorder.push_back(current);
			flagvec.push_back(0);
			readpos.push_back(readlen);
		}
	}
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



void writetofile(std::bitset<2*readlen> *read, std::vector<int> &sortedorder,std::vector<int> &revcomp,std::vector<int> &flagvec, std::vector<int>& readpos)
{
	std::ofstream fout(outfile,std::ofstream::out);
	std::ofstream foutRC(outfileRC,std::ofstream::out);
	std::ofstream foutflag(outfileflag,std::ofstream::out);
	std::ofstream foutpos(outfilepos,std::ofstream::out);
	std::vector<int>::iterator it1,it2,it3,it4;
	char s[readlen],s1[readlen];
	for (it1 = sortedorder.begin(),it2 = revcomp.begin(),it3 = flagvec.begin(),it4 = readpos.begin(); it1 != sortedorder.end(); ++it1,++it2,++it3,++it4)
	{
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
	return;
}

void bitsettostring(std::bitset<2*readlen> b,char *s)
{
	std::string s1 = b.to_string();
	char *s2 = &(s1[0]);
	for(int i = 0; i < readlen; i++)
		s[i] = inttochar[2*charinttoint[s2[2*(readlen-i-1)+1]]+charinttoint[s2[2*(readlen-i-1)]]];
	return;
}

std::bitset<2*readlen> chartobitset(char *s)
{
	char s1[2*readlen];
	for(int i = 0; i < readlen; i++)
	{
		s1[2*(readlen-i-1)+1] = chartobit[s[i]][0];
		s1[2*(readlen-i-1)] = chartobit[s[i]][1];
	}
	std::bitset<2*readlen> b(s1);
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
	char s[readlen],s1[readlen],*current;
	bitsettostring(cur,s);
	if(rev == false)
		current = s;
	else
	{
		reverse_complement(s,s1);
		current = s1;
	}

	if(resetcount == true)
	{
		std::memset(count, 0, sizeof(count[0][0]) * 4 * readlen);	
		for(int i = 0; i < readlen; i++)
		{	
			count[chartoint[current[i]]][i] = 1;
		}

	}
	else
	{
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
