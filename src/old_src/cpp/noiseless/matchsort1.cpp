//Reordering stage for noiseless reads. Uses bitset of length 2*matchlen as dictionary index. Splice function is used to 
//get substring of bitset.
//Slower than matchsort2.cpp but uses slightly less memory.

//outfile contains the reordered reads (one on each line)

#include <iostream>
#include <fstream>
#include <bitset>
#include <unordered_set>
#include <unordered_map>
#include <string>
#include <vector>
#include <algorithm>

#define infile "chrom22_50x_noRC_random.dna"
#define outfile "temp.dna"
#define readlen 100
#define maxmatch 20
#define numreads 17500000
#define matchlen 20

void stringtobitset(std::string s,std::bitset<2*readlen> &b);

std::string bitsettostring(std::bitset<2*readlen> b);

void readDnaFile(std::bitset<2*readlen> *read);

void constructdictionary(std::bitset<2*readlen> *read, std::unordered_map<std::bitset<2*matchlen>,std::vector<int>> &dict);

std::bitset<2*matchlen> splice(std::bitset<2*readlen> b, int start, int end); 

void generatemasks(std::bitset<2*readlen> *mask);

void reorder(std::bitset<2*readlen> *read, std::unordered_map<std::bitset<2*matchlen>,std::vector<int>> &dict,std::vector<int> &sortedorder);

void writetofile(std::bitset<2*readlen> *read, std::vector<int> &sortedorder);


int main()
{
	std::bitset<2*readlen> *read = new std::bitset<2*readlen> [numreads];
	std::cout << "Reading file\n";
	//using vector instead of list to save some space (at the cost of linear time required to delete elements)
	readDnaFile(read);
	std::cout << "Constructing dictionaries\n";
	std::unordered_map<std::bitset<2*matchlen>,std::vector<int>> dict;
	dict.reserve(numreads);
	constructdictionary(read,dict);
	std::vector<int> sortedorder;
	std::cout << "Reordering reads\n";
	reorder(read,dict,sortedorder);
	std::cout << "Writing to file\n";
	writetofile(read,sortedorder);	
	std::cout << "Done!\n";
	return 0;
}


void stringtobitset(std::string s,std::bitset<2*readlen> &b)
{
	int i;
	for(int i = 0; i < readlen; i++)
	{	
		switch(s[i])
		{
			case 'A':	b[2*i] = 0;
					b[2*i+1] = 0;
					break;
			case 'C':	b[2*i] = 0;
					b[2*i+1] = 1;
					break;
			case 'G':	b[2*i] = 1;
					b[2*i+1] = 0;
					break;
			case 'T':	b[2*i] = 1;
					b[2*i+1] = 1;
					break;
		}
	}
	return;
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
		stringtobitset(s,read[i]);
	}
	f.close();
	return;
}
		

void constructdictionary(std::bitset<2*readlen> *read, std::unordered_map<std::bitset<2*matchlen>,std::vector<int>> &dict)
{
	std::bitset<2*matchlen> b;
	for(int i = 0; i < numreads; i++)
	{	
		if(i%1000000 == 0)
			std::cout << i/1000000 << "M reads done"<<"\n";
		b = splice(read[i],0,2*matchlen);
		if(dict.count(b) == 1)
			dict[b].push_back(i);
		else
			dict[b] = {i};
	}
	return;
}

std::bitset<2*matchlen> splice(std::bitset<2*readlen> b, int start, int end)
{
	std::bitset<2*matchlen> b1;
	for (int i = start; i < end; i++)
		b1[i-start] = b[i];
	return b1;
}


void reorder(std::bitset<2*readlen> *read, std::unordered_map<std::bitset<2*matchlen>,std::vector<int>> &dict, std::vector<int> &sortedorder)
{	
	std::bitset<2*readlen> *mask = new std::bitset<2*readlen> [maxmatch];
	generatemasks(mask);
	std::unordered_set<int> remainingreads;

	for(int i = 0; i < numreads; i++)
		remainingreads.insert(i);
	int unmatched = 0;
	int current = 0;
	int flag;
	sortedorder.push_back(current);
	std::bitset<2*matchlen> b;
	std::bitset<2*readlen> b1;

	for(int i = 1; i < numreads; i++)
	{
		if(i%1000000 == 0)
			std::cout<<unmatched << i/1000000 << "M reads done"<<"\n";
		remainingreads.erase(current);
		b = splice(read[current],0,2*matchlen);
		int pos = std::lower_bound(dict[b].begin(),dict[b].end(),current)-dict[b].begin();
		//binary search since dict[b] is sorted
		dict[b].erase(dict[b].begin()+pos);
		if(dict[b].size() == 0)
			dict.erase(b);
		b1 = read[current];
		flag = 0;
		for(int j = 0; j < maxmatch; j++)
		{
			b = splice(b1,0,2*matchlen);
			if(dict.count(b) == 1)
			{
				for (std::vector<int>::iterator it = dict[b].begin() ; it != dict[b].end(); ++it)
				{
					int k = *it;
					if(b1 == (read[k]&mask[j]))
					{
						current = k;
						flag = 1;
						sortedorder.push_back(current);
						break;
					}			
					
				}
				if(flag == 1)
					break;
			}
			
			b1>>=2;
		}	
		if(flag == 0)
		{
			unmatched += 1;
			current = *remainingreads.begin();
			sortedorder.push_back(current);
		}
	}
	std::cout << "Reordering done, "<<unmatched<<" were unmatched";
	return;
}
void generatemasks(std::bitset<2*readlen> *mask)
{
	for(int i = 0; i < maxmatch; i++)
	{	
		for(int j = 0; j < 2*readlen - 2*i; j++)
			mask[i][j] = 1;
		for(int j = 2*readlen - 2*i; j < 2*readlen; j++)
			mask[i][j] = 0;
	}
		
	return;
}


void writetofile(std::bitset<2*readlen> *read, std::vector<int> &sortedorder)
{
	std::ofstream fout(outfile,std::ofstream::out);
	for (std::vector<int>::iterator it = sortedorder.begin() ; it != sortedorder.end(); ++it)
		fout<<bitsettostring(read[*it])<<"\n";
	fout.close();
	return;
}

std::string bitsettostring(std::bitset<2*readlen> b)
{
	std::string s;
	for(int i = 0; i < readlen; i++)
	{
		switch(b[2*i])
		{
			case 0:	switch(b[2*i+1])
				{
					case 0: s.push_back('A');
						break;
					case 1: s.push_back('C');
						break;
				}
				break;
			case 1:	switch(b[2*i+1])
				{
					case 0: s.push_back('G');
						break;
					case 1: s.push_back('T');
						break;
				}
				break;
		}
	}
	return s;
}

	
