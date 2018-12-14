//v1 mentioned in the report. Reordering without using any clean reference. 
//First union of all dictionary bins is taken and then a read in a bin picked if hamming distance is below threshold.
//generates three output files: outfile has the reordered reads (some might be reverse complemented),
//outfileRC has 0's and 1's where 0 means that read is not RC and 1 means it is, 
//outfileflag has information which can be useful to the second stage - 
//it contains 0 or 1 depending on whether the read is a matched or not
//(so second stage immediately realizes if a read is unmatched)

//Note that reads with N are currently not supported


#include <iostream>
#include <fstream>
#include <bitset>
#include <unordered_set>
#include <unordered_map>
#include <string>
#include <vector>
#include <algorithm>
#include <set>

#define infile "SRR327342_2_clean.dna"
#define outfile "temp7.dna"
#define outfileRC "tempRC7.txt"
#define outfileflag "tempflag7.txt"
#define readlen 75
#define maxmatch 23
#define numreads 14283592
#define thresh 16
#define numdict 2

void generateindexmasks(std::bitset<2*readlen> *mask1)
//function to generate dictionary index masks
//should be symmetric about readlen (e.g. for 2 dicts - if first dict is start1-end1 (both included),
//then second should be (readlen-1-end1)-(readlen-1-start1))
{
	for(int i = 0; i < numdict; i++)
		mask1[i].reset();
	for(int i = 2*23; i < 2*38; i++)
		mask1[0][i] = 1;
	for(int i = 2*37; i < 2*52; i++)
		mask1[1][i] = 1;
	return;
}

void stringtobitset(std::string s,std::bitset<2*readlen> &read, std::bitset<2*readlen> &revread);

std::string bitsettostring(std::bitset<2*readlen> b);

void readDnaFile(std::bitset<2*readlen> *read,std::bitset<2*readlen> *revread);

void constructdictionary(std::bitset<2*readlen> *read, std::unordered_map<std::bitset<2*readlen>,std::vector<int>> *dict);

void generatemasks(std::bitset<2*readlen> *mask,std::bitset<2*readlen> *revmask);

void reorder(std::bitset<2*readlen> *read, std::bitset<2*readlen> *revread, std::unordered_map<std::bitset<2*readlen>,std::vector<int>> *dict,std::vector<int> &sortedorder,std::vector<int> &revcomp,std::vector<int>& flagvec);

void writetofile(std::bitset<2*readlen> *read,std::bitset<2*readlen> *revread, std::vector<int> &sortedorder,std::vector<int> &revcomp,std::vector<int> &flagvec);


int main()
{
	std::bitset<2*readlen> *read = new std::bitset<2*readlen> [numreads];
	std::bitset<2*readlen> *revread = new std::bitset<2*readlen> [numreads];
	std::cout << "Reading file\n";
	//using vector instead of list to save some space (at the cost of linear time required to delete elements)
	readDnaFile(read,revread);
	std::cout << "Constructing dictionaries\n";
	std::unordered_map<std::bitset<2*readlen>,std::vector<int>> *dict = new std::unordered_map<std::bitset<2*readlen>,std::vector<int>> [numdict];
	constructdictionary(read,dict);
	std::vector<int> sortedorder,revcomp,flagvec;
	std::cout << "Reordering reads\n";
	reorder(read,revread,dict,sortedorder,revcomp,flagvec);
	std::cout << "Writing to file\n";
	writetofile(read,revread,sortedorder,revcomp,flagvec);	
	std::cout << "Done!\n";
	return 0;
}


void stringtobitset(std::string s,std::bitset<2*readlen> &b, std::bitset<2*readlen> &b1)
{
	int i;
	for(int i = 0; i < readlen; i++)
	{	
		switch(s[i])
		{
			case 'A':	b[2*i] = 0;
					b[2*i+1] = 0;
					b1[2*(readlen-i-1)] = 1;
					b1[2*(readlen-i-1)+1] = 1;
					break;
			case 'C':	b[2*i] = 0;
					b[2*i+1] = 1;
					b1[2*(readlen-i-1)] = 1;
					b1[2*(readlen-i-1)+1] = 0;
					break;
			case 'G':	b[2*i] = 1;
					b[2*i+1] = 0;
					b1[2*(readlen-i-1)] = 0;
					b1[2*(readlen-i-1)+1] = 1;
					break;
			case 'T':	b[2*i] = 1;
					b[2*i+1] = 1;
					b1[2*(readlen-i-1)] = 0;
					b1[2*(readlen-i-1)+1] = 0;
					break;
		}
	}
	return;
}

void readDnaFile(std::bitset<2*readlen> *read, std::bitset<2*readlen> *revread)
{
	std::ifstream f(infile, std::ifstream::in);
	f.seekg(0, f.beg);
	std::string s;
	for(int i = 0; i < numreads; i++)
	{
		if(i%1000000 == 0)
			std::cout << i/1000000 << "M reads read"<<"\n";
		f >> s;
		stringtobitset(s,read[i],revread[i]);
	}
	f.close();
	return;
}
		

void constructdictionary(std::bitset<2*readlen> *read, std::unordered_map<std::bitset<2*readlen>,std::vector<int>> *dict)
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

void reorder(std::bitset<2*readlen> *read, std::bitset<2*readlen> *revread, std::unordered_map<std::bitset<2*readlen>,std::vector<int>> *dict, std::vector<int> &sortedorder,std::vector<int> &revcomp,std::vector<int> &flagvec)
{	
	std::bitset<2*readlen> *mask = new std::bitset<2*readlen> [maxmatch];
	std::bitset<2*readlen> *revmask = new std::bitset<2*readlen> [maxmatch];
	generatemasks(mask,revmask);
	std::bitset<2*readlen> *mask1 = new std::bitset<2*readlen> [numdict];
	generateindexmasks(mask1);
	std::unordered_set<int> remainingreads;

	for(int i = 0; i < numreads; i++)
		remainingreads.insert(i);
	int unmatched = 0;
	int current = 0;
	int flag,revflag = 0;
	//flag to check if match was found or not, revflag to check if the current read is forward or reverse 
	sortedorder.push_back(current);
	revcomp.push_back(0);//for direct
	flagvec.push_back(0);//for unmatched
	std::bitset<2*readlen> *b = new std::bitset<2*readlen> [numdict];
	std::bitset<2*readlen> b1,b2;
	
	for(int i = 1; i < numreads; i++)
	{
		if(i%1000000 == 0)
			std::cout<< i/1000000 << "M reads done"<<"\n";
		remainingreads.erase(current);
		for(int l = 0; l < numdict; l++)
		{	b[l] = read[current]&mask1[l];
			int pos = std::lower_bound(dict[l][b[l]].begin(),dict[l][b[l]].end(),current)-dict[l][b[l]].begin();
			//binary search since dict[b] is sorted
			dict[l][b[l]].erase(dict[l][b[l]].begin()+pos);
			if(dict[l][b[l]].size() == 0)
				dict[l].erase(b[l]);
		}
		if (revflag == 0)	
		{	b1 = read[current];
			b2 = revread[current];
		}
		else
		{
			b1 = revread[current];
			b2 = read[current];
		}
		flag = 0;
		for(int j = 0; j < maxmatch; j++)
		{
			std::set<int> s;
			std::vector<int>::iterator it;
			for(int l = 0; l < numdict; l++)
			{
				b[l] = b1&mask1[l];
				if(dict[l].count(b[l]) == 1)
					s.insert(dict[l][b[l]].begin(),dict[l][b[l]].end());
			}
			if(s.size() > 0)
			{
				for (std::set<int>::iterator it = s.begin() ; it != s.end(); ++it)
				{
					int k = *it;
					if((b1^(read[k]&mask[j])).count()<=thresh)
					{
						current = k;
						flag = 1;
						revflag = 0;
						revcomp.push_back(revflag);
						sortedorder.push_back(current);
						flagvec.push_back(1);//for matched
						break;
					}			
					
				}
				if(flag == 1)
					break;
			}
			b1>>=2;

			std::set<int> s1;
			for(int l = 0; l < numdict; l++)
			{
				b[l] = b2&mask1[l];
				if(dict[l].count(b[l]) == 1)
					s1.insert(dict[l][b[l]].begin(),dict[l][b[l]].end());
			}
			if(s1.size() > 0)
			{
				for (std::set<int>::iterator it = s1.begin() ; it != s1.end(); ++it)
				{
					int k = *it;
					if((b2^(read[k]&revmask[j])).count()<=thresh)
					{
						current = k;
						flag = 1;
						revflag = 1;
						revcomp.push_back(revflag);
						sortedorder.push_back(current);
						flagvec.push_back(1);
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
			revflag = 0;
			revcomp.push_back(revflag);
			sortedorder.push_back(current);
			flagvec.push_back(0);
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



void writetofile(std::bitset<2*readlen> *read,std::bitset<2*readlen> *revread, std::vector<int> &sortedorder,std::vector<int> &revcomp,std::vector<int> &flagvec)
{
	std::ofstream fout(outfile,std::ofstream::out);
	std::ofstream foutRC(outfileRC,std::ofstream::out);
	std::ofstream foutflag(outfileflag,std::ofstream::out);
	std::vector<int>::iterator it1,it2,it3;
	for (it1 = sortedorder.begin(),it2 = revcomp.begin(),it3 = flagvec.begin() ; it1 != sortedorder.end(); ++it1,++it2,++it3)
	{
		foutflag << *it3;
		if(*it2 == 0)
		{
			fout<<bitsettostring(read[*it1])<<"\n";
			foutRC << 'd';
		}
		else
		{
			fout<<bitsettostring(revread[*it1])<<"\n";
			foutRC << 'r';
		}
	}
	fout.close();
	foutRC.close();
	foutflag.close();
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
