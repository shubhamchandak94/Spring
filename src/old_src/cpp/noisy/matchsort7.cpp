//Reordering for real reads
//Similar to matchsort3 but maintains a majority based clean reference read (function updaterefcount). 
//A count matrix is also maintained which the number of times each base was seen at each position in the current read.
//The reference (or its RC) is used for looking in the dictionary and for hamming distance threshold.
//The function updaterefcount is called after every read, it updates the count matrix and the clean reference.

//This does reduce the singleton reads significantly, however the effect on the overall size is small because the noise files 
//become larger and offset most of the gains in the seq file. The function has a resetcount parameter which is made true when
//we don't a matching read and pick a random read from the remaining reads. In such a situation, the count is first reset to 0
//and then at each position, the count for the corresponding base is made 1.

//The encoding stage (which also uses reference) can be easily combined with this code. However in that case we can't use 
//parallelization for the encoding stage

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

#define infile "SRR870667_1_clean.dna"
#define outfile "temp1.dna"
#define outfileRC "tempRC1.txt"
#define outfileflag "tempflag1.txt"
#define readlen 108
#define maxmatch 40
#define numreads 68266234
#define thresh 16
#define numdict 2

void generateindexmasks(std::bitset<2*readlen> *mask1)
//function to generate dictionary index masks
//should be symmetric about readlen (e.g. for 2 dicts - if first dict is start1-end1 (both included), 
//then second should be (readlen-1-end1)-(readlen-1-start1))
{
	for(int i = 0; i < numdict; i++)
		mask1[i].reset();
	for(int i = 2*34; i < 2*54; i++)
		mask1[0][i] = 1;
	for(int i = 2*54; i < 2*74; i++)
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

void updaterefcount(std::bitset<2*readlen> current, std::bitset<2*readlen> &ref, std::bitset<2*readlen> &revref, int count[][readlen], bool resetcount, int shift);


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
	updaterefcount(read[current],ref,revref,count,true,0);
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
			
		b1 = ref;
		b2 = revref;
		
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
						updaterefcount(read[current],ref,revref,count,false,j);
						revcomp.push_back(0);
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
						updaterefcount(revread[current],ref,revref,count,false,j);
						revcomp.push_back(1);
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
			revcomp.push_back(0);
			updaterefcount(read[current],ref,revref,count,true,0);
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


void updaterefcount(std::bitset<2*readlen> current, std::bitset<2*readlen> &ref, std::bitset<2*readlen> &revref, int count[][readlen], bool resetcount, int shift)
{
	if(resetcount == true)
	{
		ref = current;	
		for(int i = 0; i < readlen; i++)
		{	
			for(int j = 0; j < 4; j++)
				count[j][i] = 0;
			switch(current[2*i])
			{
				case 0:	switch(current[2*i+1])
					{
						case 0:	count[0][i] = 1;
							break;
						case 1:	count[1][i] = 1;
							break;
					}		
					break;
				case 1:	switch(current[2*i+1])
					{
						case 0:	count[2][i] = 1;
							break;
						case 1:	count[3][i] = 1;
							break;
					}		
			}
		}

	}
	else
	{
		for(int i = 0; i < readlen-shift; i++)
		{	
			for(int j = 0; j < 4; j++)
				count[j][i] = count[j][i+shift];
			switch(current[2*i])
			{
				case 0:	switch(current[2*i+1])
					{
						case 0:	count[0][i] += 1;
							break;
						case 1:	count[1][i] += 1;
							break;
					}
					break;		
				case 1:	switch(current[2*i+1])
					{
						case 0:	count[2][i] += 1;
							break;
						case 1:	count[3][i] += 1;
							break;
					}		
			}
			int max = 0,indmax = 0;
			for(int j = 0; j < 4; j++)
				if(count[j][i]>max)
				{
					max = count[j][i];
					indmax = j;
				}
			switch(indmax)
			{
				case 0: ref[2*i] = 0;
					ref[2*i+1] = 0;
					break;			
				case 1: ref[2*i] = 0;
					ref[2*i+1] = 1;
					break;			
				case 2: ref[2*i] = 1;
					ref[2*i+1] = 0;
					break;			
				case 3: ref[2*i] = 1;
					ref[2*i+1] = 1;
					break;			
			}
		}
		
		for(int i = readlen-shift; i < readlen; i++)
		{	
			ref[2*i] = current[2*i];
			ref[2*i+1] = current[2*i+1];
			for(int j = 0; j < 4; j++)
				count[j][i] = 0;
			switch(current[2*i])
			{
				case 0:	switch(current[2*i+1])
					{
						case 0:	count[0][i] = 1;
							break;
						case 1:	count[1][i] = 1;
							break;
					}		
					break;
				case 1:	switch(current[2*i+1])
					{
						case 0:	count[2][i] = 1;
							break;
						case 1:	count[3][i] = 1;
							break;
					}		
			}
		}
	}
	
	for(int j = 0; j < readlen; j++)
	{
		revref[2*j] = 1- ref[2*(readlen-j-1)];
		revref[2*j+1] = 1 - ref[2*(readlen-j-1) + 1];
	}
	return;
}
