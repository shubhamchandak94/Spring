//Reordering for real reads. First union of the bins is taken. After this for each shift, 
//the best match in the union is found, if the bestmatch is below threshold, then we pick the best matching read, 
//otherwise move to next shift/RC.
//Was intended to reduce noise for high thresholds. However the benefit is negligible.

//Two output files: outfile has the reordered files (one in each line) and outfileRC has 0 (for no RC) and 1 (RC)

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

#define infile "SRR065390_clean.dna"
#define outfile "temp2.dna"
#define outfileRC "tempRC2.txt"
#define readlen 100
#define maxmatch 20
#define numreads 67155743
#define thresh 16
#define numdict 4


void generateindexmasks(std::bitset<2*readlen> *mask1)
//function to generate dictionary index masks
//should be symmetric about readlen (e.g. for 2 dicts - if first dict is start1-end1 (both included), 
//then second should be (readlen-1-end1)-(readlen-1-start1))
{
	for(int i = 0; i < numdict; i++)
		mask1[i].reset();
	for(int i = 2*20; i < 2*35; i++)
		mask1[0][i] = 1;
	for(int i = 2*35; i < 2*50; i++)
		mask1[1][i] = 1;
	for(int i = 2*50; i < 2*65; i++)
		mask1[2][i] = 1;
	for(int i = 2*65; i < 2*80; i++)
		mask1[3][i] = 1;
	return;
}


void stringtobitset(std::string s,std::bitset<2*readlen> &read, std::bitset<2*readlen> &revread);

std::string bitsettostring(std::bitset<2*readlen> b);

void readDnaFile(std::bitset<2*readlen> *read,std::bitset<2*readlen> *revread);

void constructdictionary(std::bitset<2*readlen> *read, std::unordered_map<std::bitset<2*readlen>,std::vector<int>> *dict);

void generatemasks(std::bitset<2*readlen> *mask,std::bitset<2*readlen> *revmask);

void reorder(std::bitset<2*readlen> *read, std::bitset<2*readlen> *revread, std::unordered_map<std::bitset<2*readlen>,std::vector<int>> *dict,std::vector<int> &sortedorder,std::vector<int> &revcomp);

void writetofile(std::bitset<2*readlen> *read,std::bitset<2*readlen> *revread, std::vector<int> &sortedorder,std::vector<int> &revcomp);


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
	std::vector<int> sortedorder,revcomp;
	std::cout << "Reordering reads\n";
	reorder(read,revread,dict,sortedorder,revcomp);
	std::cout << "Writing to file\n";
	writetofile(read,revread,sortedorder,revcomp);	
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

void reorder(std::bitset<2*readlen> *read, std::bitset<2*readlen> *revread, std::unordered_map<std::bitset<2*readlen>,std::vector<int>> *dict, std::vector<int> &sortedorder,std::vector<int> &revcomp)
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
				int bestmatch = 2*readlen,bestmatchpos;
				for (std::set<int>::iterator it = s.begin() ; it != s.end(); ++it)
				{
					if((b1^(read[*it]&mask[j])).count()<bestmatch)
					{
						bestmatch = (b1^(read[*it]&mask[j])).count();
						bestmatchpos = *it;
					}
				}
				if(bestmatch <= thresh)
				{
					current = bestmatchpos;
					flag = 1;
					revflag = 0;
					revcomp.push_back(revflag);
					sortedorder.push_back(current);
					break;
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
				int bestmatch = 2*readlen,bestmatchpos;
				for (std::set<int>::iterator it = s1.begin() ; it != s1.end(); ++it)
				{
					if((b2^(read[*it]&revmask[j])).count()<bestmatch)
					{
						bestmatch = (b2^(read[*it]&revmask[j])).count();
						bestmatchpos = *it;
					}
				}
				if(bestmatch <= thresh)
				{
					current = bestmatchpos;
					flag = 1;
					revflag = 1;
					revcomp.push_back(revflag);
					sortedorder.push_back(current);
					break;
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



void writetofile(std::bitset<2*readlen> *read,std::bitset<2*readlen> *revread, std::vector<int> &sortedorder,std::vector<int> &revcomp)
{
	std::ofstream fout(outfile,std::ofstream::out);
	std::ofstream foutRC(outfileRC,std::ofstream::out);
	std::vector<int>::iterator it1,it2;
	for (it1 = sortedorder.begin(),it2 = revcomp.begin() ; it1 != sortedorder.end(); ++it1,++it2)
	{
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
