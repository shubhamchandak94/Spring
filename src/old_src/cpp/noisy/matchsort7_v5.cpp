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
#define outfile "temp0.dna"
#define outfileRC "tempRC0.txt"
#define outfileflag "tempflag0.txt"
#define readlen 100
#define maxmatch 30
#define numreads 67155743
#define thresh 24
#define numdict 2


void generateindexmasks(std::bitset<2*readlen> *mask1)
{
	for(int i = 0; i < numdict; i++)
		mask1[i].reset();
	for(int i = 2*30; i < 2*50; i++)
		mask1[0][i] = 1;
	for(int i = 2*50; i < 2*70; i++)
		mask1[1][i] = 1;
	return;
}

char inttochar[] = {'A','C','G','T'};

std::bitset<2*readlen> stringtobitset(std::string s);

std::string bitsettostring(std::bitset<2*readlen> b);

void readDnaFile(std::bitset<2*readlen> *read);

void constructdictionary(std::bitset<2*readlen> *read, std::unordered_map<std::bitset<2*readlen>,std::vector<int>> *dict);

void generatemasks(std::bitset<2*readlen> *mask,std::bitset<2*readlen> *revmask);

void reorder(std::bitset<2*readlen> *read, std::unordered_map<std::bitset<2*readlen>,std::vector<int>> *dict,std::vector<int> &sortedorder,std::vector<int> &revcomp,std::vector<int>& flagvec);

void writetofile(std::bitset<2*readlen> *read, std::vector<int> &sortedorder,std::vector<int> &revcomp,std::vector<int> &flagvec);

void updaterefcount(std::bitset<2*readlen> cur, std::bitset<2*readlen> &ref, std::bitset<2*readlen> &revref, int count[][readlen], bool resetcount, bool rev, int shift);

std::bitset<2*readlen> reverse_complement(std::bitset<2*readlen> b);

int main()
{
	std::bitset<2*readlen> *read = new std::bitset<2*readlen> [numreads];
	std::cout << "Reading file\n";
	readDnaFile(read);
	std::cout << "Constructing dictionaries\n";
	//using vector instead of list to save some space (at the cost of linear time required to delete elements)
	std::unordered_map<std::bitset<2*readlen>,std::vector<int>> *dict = new std::unordered_map<std::bitset<2*readlen>,std::vector<int>> [numdict];
	constructdictionary(read,dict);
	std::vector<int> sortedorder,revcomp,flagvec;
	std::cout << "Reordering reads\n";
	reorder(read,dict,sortedorder,revcomp,flagvec);
	std::cout << "Writing to file\n";
	writetofile(read,sortedorder,revcomp,flagvec);	
	std::cout << "Done!\n";
	return 0;
}


std::bitset<2*readlen> stringtobitset(std::string s)
{
	int i;
	std::bitset<2*readlen> b;
	for(int i = 0; i < readlen; i++)
	{	
		switch(s[i])
		{
			case 'A':	break;
			case 'C':	b[2*i+1] = 1;
					break;
			case 'G':	b[2*i] = 1;
					break;
			case 'T':	b[2*i] = 1;
					b[2*i+1] = 1;
					break;
		}
	}
	return b;
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

void reorder(std::bitset<2*readlen> *read, std::unordered_map<std::bitset<2*readlen>,std::vector<int>> *dict, std::vector<int> &sortedorder,std::vector<int> &revcomp,std::vector<int> &flagvec)
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



void writetofile(std::bitset<2*readlen> *read, std::vector<int> &sortedorder,std::vector<int> &revcomp,std::vector<int> &flagvec)
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
			fout<<bitsettostring(reverse_complement(read[*it1]))<<"\n";
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
		s.push_back(inttochar[2*b[2*i]+b[2*i+1]]);
	return s;
}

std::bitset<2*readlen> reverse_complement(std::bitset<2*readlen> b)
{
	std::bitset<2*readlen> b1;
	for(int j = 0; j < readlen; j++)
	{
		b1[2*j] = 1- b[2*(readlen-j-1)];
		b1[2*j+1] = 1 - b[2*(readlen-j-1) + 1];
	}
	return b1;
}

void updaterefcount(std::bitset<2*readlen> cur, std::bitset<2*readlen> &ref, std::bitset<2*readlen> &revref, int count[][readlen], bool resetcount, bool rev, int shift)
{
	std::bitset<2*readlen> current;
	if(rev == false)
		current = cur;
	else
		current = reverse_complement(cur);

	if(resetcount == true)
	{
		ref = current;	
		for(int i = 0; i < readlen; i++)
		{	
			for(int j = 0; j < 4; j++)
				count[j][i] = 0;
			count[2*current[2*i]+current[2*i+1]][i] = 1;
		}

	}
	else
	{
		for(int i = 0; i < readlen-shift; i++)
		{	
			for(int j = 0; j < 4; j++)
				count[j][i] = count[j][i+shift];
			count[2*current[2*i]+current[2*i+1]][i] += 1;
			
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
			count[2*current[2*i]+current[2*i+1]][i] = 1;
		}
	}

	revref = reverse_complement(ref);		
	return;
}
