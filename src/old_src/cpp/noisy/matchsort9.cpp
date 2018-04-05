#include <iostream>
#include <fstream>
#include <bitset>
#include <unordered_set>
#include <unordered_map>
#include <string>
#include <vector>
#include <algorithm>
#include <set>

#define infile "SRR065390.dna"
#define outfile "temp0.dna"
#define outfileRC "tempRC0.txt"
#define outfileflag "tempflag0.txt"
#define readlen 100
#define maxmatch 20
#define numreads 67617092
#define thresh 16
#define numdict 2

void stringtobitset(std::string s,std::bitset<4*readlen> &read, std::bitset<4*readlen> &revread);


std::string bitsettostring(std::bitset<4*readlen> b);

void readDnaFile(std::bitset<4*readlen> *read,std::bitset<4*readlen> *revread);

void constructdictionary(std::bitset<4*readlen> *read, std::unordered_map<std::bitset<4*readlen>,std::vector<int>> *dict);

void generateindexmasks(std::bitset<4*readlen> *mask1);

void generatemasks(std::bitset<4*readlen> *mask,std::bitset<4*readlen> *revmask);

void reorder(std::bitset<4*readlen> *read, std::bitset<4*readlen> *revread, std::unordered_map<std::bitset<4*readlen>,std::vector<int>> *dict,std::vector<int> &sortedorder,std::vector<int> &revcomp,std::vector<int>& flagvec);

void writetofile(std::bitset<4*readlen> *read,std::bitset<4*readlen> *revread, std::vector<int> &sortedorder,std::vector<int> &revcomp,std::vector<int> &flagvec);

void updaterefcount(std::bitset<4*readlen> current, std::bitset<4*readlen> &ref, std::bitset<4*readlen> &revref, int count[][readlen], bool resetcount, int shift);


int main()
{
	std::bitset<4*readlen> *read = new std::bitset<4*readlen> [numreads];
	std::bitset<4*readlen> *revread = new std::bitset<4*readlen> [numreads];
	std::cout << "Reading file\n";
	//using vector instead of list to save some space (at the cost of linear time required to delete elements)
	readDnaFile(read,revread);
	std::cout << "Constructing dictionaries\n";
	std::unordered_map<std::bitset<4*readlen>,std::vector<int>> *dict = new std::unordered_map<std::bitset<4*readlen>,std::vector<int>> [numdict];
	constructdictionary(read,dict);
	std::vector<int> sortedorder,revcomp,flagvec;
	std::cout << "Reordering reads\n";
	reorder(read,revread,dict,sortedorder,revcomp,flagvec);
	std::cout << "Writing to file\n";
	writetofile(read,revread,sortedorder,revcomp,flagvec);	
	std::cout << "Done!\n";
	return 0;
}


void stringtobitset(std::string s,std::bitset<4*readlen> &b, std::bitset<4*readlen> &b1)
{
	b.reset();
	b1.reset();
	for(int i = 0; i < readlen; i++)
	{	
		switch(s[i])
		{
			case 'A':	b[4*i] = 1;
					b1[4*(readlen-i-1)+3] = 1;
					break;
			case 'C':	b[4*i+1] = 1;
					b1[4*(readlen-i-1)+2] = 1;
					break;
			case 'G':	b[4*i+2] = 1;
					b1[4*(readlen-i-1)+1] = 1;
					break;
			case 'T':	b[4*i+3] = 1;
					b1[4*(readlen-i-1)] = 1;
					break;
		}
	}
	return;
}

void readDnaFile(std::bitset<4*readlen> *read, std::bitset<4*readlen> *revread)
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
		

void constructdictionary(std::bitset<4*readlen> *read, std::unordered_map<std::bitset<4*readlen>,std::vector<int>> *dict)
{
	std::bitset<4*readlen> b;
	std::bitset<4*readlen> *mask = new std::bitset<4*readlen> [numdict];
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

void reorder(std::bitset<4*readlen> *read, std::bitset<4*readlen> *revread, std::unordered_map<std::bitset<4*readlen>,std::vector<int>> *dict, std::vector<int> &sortedorder,std::vector<int> &revcomp,std::vector<int> &flagvec)
{	
	std::bitset<4*readlen> *mask = new std::bitset<4*readlen> [maxmatch];
	std::bitset<4*readlen> *revmask = new std::bitset<4*readlen> [maxmatch];
	std::bitset<4*readlen> ref,revref;
	int count[4][readlen];
	generatemasks(mask,revmask);
	std::bitset<4*readlen> *mask1 = new std::bitset<4*readlen> [numdict];
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
	std::bitset<4*readlen> *b = new std::bitset<4*readlen> [numdict];
	std::bitset<4*readlen> b1,b2;
	
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
			b1 >>= 4;

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
			b2 <<= 4;
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


void generatemasks(std::bitset<4*readlen> *mask,std::bitset<4*readlen> *revmask)
{
	for(int i = 0; i < maxmatch; i++)
	{	
		mask[i].reset();
		revmask[i].reset();
		for(int j = 0; j < 4*readlen - 4*i; j++)
			mask[i][j] = 1;
		for(int j = 4*i; j < 4*readlen; j++)
			revmask[i][j] = 1; 	
	}
	return;
}


void generateindexmasks(std::bitset<4*readlen> *mask1)
{
	for(int i = 0; i < numdict; i++)
		mask1[i].reset();
	for(int i = 4*30; i < 4*50; i++)
		mask1[0][i] = 1;
	for(int i = 4*50; i < 4*70; i++)
		mask1[1][i] = 1;
	return;
}

void writetofile(std::bitset<4*readlen> *read,std::bitset<4*readlen> *revread, std::vector<int> &sortedorder,std::vector<int> &revcomp,std::vector<int> &flagvec)
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

std::string bitsettostring(std::bitset<4*readlen> b)
{
	std::string s;
	for(int i = 0; i < readlen; i++)
	{
		if (b[4*i] == 1)
			s.push_back('A');
		else if (b[4*i+1] == 1)
			s.push_back('C');
		else if (b[4*i+2] == 1)
			s.push_back('G');
		else if (b[4*i+3] == 1)
			s.push_back('T');
		else
			s.push_back('N');
	}
	return s;
}


void updaterefcount(std::bitset<4*readlen> current, std::bitset<4*readlen> &ref, std::bitset<4*readlen> &revref, int count[][readlen], bool resetcount, int shift)
{
	if(resetcount == true)
	{
		ref = current;	
		for(int i = 0; i < readlen; i++)
		{	
			for(int j = 0; j < 4; j++)
				count[j][i] = 0;
			if (current[4*i] == 1)
				count[0][i] = 1;
			else if (current[4*i+1] == 1)
				count[1][i] = 1;
			else if (current[4*i+2] == 1)
				count[2][i] = 1;
			else if (current[4*i+3] == 1)
				count[3][i] = 1;
		}

	}
	else
	{
		ref.reset();
		for(int i = 0; i < readlen-shift; i++)
		{	
			for(int j = 0; j < 4; j++)
				count[j][i] = count[j][i+shift];
			if (current[4*i] == 1)
				count[0][i] += 1;
			else if (current[4*i+1] == 1)
				count[1][i] += 1;
			else if (current[4*i+2] == 1)
				count[2][i] += 1;
			else if (current[4*i+3] == 1)
				count[3][i] += 1;

			int max = 0,indmax = 0;
			for(int j = 0; j < 4; j++)
				if(count[j][i]>max)
				{
					max = count[j][i];
					indmax = j;
				}
			switch(indmax)
			{
				case 0: ref[4*i] = 1;
					break;			
				case 1: ref[4*i+1] = 1;
					break;			
				case 2: ref[4*i+2] = 1;
					break;			
				case 3: ref[4*i+3] = 1;
					break;			
			}
		}
		
		for(int i = readlen-shift; i < readlen; i++)
		{	
			ref[4*i] = current[4*i];
			ref[4*i+1] = current[4*i+1];
			ref[4*i+2] = current[4*i+2];
			ref[4*i+3] = current[4*i+3];
			for(int j = 0; j < 4; j++)
				count[j][i] = 0;

			if (current[4*i] == 1)
				count[0][i] = 1;
			else if (current[4*i+1] == 1)
				count[1][i] = 1;
			else if (current[4*i+2] == 1)
				count[2][i] = 1;
			else if (current[4*i+3] == 1)
				count[3][i] = 1;
		}
	}
	revref.reset();	
	for(int j = 0; j < readlen; j++)
	{	
		if(ref[4*j] == 1)
			revref[4*(readlen-j-1)+3] = 1;
		else if(ref[4*j+1] == 1)
			revref[4*(readlen-j-1)+2] = 1;
		else if(ref[4*j+2] == 1)
			revref[4*(readlen-j-1)+1] = 1;
		else if(ref[4*j+3] == 1)
			revref[4*(readlen-j-1)] = 1;
	}
	return;
}

