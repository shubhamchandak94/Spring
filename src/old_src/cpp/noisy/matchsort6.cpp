//Reordering for real reads
//Similar to matchsort3 but has a recoversingletons function and a second thresh "thresh1".
//After reordering is done, the recoversingletons function is called to try and place the singleton read before a matching read.
//However we try to make sure that the already existing order of the reads is not disrupted too much. 
//For example if we had read r1 and then r2 shifted by 4. Let r3 was a singleton which matched before r2 but with offset of 5.
//This is bad because r3 is shifted in the wrong direction with respect to r1.
//Thus we store the shifts as auxiliary information with the ordering. Since we want to modify the ordering later, 
//the sortedorder is no longer a simple vector - it is more like a doubly linked list, and with each read we store the previous
//and the next read, whether it was matched, if yes then whether it was reverse complemented while matching and the shift to its
//previous read.

//This does reduce the singleton reads significantly, however the effect on the overall size is small because the noise files 
//become larger and offset most of the gains in the seq file.

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
#define outfileflag "tempflag2.txt"
#define readlen 100
#define maxmatch 20
#define numreads 67155743
#define thresh 16
#define thresh1 24 //thresh for 2nd stage
#define numdict 2

void generateindexmasks(std::bitset<2*readlen> *mask1)
//function to generate dictionary index masks
//should be symmetric about readlen (e.g. for 2 dicts - if first dict is start1-end1 (both included), 
//then second should be (readlen-1-end1)-(readlen-1-start1))
{
	for(int i = 0; i < numdict; i++)
		mask1[i].reset();
	for(int i = 2*30; i < 2*50; i++)
		mask1[0][i] = 1;
	for(int i = 2*50; i < 2*70; i++)
		mask1[1][i] = 1;
	return;
}

void stringtobitset(std::string s,std::bitset<2*readlen> &read, std::bitset<2*readlen> &revread);

std::string bitsettostring(std::bitset<2*readlen> b);

void readDnaFile(std::bitset<2*readlen> *read,std::bitset<2*readlen> *revread);

void constructdictionary(std::bitset<2*readlen> *read, std::unordered_map<std::bitset<2*readlen>,std::vector<int>> *dict);

void generateindexmasks(std::bitset<2*readlen> *mask1);

void generatemasks(std::bitset<2*readlen> *mask,std::bitset<2*readlen> *revmask);

void reorder(std::bitset<2*readlen> *read, std::bitset<2*readlen> *revread, std::unordered_map<std::bitset<2*readlen>,std::vector<int>> *dict,std::unordered_map<int,std::vector<int>> &sortedorder);

void recoversingletons(std::bitset<2*readlen> *read, std::bitset<2*readlen> *revread, std::unordered_map<std::bitset<2*readlen>,std::vector<int>> *dict,std::unordered_map<int,std::vector<int>> &sortedorder);

void writetofile(std::bitset<2*readlen> *read,std::bitset<2*readlen> *revread, std::unordered_map<int,std::vector<int>> &sortedorder);

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
	std::unordered_map<int,std::vector<int>> sortedorder;//readno - {prev,next,matched?,revflag,offset}
	//0th read has prev -1 and last read has next 0
	std::cout << "Reordering reads\n";
	reorder(read,revread,dict,sortedorder);
	std::cout << "Trying to rematch singletons\n";
	recoversingletons(read,revread,dict,sortedorder);
	std::cout << "Writing to file\n";
	writetofile(read,revread,sortedorder);	
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

void reorder(std::bitset<2*readlen> *read, std::bitset<2*readlen> *revread, std::unordered_map<std::bitset<2*readlen>,std::vector<int>> *dict, std::unordered_map<int,std::vector<int>> &sortedorder)
{	
	std::bitset<2*readlen> *mask = new std::bitset<2*readlen> [maxmatch];
	std::bitset<2*readlen> *revmask = new std::bitset<2*readlen> [maxmatch];
	generatemasks(mask,revmask);
	std::bitset<2*readlen> *mask1 = new std::bitset<2*readlen> [numdict];
	generateindexmasks(mask1);
	std::unordered_set<int> remainingreads;

	for(int i = 0; i < numreads; i++)
		remainingreads.insert(i);
	int unmatched = 1;
	int current = 0;
	int flag,revflag = 0;
	//flag to check if match was found or not, revflag to check if the current read is forward or reverse 
	sortedorder[0] = {-1,0,0,0,0};
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
						flag = 1;
						revflag = 0;
						sortedorder[k]={current,0,1,0,j};
						sortedorder[current][1] = k;
						current = k;
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
						flag = 1;
						revflag = 1;
						sortedorder[k]={current,0,1,1,j};
						sortedorder[current][1] = k;
						current = k;
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
			int k = *remainingreads.begin();
			revflag = 0;
		        sortedorder[k] = {current,0,0,0,0};
			sortedorder[current][1] = k;	
			current = k;
		}
	}
	std::cout << "Reordering done, "<<unmatched<<" were unmatched\n";
	return;
}

void recoversingletons(std::bitset<2*readlen> *read, std::bitset<2*readlen> *revread, std::unordered_map<std::bitset<2*readlen>,std::vector<int>> *dict, std::unordered_map<int,std::vector<int>> &sortedorder)
{	
	std::bitset<2*readlen> *mask = new std::bitset<2*readlen> [maxmatch];
	std::bitset<2*readlen> *revmask = new std::bitset<2*readlen> [maxmatch];
	generatemasks(mask,revmask);
	std::bitset<2*readlen> *mask1 = new std::bitset<2*readlen> [numdict];
	generateindexmasks(mask1);
	int matched = 0;
	int flag;
	//flag to check if match was found or not, revflag to check if the current read is forward or reverse 
	std::bitset<2*readlen> *b = new std::bitset<2*readlen> [numdict];
	std::bitset<2*readlen> b1,b2;

	std::cout << "Constructing dictionaries again\n";
	for(int i = 0; i < numdict; i++)
		dict[i].clear();
	constructdictionary(read,dict);
	
	int current = 0;
	int next,prev;
	std::vector<int> currentvec;
	while(1)
	{
		currentvec = sortedorder[current];
		next = currentvec[1];
		prev = currentvec[0];
		if(prev == -1)//first read
		{
			current = next;
			continue;
		}	
		if(next == 0)//last read
			break;
		if(currentvec[2] == 0 && sortedorder[currentvec[1]][2] == 0)//singleton read
		{
			b1 = read[current];
			b2 = revread[current];
			flag = 0;
			for(int j = 0; j < maxmatch; j++)
			{
				std::set<int> s;
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
						if(k==current||(sortedorder[k][2] == 1 && sortedorder[k][4]<j)||sortedorder[k][0]==-1)
							continue;//if read k is matched with a smaller offset or is the first read, skip
						if((b1^(read[k]&mask[j])).count()<=thresh1)
						{
							flag = 1;
							if(sortedorder[k][2] == 1)//k was matched to something
							{
								sortedorder[prev][1] = next;
								sortedorder[next][0] = prev;
								sortedorder[sortedorder[k][0]][1] = current;
								sortedorder[current][0] = sortedorder[k][0];
								sortedorder[current][1] = k;
								sortedorder[current][2] = 1;
								sortedorder[current][3] = sortedorder[k][3];
								sortedorder[current][4] = sortedorder[k][4] -j;	
								sortedorder[k][0] = current;
								sortedorder[k][4] = j;
							}
							else
							{
								sortedorder[prev][1] = next;
								sortedorder[next][0] = prev;
								sortedorder[sortedorder[k][0]][1] = current;
								sortedorder[current][0] = sortedorder[k][0];
								sortedorder[current][1] = k;
								sortedorder[k][0] = current;
								sortedorder[k][2] = 1;
								sortedorder[k][3] = 0;
								sortedorder[k][4] = j;
							}
							matched += 1;
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
						if(k==current||(sortedorder[k][2] == 1 && sortedorder[k][4]<j)||sortedorder[k][0]==-1)
							continue;//if read k is matched with a smaller offset or is the first read, skip
						if((b2^(read[k]&revmask[j])).count()<=thresh1)
						{
							flag = 1;
							if(sortedorder[k][2] == 1)//k was matched to something
							{
								sortedorder[prev][1] = next;
								sortedorder[next][0] = prev;
								sortedorder[sortedorder[k][0]][1] = current;
								sortedorder[current][0] = sortedorder[k][0];
								sortedorder[current][1] = k;
								sortedorder[current][2] = 1;
								sortedorder[current][3] = 1 - sortedorder[k][3];
								sortedorder[current][4] = sortedorder[k][4] -j;	
								sortedorder[k][0] = current;
								sortedorder[k][4] = j;
							}
							else
							{
								sortedorder[prev][1] = next;
								sortedorder[next][0] = prev;
								sortedorder[sortedorder[k][0]][1] = current;
								sortedorder[current][0] = sortedorder[k][0];
								sortedorder[current][1] = k;
								sortedorder[k][0] = current;
								sortedorder[k][2] = 1;
								sortedorder[k][3] = 1;
								sortedorder[k][4] = j;
							}
							matched += 1;
							break;
						}
					}
					if(flag == 1)
						break;
				}
				b2<<=2;
			}	
					
		}
		current = next;
	}
	std::cout << matched <<" singleton reads were unmatched\n";
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


void writetofile(std::bitset<2*readlen> *read,std::bitset<2*readlen> *revread, std::unordered_map<int,std::vector<int>> &sortedorder)
{
	std::ofstream fout(outfile,std::ofstream::out);
	std::ofstream foutRC(outfileRC,std::ofstream::out);
	std::ofstream foutflag(outfileflag,std::ofstream::out);
	int current = 0;
	std::vector<int> currentvec;
	while(1)
	{	
		currentvec = sortedorder[current];
		foutflag << currentvec[2];
		if(currentvec[3] == 0)
		{
			fout<<bitsettostring(read[current])<<"\n";
			foutRC << 'd';
		}
		else
		{
			fout<<bitsettostring(revread[current])<<"\n";
			foutRC << 'r';
		}
		if (currentvec[1] == 0)//last read
			break;
		else
			current = currentvec[1];
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
