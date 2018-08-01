#include <iostream>
#include <fstream>
#include <string>
#define readlen 101

std::string infile;
std::string outfile;

void computecounts();

#define num_clusters 4
int clusterboundary[] = {10,20,30};

int main(int argc, char** argv)
{
	infile = std::string(argv[1]);
	outfile = std::string(argv[2]);

	computecounts();	
	return 0;
}

void computecounts()
{
	std::string line;
	std::ifstream myfile(infile, std::ifstream::in);
	uint64_t (*count1)[readlen][42] = new uint64_t [num_clusters][readlen][42]();
	uint64_t (*count2)[readlen-1][42][42] = new uint64_t [num_clusters][readlen-1][42][42]();
	uint64_t (*count3)[readlen-2][42][42][42] = new uint64_t [num_clusters][readlen-2][42][42][42]();
	uint64_t num_reads[num_clusters] = {};	
	while(std::getline(myfile, line))
	{
		//compute average quality value
		int total_qv = 0;
		for(int i = 0; i < readlen; i++)
			total_qv += ((unsigned int)(line[i])-33);
		double avg_qv = (double)(total_qv)/readlen;
		//find cluster number
		int cluster_num = -1;
		for(int j = 0; j < num_clusters-1; j++)
			if(avg_qv < clusterboundary[j])
				cluster_num = j;

		if(cluster_num == -1)
			cluster_num = num_clusters-1;
		num_reads[cluster_num]++;
		for(int i = 0; i < readlen-2; i++)
		{
			count1[cluster_num][i][(unsigned int)line[i]-33] += 1;
			count2[cluster_num][i][(unsigned int)line[i]-33][(unsigned int)line[i+1]-33] += 1;
			count3[cluster_num][i][(unsigned int)line[i]-33][(unsigned int)line[i+1]-33][(unsigned int)line[i+2]-33] += 1;
		}
		int i = readlen - 2;
		count1[cluster_num][i][(unsigned int)line[i]-33] += 1;
		count2[cluster_num][i][(unsigned int)line[i]-33][(unsigned int)line[i+1]-33] += 1;
		i = readlen - 1;
		count1[cluster_num][i][(unsigned int)line[i]-33] += 1;
	}	
	myfile.close();	
		
	std::ofstream f_out(outfile, std::ios::binary);
	for(int j = 0; j < num_clusters; j++)
	{
		f_out.write((char*)&num_reads[j],sizeof(uint64_t));
		for(int i = 0; i < readlen; i++)
			for(int k = 0; k < 42; k++)
				f_out.write((char*)&count1[j][i][k],sizeof(uint64_t));
		for(int i = 0; i < readlen-1; i++)
			for(int k = 0; k < 42; k++)
				for(int k1 = 0; k1 < 42; k1++)
					f_out.write((char*)&count2[j][i][k][k1],sizeof(uint64_t));
		for(int i = 0; i < readlen-2; i++)
			for(int k = 0; k < 42; k++)
				for(int k1 = 0; k1 < 42; k1++)
					for(int k2 = 0; k2 < 42; k2++)
					f_out.write((char*)&count3[j][i][k][k1][k2],sizeof(uint64_t));
	}
	f_out.close();	
	return;	
}
