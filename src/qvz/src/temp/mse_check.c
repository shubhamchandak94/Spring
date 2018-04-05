#include <stdlib.h>
#include <stdio.h>

int main( int argc, char **argv){
	
	FILE *f1, *f2;
	char *line1, *line2;
	int error = 0, i = 0;
	unsigned int num_lines = 0, columns = 36, num_lines_total = 500000;
	double distortion = 0.0;

	line1 = (char*)calloc(4096, 1);
	line2 = (char*)calloc(4096, 1);
	f1 = fopen(argv[1], "r");
	f2 = fopen(argv[2], "r");
	
	fgets(line1, columns+2, f1);
	fgets(line2, columns+2, f2);
	do{
		num_lines++;
		error = 0;	
		if(num_lines%1000000 == 0) printf("Line: %dM\n", num_lines/1000000);
		for(i = 0; i<columns;i++){
			error += ((int)line1[i] - (int)line2[i])*((int)line1[i] - (int)line2[i]);
		}
		distortion += error / ( (double)columns);
		fgets(line1, columns+2, f1);
		fgets(line2, columns+2, f2);
	}while (num_lines < num_lines_total);
	distortion = distortion / (double)num_lines;
	printf ("MSE:%f\n", distortion);

}
