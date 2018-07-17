#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "distortion.h"
#include "util.h"

/**
 * Allocates memory for a distortion matrix
 */
struct distortion_t *alloc_distortion_matrix(uint8_t symbols) {
	struct distortion_t *rtn = (struct distortion_t *) calloc(1, sizeof(struct distortion_t));
	rtn->symbols = symbols;
	rtn->distortion = (double *) calloc(symbols*symbols, sizeof(double));
	return rtn;
}

/**
 * Deallocates memory from a distortion matrix
 */
void free_distortion_matrix(struct distortion_t *d) {
	free(d->distortion);
	free(d);
}

/**
 * Public facing method for allocating distortion matrices
 */
struct distortion_t *generate_distortion_matrix(uint8_t symbols, int type) {
	switch (type) {
		case DISTORTION_MANHATTAN:
			return gen_manhattan_distortion(symbols);
		case DISTORTION_MSE:
			return gen_mse_distortion(symbols);
		case DISTORTION_LORENTZ:
			return gen_lorentzian_distortion(symbols);
		case DISTORTION_CUSTOM:
			printf("Custom distortion matrices should be allocated with gen_custom_distortion() instead.\n");
			exit(1);
		default:
			printf("Invalid distortion type %d specified.\n", type);
			exit(1);
	}
}

/**
 * Generate a distortion matrix according to the Manhattan distance (L1) metric
 */
struct distortion_t *gen_manhattan_distortion(uint8_t symbols) {
	struct distortion_t *rtn = alloc_distortion_matrix(symbols);
	uint8_t x, y;

	for (x = 0; x < symbols; ++x) {
		for (y = 0; y < symbols; ++y) {
			rtn->distortion[x + y*symbols] = abs(x - y);
		}
	}

	return rtn;
}

/**
 * Generates a distortion matrix according to the MSE (L2) metric
 */
struct distortion_t *gen_mse_distortion(uint8_t symbols) {
	struct distortion_t *rtn = alloc_distortion_matrix(symbols);
	uint8_t x, y;

	for (x = 0; x < symbols; ++x) {
		for (y = 0; y < symbols; ++y) {
			rtn->distortion[x + y*symbols] = (x - y)*(x - y);
		}
	}

	return rtn;
}

/**
 * Generates a distortion matrix according to the lorentzian (log-L1) metric
 */
struct distortion_t *gen_lorentzian_distortion(uint8_t symbols) {
	struct distortion_t *rtn = alloc_distortion_matrix(symbols);
	uint8_t x, y;

	for (x = 0; x < symbols; ++x) {
		for (y = 0; y < symbols; ++y) {
			rtn->distortion[x + y*symbols] = log2( 1.0 + (double)(abs(x - y)) );
		}
	}

	return rtn;
}

/**
 * Reads in a custom distortion matrix specified in the given file
 * The file format is S rows of S columns containing double valued distortions
 * separated by commas. Lines beginning with a # are ignored as comments
 */
struct distortion_t *gen_custom_distortion(uint8_t symbols, const char *filename) {
	struct distortion_t *dist = alloc_distortion_matrix(symbols);
	uint8_t x, y;
	FILE *fp;
	char line[1024];
	char *field;
	uint8_t missing;

	fp = fopen(filename, "rt");
	if (!fp) {
		perror("Unable to open distortion definition file");
		exit(1);
	}

	x = 0;
	while (x < symbols && fgets(line, 1024, fp) != NULL) {
		missing = 0;
		field = line - 1;
		y = 0;

		if (line[0] == '#')
			continue;

		while (y < symbols && field != NULL) {
			field += 1;
			dist->distortion[x + symbols*y] = atof(field);
			field = strchr(field, ',');
			y += 1;
		}

		while (y < symbols) {
			missing = 1;
			dist->distortion[x + symbols*y] = 0.0;
		}

		if (missing) {
			printf("Warning: one or more entries in the distortion matrix on line %d were missing", x);
			printf(" they have been filled with 0.0\n");
		}

		x += 1;
	}

	fclose(fp);
	return dist;
}

/**
 * Retrieve the distortion for a pair (x, y). Generally x is the true value and
 * y is the reconstructed value. Handles the matrix->linear array indexing
 */
double get_distortion(struct distortion_t *dist, uint8_t x, uint8_t y) {
	return dist->distortion[x + dist->symbols*y];
}

/**
 * Print a distortion matrix to stdout for debuggin
 */
void print_distortion(struct distortion_t *dist) {
	uint8_t x, y;

	printf("    |");
	for (y = 0; y < dist->symbols; ++y) {
		printf(" %2d |", y);
	}
	printf("\n");

	printf("----+");
	for (y = 0; y < dist->symbols; ++y) {
		printf("----+");
	}
	printf("\n");

	for (x = 0; x < dist->symbols; ++x) {
		printf(" %2d |", x);
		for (y = 0; y < dist->symbols; ++y) {
			printf("%2.2f|", dist->distortion[x + y*dist->symbols]);
		}
		printf("\n");
	}
}
