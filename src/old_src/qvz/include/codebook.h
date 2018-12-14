#ifndef _CODEBOOK_H_
#define _CODEBOOK_H_
/**
 * Functions and definitions relating to reading codebooks from files, used
 * for both the encoder and decoder code
 */

#include "util.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "well.h"
#include "pmf.h"
#include "distortion.h"
#include "quantizer.h"
#include "lines.h"

#define MODE_RATIO		0	// Traditional implementation, output bitrate is scaled from input
#define MODE_FIXED		1	// Fixed rate per symbol
#define MODE_FIXED_MSE	2	// Fixed average MSE per column

/**
 * Options for the compression process
 */
struct qv_options_t {
	uint8_t verbose;
	uint8_t stats;
	uint8_t mode;
	uint8_t clusters;
    uint8_t uncompressed;
    uint8_t distortion;
	char *dist_file;
    char *uncompressed_name;
	double ratio;		// Used for parameter to all modes
	double e_dist;		// Expected distortion as calculated during optimization
	double cluster_threshold;
};

/**
 * Stores an array of conditional PMFs for the current column given the previous
 * column. PMF pointers are stored in a flat array so don't try to find the PMF you
 * want directly--use the accessor
 */
struct cond_pmf_list_t {
	uint32_t columns;
	const struct alphabet_t *alphabet;
	struct pmf_t **pmfs;
	struct pmf_list_t *marginal_pmfs;
};

/**
 * Stores an array of quantizer pointers for the column for all possible left context
 * values. Unused ones are left as null pointers. This is also stored as a flat array
 * so the accessor must be used to look up the correct quantizer
 * The dreaded triple pointer is used to store an array of (different length) arrays
 * of pointers to quantizers
 */
struct cond_quantizer_list_t {
	uint32_t columns;
	uint32_t lines;
	struct alphabet_t **input_alphabets;
	struct quantizer_t ***q;
	double **ratio;				// Raw ratio
	uint8_t **qratio;			// Quantized ratio
	struct qv_options_t *options;
};

// Memory management
struct cond_pmf_list_t *alloc_conditional_pmf_list(const struct alphabet_t *alphabet, uint32_t columns);
struct cond_quantizer_list_t *alloc_conditional_quantizer_list(uint32_t columns);
void free_conditional_pmf_list(struct cond_pmf_list_t *);
void free_cond_quantizer_list(struct cond_quantizer_list_t *);

// Per-column initializer for conditional quantizer list
void cond_quantizer_init_column(struct cond_quantizer_list_t *list, uint32_t column, const struct alphabet_t *input_union);

// Accessors
struct pmf_t *get_cond_pmf(struct cond_pmf_list_t *list, uint32_t column, symbol_t prev);
struct quantizer_t *get_cond_quantizer_indexed(struct cond_quantizer_list_t *list, uint32_t column, uint32_t index);
struct quantizer_t *get_cond_quantizer(struct cond_quantizer_list_t *list, uint32_t column, symbol_t prev);
void store_cond_quantizers(struct quantizer_t *restrict lo, struct quantizer_t *restrict hi, double ratio, struct cond_quantizer_list_t *list, uint32_t column, symbol_t prev);
void store_cond_quantizers_indexed(struct quantizer_t *restrict lo, struct quantizer_t *restrict hi, double ratio, struct cond_quantizer_list_t *list, uint32_t column, uint32_t index);
struct quantizer_t *choose_quantizer(struct cond_quantizer_list_t *list, struct well_state_t *well, uint32_t column, symbol_t prev, uint32_t *q_idx);
uint32_t find_state_encoding(struct quantizer_t *codebook, symbol_t value);

// Meat of the implementation
void calculate_statistics(struct quality_file_t *);
double optimize_for_entropy(struct pmf_t *pmf, struct distortion_t *dist, double target, struct quantizer_t **lo, struct quantizer_t **hi);
void generate_codebooks(struct quality_file_t *info);

// Master functions to handle codebooks in the output file
void write_codebooks(FILE *fp, struct quality_file_t *info);
void write_codebook(FILE *fp, struct cond_quantizer_list_t *quantizers);
void read_codebooks(FILE *fp, struct quality_file_t *info);
struct cond_quantizer_list_t *read_codebook(FILE *fp, struct quality_file_t *info);

#define MAX_CODEBOOK_LINE_LENGTH 3366
#define COPY_Q_TO_LINE(line, q, i, size) for (i = 0; i < size; ++i) { line[i] = q[i] + 33; }
#define COPY_Q_FROM_LINE(line, q, i, size) for (i = 0; i < size; ++i) { q[i] = line[i] - 33; }

void print_codebook(struct cond_quantizer_list_t *);

#endif
