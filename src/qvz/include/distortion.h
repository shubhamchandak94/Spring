#ifndef _DISTORTION_H_
#define _DISTORTION_H_

#include <stdint.h>

// Constants for generating distortion matrices
#define DISTORTION_MANHATTAN		1
#define DISTORTION_MSE				2
#define DISTORTION_LORENTZ			3
#define DISTORTION_CUSTOM			4

/**
 * Used to store distortion matrix information that is used during quantizer generation
 */
struct distortion_t {
	double *distortion;
	uint8_t symbols;
};

// Memory management functions
struct distortion_t *alloc_distortion_matrix(uint8_t symbols);
void free_distortion_matrix(struct distortion_t *);

// Methods for generating distortion matrices of different types
struct distortion_t *generate_distortion_matrix(uint8_t symbols, int type);
struct distortion_t *gen_mse_distortion(uint8_t symbols);
struct distortion_t *gen_manhattan_distortion(uint8_t symbols);
struct distortion_t *gen_lorentzian_distortion(uint8_t symbols);
struct distortion_t *gen_custom_distortion(uint8_t symbols, const char *filename);

// Accessors
double get_distortion(struct distortion_t *dist, uint8_t x, uint8_t y);

void print_distortion(struct distortion_t *dist);

#endif
