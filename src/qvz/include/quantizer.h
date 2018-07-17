#ifndef SPRING_QVZ_QUANTIZER_H_
#define SPRING_QVZ_QUANTIZER_H_

#include <stdint.h>

#include "algorithms/SPRING/qvz/include/distortion.h"
#include "algorithms/SPRING/qvz/include/pmf.h"
#include "algorithms/SPRING/qvz/include/util.h"

#define QUANTIZER_MAX_ITER 100

namespace spring {
namespace qvz {

/**
 * Structure holding information about a quantizer, which just maps input
 * symbols
 * to output symbols for a specific alphabet
 */
struct quantizer_t {
  const struct alphabet_t *restrict alphabet;
  struct alphabet_t *restrict output_alphabet;
  symbol_t *restrict q;
  double ratio;
  double mse;
};

// Memory management
struct quantizer_t *alloc_quantizer(const struct alphabet_t *);
void free_quantizer(struct quantizer_t *);

// Generates a quantizer via optimization
struct quantizer_t *generate_quantizer(struct pmf_t *restrict pmf,
                                       struct distortion_t *restrict dist,
                                       uint32_t states);

// Calculate the output pmf when the quantizer is applied to the input pmf
struct pmf_t *apply_quantizer(struct quantizer_t *restrict q,
                              struct pmf_t *restrict pmf,
                              struct pmf_t *restrict output);

// Find the output alphabet of a quantizer
void find_output_alphabet(struct quantizer_t *);

// Display/debugging
void print_quantizer(struct quantizer_t *);

} // namespace qvz
} // naemspace spring

#endif
