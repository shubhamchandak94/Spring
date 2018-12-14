#ifndef SPRING_QVZ_WELL_H_
#define SPRING_QVZ_WELL_H_
/**
 * Declarations necessary for using the WELL-1024a PRNG in other code
 */

#include <stdint.h>

namespace spring {
namespace qvz {

struct well_state_t {
  uint32_t state[32];
  uint32_t n;
  uint32_t bit_output;
  uint32_t bits_left;
};

uint32_t well_1024a(struct well_state_t *state);
uint32_t well_1024a_bits(struct well_state_t *state, uint8_t bits);

} // namespace qvz
} // namespace spring

#endif
