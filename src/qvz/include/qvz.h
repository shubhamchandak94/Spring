#ifndef SPRING_QVZ_QVZ_H_
#define SPRING_QVZ_QVZ_H_

#include <string>

#include "qvz/include/codebook.h"

namespace spring {
namespace qvz {

#define ALPHABET_SIZE 72
/**
 *
 */
void encode(struct qv_options_t *opts, uint32_t max_readlen, uint32_t numreads, std::string *quality_string_array, uint32_t *str_len_array);
/**
 *
 */

} // namespace qvz
} // namespace spring

#endif // SPRING_QVZ_QVZ_H_

