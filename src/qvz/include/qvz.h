#ifndef SPRING_QVZ_QVZ_H_
#define SPRING_QVZ_QVZ_H_

#include <stdio.h>
#include <string>

#include "qvz/include/codebook.h"

namespace spring {
namespace qvz {

#define ALPHABET_SIZE 72
/**
 *
 */
void encode_lossless(const char * outfile_name, struct qv_options_t *opts, uint32_t max_readlen, uint32_t numreads, std::string *quality_string_array);

void decode_lossless(const char * infile_name, struct qv_options_t *opts, uint32_t max_readlen, uint32_t numreads, std::string *quality_string_array, uint16_t *read_lengths);
/**
 *
 */

} // namespace qvz
} // namespace spring

#endif // SPRING_QVZ_MAIN_H_
