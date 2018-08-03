#ifndef SPRING_QVZ_MAIN_H_
#define SPRING_QVZ_MAIN_H_

#include "qvz/include/util.h"

#include <stdio.h>
#include <string>

#include "qvz/include/cluster.h"
#include "qvz/include/codebook.h"
#include "qvz/include/qv_compressor.h"

namespace spring {
namespace qvz {

#define ALPHABET_SIZE 72
/**
 *
 */
void encode(FILE *fout, struct qv_options_t *opts, uint32_t max_readlen,
            uint32_t numreads, char *quality_array, uint16_t *read_lengths,
            std::string &infile_order, uint64_t startpos);

void encode_lossless(FILE *fout, struct qv_options_t *opts, uint32_t max_readlen, uint32_t numreads, std::string *quality_string_array);
/**
 *
 */
void decode(char *input_file, char *output_file, struct qv_options_t *opts,
            uint16_t *read_lengths);

/**
 * Displays a usage name
 * @param name Program name string
 */
void usage(char *name);

} // namespace qvz
} // namespace spring

#endif // SPRING_QVZ_MAIN_H_
