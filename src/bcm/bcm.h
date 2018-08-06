#ifndef SPRING_BCM_H_
#define SPRING_BCM_H_

#include <string>
#include "params.h"

namespace spring {
namespace bcm {

int bcm_compress(const char *infile, const char *outfile, int bsize = BCM_BLOCK_SIZE);

int bcm_decompress(const char *infile, const char *outfile);

int bcm_str_array_compress(const char *outfile, std::string *str_array_param, uint32_t size_str_array_param, uint32_t *str_lengths_param, int bsize = BCM_BLOCK_SIZE);

int bcm_str_array_decompress(const char *infile, std::string *str_array_param, uint32_t size_str_array_param, uint32_t *str_lengths_param);

} // namespace bcm
} // namespace spring

#endif // SPRING_BCM_H_
