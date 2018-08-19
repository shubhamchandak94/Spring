#ifndef SPRING_DECOMPRESS_H_
#define SPRING_DECOMPRESS_H_

#include <string>
#include "util.h"

namespace spring
{

void decompress(const std::string &temp_dir, const std::string &outfile_1,
const std::string &outfile_2, const compression_params &cp, const int &num_thr);

void decompress_long(const std::string &temp_dir, const std::string &outfile_1,
const std::string &outfile_2, const compression_params &cp, const int &num_thr);

} // namespace spring

#endif // SPRING_DECOMPRESS_H_
