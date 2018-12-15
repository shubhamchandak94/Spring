#ifndef SPRING_SPRING_H_
#define SPRING_SPRING_H_

#include <string>
#include "util.h"

namespace spring {

void compress(const std::string &temp_dir, const std::vector<std::string> &infile_vec,
              const std::vector<std::string> &outfile_vec, const int &num_thr,
              const bool &pairing_only_flag, const bool &no_quality_flag, const bool &no_ids_flag,
              const std::vector<std::string> &quality_opts, const bool &long_flag, const bool &gzip_flag);

void decompress(const std::string &temp_dir, const std::vector<std::string> &infile_vec,
                const std::vector<std::string> &outfile_vec, const int &num_thr, const std::vector<uint64_t> &decompress_range_vec);

void call_reorder(const std::string &temp_dir, compression_params &cp);

void call_encoder(const std::string &temp_dir, compression_params &cp);

std::string random_string(size_t length);

}  // namespace spring

#endif  // SPRING_SPRING_H_
