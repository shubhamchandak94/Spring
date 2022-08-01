/*
* Copyright 2018 University of Illinois Board of Trustees and Stanford
University. All Rights Reserved.
* Licensed under the “Non-exclusive Research Use License for SPRING Software”
license (the "License");
* You may not use this file except in compliance with the License.
* The License is included in the distribution as license.pdf file.

* Software distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
limitations under the License.
*/

#ifndef SPRING_DECOMPRESS_H_
#define SPRING_DECOMPRESS_H_

#include <string>
#include "util.h"

namespace spring {

void decompress_short(const std::string &temp_dir, const std::string &outfile_1,
                      const std::string &outfile_2,
                      const compression_params &cp, const int &num_thr,
                      const uint64_t &start_num, const uint64_t &end_num,
                      const bool &gzip_flag, const int &gzip_level);

void decompress_long(const std::string &temp_dir, const std::string &outfile_1,
                     const std::string &outfile_2, const compression_params &cp,
                     const int &num_thr, const uint64_t &start_num,
                     const uint64_t &end_num, const bool &gzip_flag,
                     const int &gzip_level);

void decompress_unpack_seq(const std::string &infile_seq, const int &num_thr_e,
                           const int &num_thr);

}  // namespace spring

#endif  // SPRING_DECOMPRESS_H_
