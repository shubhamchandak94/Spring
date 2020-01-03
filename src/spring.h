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

#ifndef SPRING_SPRING_H_
#define SPRING_SPRING_H_

#include <string>
#include "util.h"

namespace spring {

void compress(const std::string &temp_dir,
              const std::vector<std::string> &infile_vec,
              const std::vector<std::string> &outfile_vec, const int &num_thr,
              const bool &pairing_only_flag, const bool &no_quality_flag,
              const bool &no_ids_flag,
              const std::vector<std::string> &quality_opts,
              const bool &long_flag, const bool &gzip_flag, const bool &fasta_flag);

void decompress(const std::string &temp_dir,
                const std::vector<std::string> &infile_vec,
                const std::vector<std::string> &outfile_vec, const int &num_thr,
                const std::vector<uint64_t> &decompress_range_vec,
                const bool &gzip_flag);

std::string random_string(size_t length);

}  // namespace spring

#endif  // SPRING_SPRING_H_
