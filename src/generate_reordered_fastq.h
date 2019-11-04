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

#ifndef SPRING_GENERATE_REORDERED_FASTQ_H_
#define SPRING_GENERATE_REORDERED_FASTQ_H_

#include <string>
#include <vector>

#include "util.h"

namespace spring {

void generate_reordered_fastq(const std::string &temp_dir,
                              const compression_params &cp,
                              const std::vector<std::string> &infile_vector,
                              const std::vector<std::string> &outfile_vector,
                              const bool gzipped_output_flag,
                              const bool gzipped_input_flag);



void generate_order_pe(const std::string &file_order, uint32_t *order_array,
                       const uint32_t &numreads);

void generate_order_se(const std::string &file_order, uint32_t *order_array,
                       const uint32_t &numreads);

} // namespace spring

# endif // SPRING_GENERATE_REORDERED_FASTQ_H_
