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

#include <algorithm>
#include <boost/filesystem.hpp>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>  // std::setw
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "call_template_functions.h"
#include "encoder.h"
#include "generate_reordered_fastq.h"
#include "params.h"
#include "pe_encode.h"
#include "preprocess.h"
#include "reorder.h"
#include "spring.h"
#include "util.h"

namespace spring {

void spring_reorder(const std::string &temp_dir,
              const std::vector<std::string> &infile_vec,
              const std::vector<std::string> &outfile_vec, const int &num_thr,
              const bool &gzipped_input_flag, const bool &gzipped_output_flag) {
  //
  // Ensure that omp parallel regions are executed with the requested
  // #threads.
  //
  int omp_dyn_var = omp_get_dynamic();
  omp_set_dynamic(0);
  std::cout << "Starting compression...\n";
  auto compression_start = std::chrono::steady_clock::now();

  std::string infile_1, infile_2, outfile_1, outfile_2;
  bool paired_end;
  // Check options
  switch (infile_vec.size()) {
    case 0:
      throw std::runtime_error("No input file specified");
      break;
    case 1:
      paired_end = false;
      infile_1 = infile_vec[0];
      break;
    case 2:
      paired_end = true;
      infile_1 = infile_vec[0];
      infile_2 = infile_vec[1];
      break;
    default:
      throw std::runtime_error("Too many (>2) input files specified");
  }

  switch (outfile_vec.size()) {
    case 0:
      throw std::runtime_error("No output file specified");
      break;
    case 1:
      if (!paired_end) outfile_1 = outfile_vec[0];
      if (paired_end) {
        outfile_1 = outfile_vec[0] + ".1";
        outfile_2 = outfile_vec[0] + ".2";
      }
      break;
    case 2:
      if (!paired_end) {
        std::cerr << "WARNING: Two output files provided for single end data. "
                     "Output will be written to the first file provided.\n";
        outfile_1 = outfile_vec[0];
      } else {
        outfile_1 = outfile_vec[0];
        outfile_2 = outfile_vec[1];
      }
      break;
    default:
      throw std::runtime_error("Too many (>2) output files specified");
  }

  compression_params *cp_ptr = new compression_params;
  memset(cp_ptr, 0, sizeof(compression_params));  // remove valgrind error
  compression_params &cp = *cp_ptr;
  cp.paired_end = paired_end;
  cp.num_reads_per_block = NUM_READS_PER_BLOCK;
  cp.num_thr = num_thr;

  std::cout << "Preprocessing ...\n";
  auto preprocess_start = std::chrono::steady_clock::now();
  preprocess(infile_1, infile_2, temp_dir, cp, gzipped_input_flag);
  auto preprocess_end = std::chrono::steady_clock::now();
  std::cout << "Preprocessing done!\n";
  std::cout << "Time for this step: "
            << std::chrono::duration_cast<std::chrono::seconds>(
                   preprocess_end - preprocess_start)
                   .count()
            << " s\n";

  std::cout << "Reordering ...\n";
  auto reorder_start = std::chrono::steady_clock::now();
  call_reorder(temp_dir, cp);
  auto reorder_end = std::chrono::steady_clock::now();
  std::cout << "Reordering done!\n";
  std::cout << "Time for this step: "
            << std::chrono::duration_cast<std::chrono::seconds>(reorder_end -
                                                                reorder_start)
                   .count()
            << " s\n";

  std::cout << "Realigning singletons ...\n";
  auto encoder_start = std::chrono::steady_clock::now();
  call_encoder(temp_dir, cp);
  auto encoder_end = std::chrono::steady_clock::now();
  std::cout << "Realigning singletons done!\n";
  std::cout << "Time for this step: "
            << std::chrono::duration_cast<std::chrono::seconds>(encoder_end -
                                                                encoder_start)
                   .count()
            << " s\n";

  std::cout << "Generating reordered FASTQ ...\n";
  auto grf_start = std::chrono::steady_clock::now();
  std::vector<std::string> infile_vector = {infile_1, infile_2};
  std::vector<std::string> outfile_vector = {outfile_1, outfile_2};
  generate_reordered_fastq(temp_dir, cp, infile_vector, outfile_vector,
                           gzipped_output_flag, gzipped_input_flag);
  auto grf_end = std::chrono::steady_clock::now();
  std::cout << "Generating reordered FASTQ done!\n";
  std::cout << "Time for this step: "
            << std::chrono::duration_cast<std::chrono::seconds>(grf_end -
                                                                grf_start)
                   .count()
            << " s\n";

  delete cp_ptr;
  auto compression_end = std::chrono::steady_clock::now();
  std::cout << "Done!\n";
  std::cout << "Total time: "
            << std::chrono::duration_cast<std::chrono::seconds>(
                   compression_end - compression_start)
                   .count()
            << " s\n";

  return;
}

std::string random_string(size_t length) {
  auto randchar = []() -> char {
    const char charset[] =
        "0123456789"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz";
    const size_t max_index = (sizeof(charset) - 1);
    return charset[rand() % max_index];
  };
  std::string str(length, 0);
  std::generate_n(str.begin(), length, randchar);
  return str;
}

}  // namespace spring
