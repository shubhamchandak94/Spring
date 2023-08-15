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
#include "decompress.h"
#include "encoder.h"
#include "params.h"
#include "pe_encode.h"
#include "preprocess.h"
#include "reorder.h"
#include "reorder_compress_quality_id.h"
#include "reorder_compress_streams.h"
#include "spring.h"
#include "util.h"

namespace spring {

void compress(const std::string &temp_dir,
              const std::vector<std::string> &infile_vec,
              const std::vector<std::string> &outfile_vec, const int &num_thr,
              const bool &pairing_only_flag, const bool &no_quality_flag,
              const bool &no_ids_flag,
              const std::vector<std::string> &quality_opts,
              const bool &long_flag, const bool &gzip_flag, const bool &fasta_flag) {
  //
  // Ensure that omp parallel regions are executed with the requested
  // #threads.
  //
  omp_set_dynamic(0);

  std::cout << __FILE__ << ":" __LINE__ << ":DEBUG Starting compression...\n";
  auto compression_start = std::chrono::steady_clock::now();

  std::string infile_1, infile_2, outfile;
  bool paired_end, preserve_quality, preserve_id, preserve_order;
  // Check options
  preserve_order = !pairing_only_flag;
  preserve_id = !no_ids_flag;
  preserve_quality = !no_quality_flag;
  if (fasta_flag)
    preserve_quality = false;
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
  if (outfile_vec.size() == 1)
    outfile = outfile_vec[0];
  else {
     std::cerr << __LINE__ << ": For compression you must provide single output file\n";
    throw std::runtime_error("Number of output files not equal to 1");
  }

  compression_params *cp_ptr = new compression_params;
  memset(cp_ptr, 0, sizeof(compression_params));  // remove valgrind error
  compression_params &cp = *cp_ptr;
  cp.paired_end = paired_end;
  cp.preserve_order = preserve_order;
  cp.preserve_id = preserve_id;
  cp.preserve_quality = preserve_quality;
  cp.long_flag = long_flag;
  cp.num_reads_per_block = NUM_READS_PER_BLOCK;
  cp.num_reads_per_block_long = NUM_READS_PER_BLOCK_LONG;
  cp.num_thr = num_thr;

  if (preserve_quality) {
    if (quality_opts.empty()) {
      cp.qvz_flag = cp.ill_bin_flag = cp.bin_thr_flag = false;
    } else if (quality_opts[0] == "lossless") {
      cp.qvz_flag = cp.ill_bin_flag = cp.bin_thr_flag = false;
    } else if (quality_opts[0] == "qvz") {
      if (quality_opts.size() != 2) {
        throw std::runtime_error("Invalid quality options.");
      } else {
        cp.qvz_ratio = atof(quality_opts[1].c_str());
        if (cp.qvz_ratio == 0.0) {
          throw std::runtime_error("Invalid qvz ratio provided.");
        }
      }
      cp.qvz_flag = true;
      cp.ill_bin_flag = cp.bin_thr_flag = false;
    } else if (quality_opts[0] == "ill_bin") {
      cp.ill_bin_flag = true;
      cp.qvz_flag = cp.bin_thr_flag = false;
    } else if (quality_opts[0] == "binary") {
      if (quality_opts.size() != 4) {
        throw std::runtime_error("Invalid quality options.");
      } else {
        cp.bin_thr_thr = atoi(quality_opts[1].c_str());
        cp.bin_thr_high = atoi(quality_opts[2].c_str());
        cp.bin_thr_low = atoi(quality_opts[3].c_str());
        if (cp.bin_thr_high < cp.bin_thr_thr ||
            cp.bin_thr_low > cp.bin_thr_thr ||
            cp.bin_thr_high < cp.bin_thr_low) {
          throw std::runtime_error(
              "Options do not satisfy low <= thr <= high.");
        }
      }
      cp.qvz_flag = cp.ill_bin_flag = false;
      cp.bin_thr_flag = true;
    } else {
      throw std::runtime_error("Invalid quality options.");
    }
  }

  std::cout << __LINE__ << ": Preprocessing ...\n";
  auto preprocess_start = std::chrono::steady_clock::now();
  preprocess(infile_1, infile_2, temp_dir, cp, gzip_flag, fasta_flag);
  auto preprocess_end = std::chrono::steady_clock::now();
  std::cout << "Preprocessing done!\n";
  std::cout << "Time for this step: "
            << std::chrono::duration_cast<std::chrono::seconds>(
                   preprocess_end - preprocess_start)
                   .count()
            << " s\n";
  std::cout << "Temporary directory size: " << get_directory_size(temp_dir) << "\n";

  if (!long_flag) {
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

    std::cout << "temp_dir size: " << get_directory_size(temp_dir) << "\n";

    std::cout << "Encoding ...\n";
    auto encoder_start = std::chrono::steady_clock::now();
    call_encoder(temp_dir, cp);
    auto encoder_end = std::chrono::steady_clock::now();
    std::cout << "Encoding done!\n";
    std::cout << "Time for this step: "
              << std::chrono::duration_cast<std::chrono::seconds>(encoder_end -
                                                                  encoder_start)
                     .count()
              << " s\n";
    std::cout << "Temporary directory size: " << get_directory_size(temp_dir) << "\n";

    if (!preserve_order && (preserve_quality || preserve_id)) {
      std::cout << "Reordering and compressing quality and/or ids ...\n";
      auto rcqi_start = std::chrono::steady_clock::now();
      reorder_compress_quality_id(temp_dir, cp);
      auto rcqi_end = std::chrono::steady_clock::now();
      std::cout << "Reordering and compressing quality and/or ids done!\n";
      std::cout << "Time for this step: "
                << std::chrono::duration_cast<std::chrono::seconds>(rcqi_end -
                                                                    rcqi_start)
                       .count()
                << " s\n";
      std::cout << "Temporary directory size: " << get_directory_size(temp_dir) << "\n";
    }

    if (!preserve_order && paired_end) {
      std::cout << "Encoding pairing information ...\n";
      auto pe_encode_start = std::chrono::steady_clock::now();
      pe_encode(temp_dir, cp);
      auto pe_encode_end = std::chrono::steady_clock::now();
      std::cout << "Encoding pairing information done!\n";
      std::cout << "Time for this step: "
                << std::chrono::duration_cast<std::chrono::seconds>(
                       pe_encode_end - pe_encode_start)
                       .count()
                << " s\n";
      std::cout << "Temporary directory size: " << get_directory_size(temp_dir) << "\n";
    }

    std::cout << "Reordering and compressing streams ...\n";
    auto rcs_start = std::chrono::steady_clock::now();
    reorder_compress_streams(temp_dir, cp);
    auto rcs_end = std::chrono::steady_clock::now();
    std::cout << "Reordering and compressing streams done!\n";
    std::cout << "Time for this step: "
              << std::chrono::duration_cast<std::chrono::seconds>(rcs_end -
                                                                  rcs_start)
                     .count()
              << " s\n";
    std::cout << "Temporary directory size: " << get_directory_size(temp_dir) << "\n";
  }

  // Write compression params to a file
  std::string compression_params_file = temp_dir + "/cp.bin";
  std::ofstream f_cp(compression_params_file, std::ios::binary);
  f_cp.write((char *)&cp, sizeof(compression_params));
  f_cp.close();

  // Print out sizes of reads, quality and id after compression
  namespace fs = boost::filesystem;
  uint64_t size_read = 0;
  uint64_t size_quality = 0;
  uint64_t size_id = 0;
  fs::path p{temp_dir};
  fs::directory_iterator itr{p};
  for (; itr != fs::directory_iterator{}; ++itr) {
    std::string current_file = itr->path().filename().string();
    switch (current_file[0]) {
      case 'r':
        size_read += fs::file_size(itr->path());
        break;
      case 'q':
        size_quality += fs::file_size(itr->path());
        break;
      case 'i':
        size_id += fs::file_size(itr->path());
        break;
    }
  }
  std::cout << "\n";
  std::cout << "Sizes of streams after compression: \n";
  std::cout << "Reads:      " << std::setw(12) << size_read << " bytes\n";
  std::cout << "Quality:    " << std::setw(12) << size_quality << " bytes\n";
  std::cout << "ID:         " << std::setw(12) << size_id << " bytes\n";

  auto tar_start = std::chrono::steady_clock::now();
  std::cout << "Creating tar archive ...";
  std::string tar_command = "tar -cf " + outfile + " -C " + temp_dir + " . ";
  int tar_status = std::system(tar_command.c_str());
  if (tar_status != 0)
    throw std::runtime_error("Error occurred during tar archive generation.");
  std::cout << "Tar archive done!\n";
  auto tar_end = std::chrono::steady_clock::now();
  std::cout << "Time for this step: "
            << std::chrono::duration_cast<std::chrono::seconds>(tar_end -
                                                                tar_start)
                   .count()
            << " s\n";

  delete cp_ptr;
  auto compression_end = std::chrono::steady_clock::now();
  std::cout << "Compression done!\n";
  std::cout << "Total time for compression: "
            << std::chrono::duration_cast<std::chrono::seconds>(
                   compression_end - compression_start)
                   .count()
            << " s\n";

  fs::path p1{outfile};
  std::cout << "\n";
  std::cout << "Total size: " << std::setw(12) << fs::file_size(p1)
            << " bytes\n";
  return;
}

void decompress(const std::string &temp_dir,
                const std::vector<std::string> &infile_vec,
                const std::vector<std::string> &outfile_vec, const int &num_thr,
                const std::vector<uint64_t> &decompress_range_vec,
                const bool &gzip_flag, const int &gzip_level) {
  //
  // Ensure that omp parallel regions are executed with the requested
  // #threads.
  //
  omp_set_dynamic(0);

  std::cout << "Starting decompression...\n";
  auto decompression_start = std::chrono::steady_clock::now();
  compression_params *cp_ptr = new compression_params;
  compression_params &cp = *cp_ptr;

  std::string infile, outfile_1, outfile_2;

  if (infile_vec.size() == 1)
    infile = infile_vec[0];
  else
    throw std::runtime_error("Number of input files not equal to 1");

  std::cout << "Untarring tar archive ...\n";
  std::string untar_command = "tar -xf " + infile + " -C " + temp_dir;
  int untar_status = std::system(untar_command.c_str());
  if (untar_status != 0)
    throw std::runtime_error("Error occurred during untarring.");
  std::cout << "Untarring archive done!\n";

  // Read compression params
  std::string compression_params_file = temp_dir + "/cp.bin";
  std::ifstream f_cp(compression_params_file, std::ios::binary);
  if (!f_cp.is_open()) throw std::runtime_error("Can't open parameter file.");
  f_cp.read((char *)&cp, sizeof(compression_params));
  if (!f_cp.good())
    throw std::runtime_error("Can't read compression parameters.");
  f_cp.close();

  bool paired_end = cp.paired_end;
  bool long_flag = cp.long_flag;

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
                     "Output will be written to the first file provided.";
        outfile_1 = outfile_vec[0];
      } else {
        outfile_1 = outfile_vec[0];
        outfile_2 = outfile_vec[1];
      }
      break;
    default:
      throw std::runtime_error("Too many (>2) output files specified");
  }
  uint64_t num_read_pairs = paired_end ? cp.num_reads / 2 : cp.num_reads;
  uint64_t start_num = 0;
  uint64_t end_num = num_read_pairs;
  if (decompress_range_vec.size() != 0) {
    if (decompress_range_vec.size() != 2)
      throw std::runtime_error("Invalid decompression range parameters.");
    if (decompress_range_vec[0] == 0 ||
        decompress_range_vec[0] > decompress_range_vec[1] ||
        decompress_range_vec[0] > num_read_pairs ||
        decompress_range_vec[1] > num_read_pairs)
      throw std::runtime_error("Invalid decompression range parameters.");
    start_num = decompress_range_vec[0] - 1;
    end_num = decompress_range_vec[1];
  }

  std::cout << "Decompressing ...\n";
  if (long_flag)
    decompress_long(temp_dir, outfile_1, outfile_2, cp, num_thr, start_num,
                    end_num, gzip_flag, gzip_level);
  else
    decompress_short(temp_dir, outfile_1, outfile_2, cp, num_thr, start_num,
                     end_num, gzip_flag, gzip_level);

  delete cp_ptr;
  auto decompression_end = std::chrono::steady_clock::now();
  std::cout << "Decompression done!\n";
  std::cout << "Total time for decompression: "
            << std::chrono::duration_cast<std::chrono::seconds>(
                   decompression_end - decompression_start)
                   .count()
            << " s\n";
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
