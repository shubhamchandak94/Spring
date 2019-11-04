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

#include "preprocess.h"
#include <omp.h>
#include <algorithm>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include "libbsc/bsc.h"
#include "params.h"
#include "util.h"

namespace spring {

void preprocess(const std::string &infile_1, const std::string &infile_2,
                const std::string &temp_dir, compression_params &cp,
                const bool &gzip_flag) {
  std::string infile[2] = {infile_1, infile_2};
  std::string outfileclean[2];
  std::string outfileN[2];
  std::string outfileorderN[2];
  std::string outfileid[2];
  std::string outfilequality[2];
  std::string outfileread[2];
  std::string outfilereadlength[2];
  std::string basedir = temp_dir;
  outfileclean[0] = basedir + "/input_clean_1.dna";
  outfileclean[1] = basedir + "/input_clean_2.dna";
  outfileN[0] = basedir + "/input_N.dna";
  outfileN[1] = basedir + "/input_N.dna.2";
  outfileorderN[0] = basedir + "/read_order_N.bin";
  outfileorderN[1] = basedir + "/read_order_N.bin.2";
  outfileid[0] = basedir + "/id_1";
  outfileid[1] = basedir + "/id_2";
  outfilequality[0] = basedir + "/quality_1";
  outfilequality[1] = basedir + "/quality_2";
  outfileread[0] = basedir + "/read_1";
  outfileread[1] = basedir + "/read_2";
  outfilereadlength[0] = basedir + "/readlength_1";
  outfilereadlength[1] = basedir + "/readlength_2";

  std::ifstream fin_f[2];
  std::ofstream fout_clean[2];
  std::ofstream fout_N[2];
  std::ofstream fout_order_N[2];
  std::ofstream fout_id[2];
  std::ofstream fout_quality[2];
  std::istream *fin[2] = {&fin_f[0], &fin_f[1]};
  boost::iostreams::filtering_streambuf<boost::iostreams::input> *inbuf[2];

  for (int j = 0; j < 2; j++) {
    if (j == 1 && !cp.paired_end) continue;
    if (gzip_flag) {
      fin_f[j].open(infile[j], std::ios_base::binary);
      inbuf[j] =
          new boost::iostreams::filtering_streambuf<boost::iostreams::input>;
      inbuf[j]->push(boost::iostreams::gzip_decompressor());
      inbuf[j]->push(fin_f[j]);
      fin[j] = new std::istream(inbuf[j]);
    } else {
      fin_f[j].open(infile[j]);
    }
    if (!cp.long_flag) {
      fout_clean[j].open(outfileclean[j]);
      fout_N[j].open(outfileN[j]);
      fout_order_N[j].open(outfileorderN[j], std::ios::binary);
      if (!cp.preserve_order) {
        if (cp.preserve_id) fout_id[j].open(outfileid[j]);
        if (cp.preserve_quality) fout_quality[j].open(outfilequality[j]);
      }
    }
  }

  uint32_t max_readlen = 0;
  uint64_t num_reads[2] = {0, 0};
  uint64_t num_reads_clean[2] = {0, 0};
  uint32_t num_reads_per_block;
  if (!cp.long_flag)
    num_reads_per_block = cp.num_reads_per_block;
  else
    num_reads_per_block = cp.num_reads_per_block_long;
  uint8_t paired_id_code = 0;
  bool paired_id_match = false;

  char *quality_binning_table = new char[128];
  if (cp.ill_bin_flag) generate_illumina_binning_table(quality_binning_table);
  if (cp.bin_thr_flag)
    generate_binary_binning_table(quality_binning_table, cp.bin_thr_thr,
                                  cp.bin_thr_high, cp.bin_thr_low);

  // Check that we were able to open the input files and also look for
  // paired end matching ids if relevant
  if (!fin_f[0].is_open()) throw std::runtime_error("Error opening input file");
  if (cp.paired_end) {
    if (!fin_f[1].is_open())
      throw std::runtime_error("Error opening input file");
    if (cp.preserve_id) {
      std::string id_1, id_2;
      std::getline(*fin[0], id_1);
      std::getline(*fin[1], id_2);
      paired_id_code = find_id_pattern(id_1, id_2);
      if (paired_id_code != 0) paired_id_match = true;
      // bring back to start
      if (gzip_flag) {
        for (int j = 0; j < 2; j++) {
          fin_f[j].close();
          delete fin[j];
          delete inbuf[j];
          fin_f[j].open(infile[j], std::ios_base::binary);
          inbuf[j] = new boost::iostreams::filtering_streambuf<
              boost::iostreams::input>;
          inbuf[j]->push(boost::iostreams::gzip_decompressor());
          inbuf[j]->push(fin_f[j]);
          fin[j] = new std::istream(inbuf[j]);
        }
      } else {
        fin_f[0].seekg(0);
        fin_f[1].seekg(0);
      }
    }
  }
  uint64_t num_reads_per_step = (uint64_t)cp.num_thr * num_reads_per_block;
  std::string *read_array = new std::string[num_reads_per_step];
  std::string *id_array_1 = new std::string[num_reads_per_step];
  std::string *id_array_2 = new std::string[num_reads_per_step];
  std::string *quality_array = new std::string[num_reads_per_step];
  bool *read_contains_N_array = new bool[num_reads_per_step];
  uint32_t *read_lengths_array = new uint32_t[num_reads_per_step];
  bool *paired_id_match_array = new bool[cp.num_thr];

  omp_set_num_threads(cp.num_thr);

  uint32_t num_blocks_done = 0;

  while (true) {
    bool done[2] = {true, true};
    for (int j = 0; j < 2; j++) {
      if (j == 1 && !cp.paired_end) continue;
      done[j] = false;
      std::string *id_array = (j == 0) ? id_array_1 : id_array_2;
      uint32_t num_reads_read = read_fastq_block(
          fin[j], id_array, read_array, quality_array, num_reads_per_step);
      if (num_reads_read < num_reads_per_step) done[j] = true;
      if (num_reads_read == 0) continue;
      if (num_reads[0] + num_reads[1] + num_reads_read > MAX_NUM_READS) {
        std::cerr << "Max number of reads allowed is " << MAX_NUM_READS << "\n";
        throw std::runtime_error("Too many reads.");
      }
#pragma omp parallel
      {
        bool done = false;
        uint64_t tid = omp_get_thread_num();
        if (j == 1) paired_id_match_array[tid] = paired_id_match;
        if (tid * num_reads_per_block >= num_reads_read) done = true;
        uint32_t num_reads_thr = std::min((uint64_t)num_reads_read,
                                          (tid + 1) * num_reads_per_block) -
                                 tid * num_reads_per_block;
        std::ofstream fout_readlength;
        if (!done) {
          if (cp.long_flag)
            fout_readlength.open(outfilereadlength[j] + "." +
                                     std::to_string(num_blocks_done + tid),
                                 std::ios::binary);
          // check if reads and qualities have equal lengths
          for (uint32_t i = tid * num_reads_per_block;
               i < tid * num_reads_per_block + num_reads_thr; i++) {
            size_t len = read_array[i].size();
            if (cp.long_flag && len > MAX_READ_LEN_LONG) {
              std::cerr << "Max read length for long mode is "
                        << MAX_READ_LEN_LONG << ", but found read of length "
                        << len << "\n";
              throw std::runtime_error("Too long read length");
            }
            if (!cp.long_flag && len > MAX_READ_LEN) {
              std::cerr << "Max read length without long mode is "
                        << MAX_READ_LEN << ", but found read of length " << len
                        << "\n";
              throw std::runtime_error(
                  "Too long read length (please try --long/-l flag).");
            }
            if (cp.preserve_quality && (quality_array[i].size() != len))
              throw std::runtime_error(
                  "Read length does not match quality length.");
            if (cp.preserve_id && (id_array[i].size() == 0))
              throw std::runtime_error("Identifier of length 0 detected.");
            read_lengths_array[i] = (uint32_t)len;

            // mark reads with N
            if (!cp.long_flag)
              read_contains_N_array[i] =
                  (read_array[i].find('N') != std::string::npos);

            // Write read length to a file (for long mode)
            if (cp.long_flag)
              fout_readlength.write((char *)&read_lengths_array[i],
                                    sizeof(uint32_t));

            if (j == 1 && paired_id_match_array[tid])
              paired_id_match_array[tid] = check_id_pattern(
                  id_array_1[i], id_array_2[i], paired_id_code);
          }
          if (cp.long_flag) fout_readlength.close();
          // apply binning (if asked to do so)
          if (cp.preserve_quality && (cp.ill_bin_flag || cp.bin_thr_flag))
            quantize_quality(quality_array + tid * num_reads_per_block,
                             num_reads_thr, quality_binning_table);

          // apply qvz quantization (if flag specified and order preserving)
          if (cp.preserve_quality && cp.qvz_flag && cp.preserve_order)
            quantize_quality_qvz(
                quality_array + tid * num_reads_per_block, num_reads_thr,
                read_lengths_array + tid * num_reads_per_block, cp.qvz_ratio);
          if (!cp.long_flag) {
            if (cp.preserve_order) {
              // Compress ids
              if (cp.preserve_id) {
                std::string outfile_name =
                    outfileid[j] + "." + std::to_string(num_blocks_done + tid);
                compress_id_block(outfile_name.c_str(),
                                  id_array + tid * num_reads_per_block,
                                  num_reads_thr);
              }
              // Compress qualities
              if (cp.preserve_quality) {
                std::string outfile_name =
                    outfilequality[j] + "." +
                    std::to_string(num_blocks_done + tid);
                bsc::BSC_str_array_compress(
                    outfile_name.c_str(),
                    quality_array + tid * num_reads_per_block, num_reads_thr,
                    read_lengths_array + tid * num_reads_per_block);
              }
            }
          } else {
            // Compress read lengths file
            std::string infile_name = outfilereadlength[j] + "." +
                                      std::to_string(num_blocks_done + tid);
            std::string outfile_name = outfilereadlength[j] + "." +
                                       std::to_string(num_blocks_done + tid) +
                                       ".bsc";
            bsc::BSC_compress(infile_name.c_str(), outfile_name.c_str());
            remove(infile_name.c_str());
            // Compress ids
            if (cp.preserve_id) {
              std::string outfile_name =
                  outfileid[j] + "." + std::to_string(num_blocks_done + tid);
              compress_id_block(outfile_name.c_str(),
                                id_array + tid * num_reads_per_block,
                                num_reads_thr);
            }
            // Compress qualities
            if (cp.preserve_quality) {
              std::string outfile_name = outfilequality[j] + "." +
                                         std::to_string(num_blocks_done + tid);
              bsc::BSC_str_array_compress(
                  outfile_name.c_str(),
                  quality_array + tid * num_reads_per_block, num_reads_thr,
                  read_lengths_array + tid * num_reads_per_block);
            }
            // Compress reads
            outfile_name =
                outfileread[j] + "." + std::to_string(num_blocks_done + tid);
            bsc::BSC_str_array_compress(
                outfile_name.c_str(), read_array + tid * num_reads_per_block,
                num_reads_thr, read_lengths_array + tid * num_reads_per_block);
          }
        }  // if(!done)
      }    // omp parallel
      // if id match not found in any thread, set to false
      if (cp.paired_end && (j == 1)) {
        if (paired_id_match)
          for (int tid = 0; tid < cp.num_thr; tid++)
            paired_id_match &= paired_id_match_array[tid];
        if (!paired_id_match) paired_id_code = 0;
      }
      if (!cp.long_flag) {
        // write reads and read_order_N to respective files
        for (uint32_t i = 0; i < num_reads_read; i++) {
          if (!read_contains_N_array[i]) {
            fout_clean[j] << read_array[i] << "\n";
            num_reads_clean[j]++;
          } else {
            uint32_t pos_N = num_reads[j] + i;
            fout_order_N[j].write((char *)&pos_N, sizeof(uint32_t));
            fout_N[j] << read_array[i] << "\n";
          }
        }
        if (!cp.preserve_order) {
          // write qualities to file
          for (uint32_t i = 0; i < num_reads_read; i++)
            fout_quality[j] << quality_array[i] << "\n";
          // write ids to file
          for (uint32_t i = 0; i < num_reads_read; i++)
            fout_id[j] << id_array[i] << "\n";
        }
      }
      num_reads[j] += num_reads_read;
      max_readlen =
          std::max(max_readlen,
                   *(std::max_element(read_lengths_array,
                                      read_lengths_array + num_reads_read)));
    }
    if (cp.paired_end)
      if (num_reads[0] != num_reads[1])
        throw std::runtime_error(
            "Number of reads in paired files do not match.");
    if (done[0] && done[1]) break;
    num_blocks_done += cp.num_thr;
  }

  delete[] read_array;
  delete[] id_array_1;
  delete[] id_array_2;
  delete[] quality_array;
  delete[] read_contains_N_array;
  delete[] read_lengths_array;
  delete[] quality_binning_table;
  delete[] paired_id_match_array;
  if (gzip_flag) {
    for (int j = 0; j < 2; j++) {
      if (j == 1 && !cp.paired_end) continue;
      delete fin[j];
      delete inbuf[j];
    }
  }
  // close files
  if (cp.long_flag) {
    fin_f[0].close();
    if (cp.paired_end) fin_f[1].close();
  } else {
    for (int j = 0; j < 2; j++) {
      if (j == 1 && !cp.paired_end) continue;
      fin_f[j].close();
      fout_clean[j].close();
      fout_N[j].close();
      fout_order_N[j].close();
      if (!cp.preserve_order) {
        if (cp.preserve_id) fout_id[j].close();
        if (cp.preserve_quality) fout_quality[j].close();
      }
    }
  }
  if (num_reads[0] == 0) throw std::runtime_error("No reads found.");

  if (!cp.long_flag && cp.paired_end) {
    // merge input_N and input_order_N for the two files
    std::ofstream fout_N(outfileN[0], std::ios::app);
    std::ifstream fin_N(outfileN[1]);
    fout_N << fin_N.rdbuf();
    fout_N.close();
    fin_N.close();
    remove(outfileN[1].c_str());
    std::ofstream fout_order_N(outfileorderN[0],
                               std::ios::app | std::ios::binary);
    std::ifstream fin_order_N(outfileorderN[1], std::ios::binary);
    uint32_t num_N_file_2 = num_reads[1] - num_reads_clean[1];
    uint32_t order_N;
    for (uint32_t i = 0; i < num_N_file_2; i++) {
      fin_order_N.read((char *)&order_N, sizeof(uint32_t));
      order_N += num_reads[0];
      fout_order_N.write((char *)&order_N, sizeof(uint32_t));
    }
    fin_order_N.close();
    fout_order_N.close();
    remove(outfileorderN[1].c_str());
  }

  if (cp.paired_end && paired_id_match) {
    // delete id files for second file since we found a pattern
    if (!cp.long_flag && !cp.preserve_order) {
      remove(outfileid[1].c_str());
    } else {
      uint32_t num_blocks =
          1 +
          (num_reads[0] - 1) /
              num_reads_per_block;  // ceil of num_reads[0]/num_reads_per_block
      for (uint32_t i = 0; i < num_blocks; i++)
        remove((outfileid[1] + "." + std::to_string(i)).c_str());
    }
  }
  cp.paired_id_code = paired_id_code;
  cp.paired_id_match = paired_id_match;
  cp.num_reads = num_reads[0] + num_reads[1];
  cp.num_reads_clean[0] = num_reads_clean[0];
  cp.num_reads_clean[1] = num_reads_clean[1];
  cp.max_readlen = max_readlen;

  std::cout << "Max Read length: " << cp.max_readlen << "\n";
  std::cout << "Total number of reads: " << cp.num_reads << "\n";

  if (!cp.long_flag)
    std::cout << "Total number of reads without N: "
              << cp.num_reads_clean[0] + cp.num_reads_clean[1] << "\n";
  if (cp.preserve_id && cp.paired_end)
    std::cout << "Paired id match code: " << (int)cp.paired_id_code << "\n";
}

}  // namespace spring
