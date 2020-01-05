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

#include "params.h"
#include "util.h"

namespace spring {

void preprocess(const std::string &infile_1, const std::string &infile_2,
                const std::string &temp_dir, compression_params &cp,
                const bool &gzip_flag, const bool &fasta_flag) {
  std::string infile[2] = {infile_1, infile_2};
  std::string outfileclean[2];
  std::string outfileN[2];
  std::string outfileorderN[2];
  std::string basedir = temp_dir;
  outfileclean[0] = basedir + "/input_clean_1.dna";
  outfileclean[1] = basedir + "/input_clean_2.dna";
  outfileN[0] = basedir + "/input_N.dna";
  outfileN[1] = basedir + "/input_N.dna.2";
  outfileorderN[0] = basedir + "/read_order_N.bin";
  outfileorderN[1] = basedir + "/read_order_N.bin.2";

  std::ifstream fin_f[2];
  std::ofstream fout_clean[2];
  std::ofstream fout_N[2];
  std::ofstream fout_order_N[2];
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
    fout_clean[j].open(outfileclean[j],std::ios::binary);
    fout_N[j].open(outfileN[j],std::ios::binary);
    fout_order_N[j].open(outfileorderN[j], std::ios::binary);
  }

  uint32_t max_readlen = 0;
  uint64_t num_reads[2] = {0, 0};
  uint64_t num_reads_clean[2] = {0, 0};
  uint32_t num_reads_per_block;
  num_reads_per_block = cp.num_reads_per_block;

  // Check that we were able to open the input files and also look for
  // paired end matching ids if relevant
  if (!fin_f[0].is_open()) throw std::runtime_error("Error opening input file");
  if (cp.paired_end) {
    if (!fin_f[1].is_open())
      throw std::runtime_error("Error opening input file");
  }

  uint64_t num_reads_per_step = (uint64_t)cp.num_thr * num_reads_per_block;
  std::string *read_array = new std::string[num_reads_per_step];
  std::string *id_array = new std::string[num_reads_per_step];
  std::string *quality_array = new std::string[num_reads_per_step];
  bool *read_contains_N_array = new bool[num_reads_per_step];
  uint32_t *read_lengths_array = new uint32_t[num_reads_per_step];

  omp_set_num_threads(cp.num_thr);

  uint32_t num_blocks_done = 0;

  while (true) {
    bool done[2] = {true, true};
    for (int j = 0; j < 2; j++) {
      if (j == 1 && !cp.paired_end) continue;
      done[j] = false;
      uint32_t num_reads_read = read_fastq_block(
          fin[j], id_array, read_array, quality_array, num_reads_per_step, fasta_flag);
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
        if (tid * num_reads_per_block >= num_reads_read) done = true;
        uint32_t num_reads_thr = std::min((uint64_t)num_reads_read,
                                          (tid + 1) * num_reads_per_block) -
                                 tid * num_reads_per_block;
        std::ofstream fout_readlength;
        if (!done) {
          // check if reads and qualities have equal lengths
          for (uint32_t i = tid * num_reads_per_block;
               i < tid * num_reads_per_block + num_reads_thr; i++) {
            size_t len = read_array[i].size();
            if (len > MAX_READ_LEN) {
              std::cerr << "Max read length is "
                        << MAX_READ_LEN << ", but found read of length " << len
                        << "\n";
              throw std::runtime_error(
                  "Too long read length.");
            }
            if (!fasta_flag && quality_array[i].size() != len)
              throw std::runtime_error(
                  "Read length does not match quality length.");
            read_lengths_array[i] = (uint32_t)len;

            // mark reads with N
            read_contains_N_array[i] =
                (read_array[i].find('N') != std::string::npos);
          }
        }  // if(!done)
      }    // omp parallel
      // write reads and read_order_N to respective files
      for (uint32_t i = 0; i < num_reads_read; i++) {
        if (!read_contains_N_array[i]) {
          write_dna_in_bits(read_array[i],fout_clean[j]);
          num_reads_clean[j]++;
        } else {
          uint32_t pos_N = num_reads[j] + i;
          fout_order_N[j].write((char *)&pos_N, sizeof(uint32_t));
          write_dnaN_in_bits(read_array[i],fout_N[j]);
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
  delete[] id_array;
  delete[] read_contains_N_array;
  delete[] read_lengths_array;
  if (gzip_flag) {
    for (int j = 0; j < 2; j++) {
      if (j == 1 && !cp.paired_end) continue;
      delete fin[j];
      delete inbuf[j];
    }
  }
  for (int j = 0; j < 2; j++) {
    if (j == 1 && !cp.paired_end) continue;
    fin_f[j].close();
    fout_clean[j].close();
    fout_N[j].close();
    fout_order_N[j].close();
  }

  if (num_reads[0] == 0) throw std::runtime_error("No reads found.");

  if (cp.paired_end) {
    // merge input_N and input_order_N for the two files
    std::ofstream fout_N(outfileN[0], std::ios::app|std::ios::binary);
    std::ifstream fin_N(outfileN[1], std::ios::binary);
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

  cp.num_reads = num_reads[0] + num_reads[1];
  cp.num_reads_clean[0] = num_reads_clean[0];
  cp.num_reads_clean[1] = num_reads_clean[1];
  cp.max_readlen = max_readlen;

  std::cout << "Max Read length: " << cp.max_readlen << "\n";
  std::cout << "Total number of reads: " << cp.num_reads << "\n";

  std::cout << "Total number of reads without N: "
            << cp.num_reads_clean[0] + cp.num_reads_clean[1] << "\n";
}

}  // namespace spring
