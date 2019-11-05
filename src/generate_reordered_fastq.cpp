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

#include <stdexcept>
#include <fstream>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>

#include "generate_reordered_fastq.h"

namespace spring {

void generate_reordered_fastq(const std::string &temp_dir,
                              const compression_params &cp,
                              const std::vector<std::string> &infile_vector,
                              const std::vector<std::string> &outfile_vector,
                              const bool gzipped_output_flag,
                              const bool gzipped_input_flag) {
  uint32_t num_pairs = cp.paired_end?cp.num_reads/2:cp.num_reads;
  uint8_t num_files = cp.paired_end?2:1;
  uint32_t *order_array = new uint32_t[num_pairs];
  std::string file_order = temp_dir + "/read_order.bin";
  if (cp.paired_end)
    generate_order_pe(file_order, order_array, cp.num_reads);
  else
    generate_order_se(file_order, order_array, cp.num_reads);
  if (!is_permutation(order_array, num_pairs)) throw std::runtime_error("order_array not permutation of 1...num_pairs.");

  uint32_t str_array_size = cp.paired_end ? num_pairs / 5 + 1 ? num_pairs / 10 + 1;
  // chosen so that roughly these many FASTQ records can be stored in
  // memory without exceeding the RAM consumption of reordering stage

  std::string *read_array = new std::string[str_array_size];
  std::string *id_array = new std::string[str_array_size];
  std::string *quality_array = new std::string[str_array_size];
  std::string cur_read, cur_id, cur_quality, cur_comment;

  std::ofstream fout_fastq;
  std::ifstream fin_f[2];
  std::istream *fin[2] = {&fin_f[0], &fin_f[1]};
  boost::iostreams::filtering_streambuf<boost::iostreams::input> *inbuf[2];
  for (int j = 0; j < 2; j++) {
    if (j == 1 && !cp.paired_end) continue;
    if (gzipped_input_flag) {
      fin_f[j].open(infile_vector[j], std::ios_base::binary);
      inbuf[j] =
          new boost::iostreams::filtering_streambuf<boost::iostreams::input>;
      inbuf[j]->push(boost::iostreams::gzip_decompressor());
      inbuf[j]->push(fin_f[j]);
      fin[j] = new std::istream(inbuf[j]);
    } else {
      fin_f[j].open(infile_vector[j]);
    }
  }

  for (uint8_t file_num = 0; file_num < num_files; file_num++) {
    if (gzipped_output_flag)
      fout_fastq.open(outfile_vector[file_num], std::ios::binary);
    else
      fout_fastq.open(outfile_vector[file_num]);
    if (!fout_fastq.is_open()) throw std::runtime_error("Error opening output file");

    for (uint32_t i = 0; i <= num_pairs / str_array_size; i++) {
      uint32_t num_reads_bin = str_array_size;
      if (i == num_pairs / str_array_size) num_reads_bin = num_pairs % str_array_size;
      if (num_reads_bin == 0) break;
      uint32_t start_read_bin = i * str_array_size;
      uint32_t end_read_bin = i * str_array_size + num_reads_bin;
      // Read the file and pick up lines corresponding to this bin

      // file seek to 0
      if (gzipped_input_flag) {
        fin_f[file_num].close();
        delete fin[file_num];
        delete inbuf[file_num];
        fin_f[file_num].open(infile_vector[file_num], std::ios_base::binary);
        inbuf[file_num] = new boost::iostreams::filtering_streambuf<
            boost::iostreams::input>;
        inbuf[file_num]->push(boost::iostreams::gzip_decompressor());
        inbuf[file_num]->push(fin_f[file_num]);
        fin[file_num] = new std::istream(inbuf[file_num]);
      }
      else {
        fin_f[file_num].seekg(0);
      }

      for (uint32_t j = 0; j < num_pairs; j++) {
          std::getline(*fin[file_num], cur_id);
          std::getline(*fin[file_num], cur_read);
          std::getline(*fin[file_num], cur_comment);
          std::getline(*fin[file_num], cur_quality);
          if (order_array[j] >= start_read_bin && order_array[j] < end_read_bin) {
              uint32_t index = order_array[j] - start_read_bin;
              read_array[index] = cur_read;
              id_array[index] = cur_id;
              quality_array[index] = cur_quality;
          }
      }

      write_fastq_block(fout_fastq, id_array, read_array, quality_array,
                        num_reads_bin, true, cp.num_thr, gzipped_output_flag);
    }
    fout_fastq.close();
  }

  delete[] read_array;
  delete[] id_array;
  delete[] quality_array;
  delete[] order_array;
  if (gzipped_input_flag) {
    for (int j = 0; j < 2; j++) {
      if (j == 1 && !cp.paired_end) continue;
      delete fin[j];
      delete inbuf[j];
    }
  }
  for (int j = 0; j < 2; j++) {
    if (j == 1 && !cp.paired_end) continue;
    fin_f[j].close();
  }
}

void generate_order_pe(const std::string &file_order, uint32_t *order_array,
                       const uint32_t &numreads) {
  std::ifstream fin_order(file_order, std::ios::binary);
  uint32_t order;
  uint32_t pos_after_reordering = 0;
  uint32_t numreads_by_2 = numreads / 2;
  for (uint32_t i = 0; i < numreads; i++) {
    fin_order.read((char *)&order, sizeof(uint32_t));
    if (order < numreads_by_2) {
      order_array[order] = pos_after_reordering++;
    }
  }
  fin_order.close();
}

void generate_order_se(const std::string &file_order, uint32_t *order_array,
                       const uint32_t &numreads) {
  std::ifstream fin_order(file_order, std::ios::binary);
  uint32_t order;
  for (uint32_t i = 0; i < numreads; i++) {
    fin_order.read((char *)&order, sizeof(uint32_t));
    order_array[order] = i;
  }
  fin_order.close();
}

} // namespace spring
