#include <omp.h>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

#include "util.h"
#include "libbsc/bsc.h"
#include "ID_compression/include/sam_block.h"
#include "reorder_compress_quality_id.h"

namespace spring {

void reorder_compress_quality_id(const std::string &temp_dir, const compression_params &cp) {
  // Read some parameters
  uint32_t numreads = cp.num_reads;
  int num_thr = cp.num_thr;
  bool preserve_id = cp.preserve_id;
  bool preserve_quality = cp.preserve_quality;
  bool paired_end = cp.paired_end;
  uint32_t num_reads_per_block = cp.num_reads_per_block;
  bool paired_id_match = cp.paired_id_match;

  std::string basedir = temp_dir;

  std::string file_order = basedir + "/read_order.bin";
  std::string file_id[2];
  std::string file_quality[2];
  file_id[0] = basedir + "/id_1";
  file_id[1] = basedir + "/id_2";
  file_quality[0] = basedir + "/quality_1";
  file_quality[1] = basedir + "/quality_2";

  uint32_t *order_array;
  // array containing index mapping position in original fastq to
  // position after reordering
  if(paired_end) {
    order_array = new uint32_t[numreads/2];
    generate_order_pe(file_order, order_array, numreads);
  }
  else {
    order_array = new uint32_t[numreads];
    generate_order_se(file_order, order_array, numreads);
  }

  omp_set_num_threads(num_thr);

  uint32_t str_array_size = (1 + (numreads/4 - 1)/num_reads_per_block)*num_reads_per_block;
  // smallest multiple of num_reads_per_block bigger than numreads/4
  // numreads/4 chosen so that these many qualities/ids can be stored in
  // memory without exceeding the RAM consumption of reordering stage
  std::string *str_array = new std::string[str_array_size];
  // array to load ids and/or qualities into

  if(preserve_quality) {
    std::cout << "Compressing qualities\n";
    for(int j = 0; j < 2; j++) {
      if(!paired_end && j==1)
        break;
      uint32_t num_reads_per_file = paired_end?numreads/2:numreads;
      reorder_compress(file_quality[j], num_reads_per_file, num_thr, num_reads_per_block, str_array, str_array_size, order_array, "quality");
      remove(file_quality[j].c_str());
    }
  }
  if(preserve_id) {
    std::cout << "Compressing ids\n";
    for(int j = 0; j < 2; j++) {
      if(!paired_end && j==1)
        break;
      if(j == 1 && paired_id_match)
        break;
      uint32_t num_reads_per_file = paired_end?numreads/2:numreads;
      reorder_compress(file_id[j], num_reads_per_file, num_thr, num_reads_per_block, str_array, str_array_size, order_array, "id");
      remove(file_id[j].c_str());
    }
  }

  delete[] order_array;
  delete[] str_array;
  return;
}

void generate_order_pe(const std::string &file_order, uint32_t *order_array, const uint32_t &numreads) {
  std::ifstream fin_order(file_order, std::ios::binary);
  uint32_t order;
  uint32_t pos_after_reordering = 0;
  uint32_t numreads_by_2 = numreads/2;
  for (uint32_t i = 0; i < numreads; i++) {
    fin_order.read((char *)&order, sizeof(uint32_t));
    if (order < numreads_by_2) {
      order_array[order] = pos_after_reordering++;
    }
  }
  fin_order.close();
}

void generate_order_se(const std::string &file_order, uint32_t *order_array, const uint32_t &numreads) {
  std::ifstream fin_order(file_order, std::ios::binary);
  uint32_t order;
  for (uint32_t i = 0; i < numreads; i++) {
    fin_order.read((char *)&order, sizeof(uint32_t));
    order_array[order] = i;
  }
  fin_order.close();
}

void reorder_compress(const std::string &file_name, const uint32_t &num_reads_per_file, const int &num_thr, const uint32_t &num_reads_per_block, std::string *str_array, const uint32_t &str_array_size, uint32_t *order_array, const std::string &mode) {
  uint32_t *read_lengths_array;
  // Allocate array of lengths if mode is "quality"
  if(mode == "quality") {
    uint64_t num_reads_per_step = (uint64_t)num_thr*num_reads_per_block;
    // allocate less if str_array_size is very small
    if(num_reads_per_step > str_array_size)
      num_reads_per_step = str_array_size;
    read_lengths_array = new uint32_t [num_reads_per_step];
  }

  for(uint32_t i = 0; i <= num_reads_per_file/str_array_size; i++) {
    uint32_t num_reads_bin = str_array_size;
    if(i == num_reads_per_file/str_array_size)
      num_reads_bin = num_reads_per_file%str_array_size;
    if(num_reads_bin == 0)
      break;
    uint32_t start_read_bin = i*str_array_size;
    uint32_t end_read_bin = i*str_array_size + num_reads_bin;
    // Read the file and pick up lines corresponding to this bin
    std::ifstream f_in(file_name);
    std::string temp_str;
    for(uint32_t i = 0; i < num_reads_per_file; i++) {
      std::getline(f_in, temp_str);
      if(order_array[i] >= start_read_bin && order_array[i] < end_read_bin)
        str_array[order_array[i] - start_read_bin] = temp_str;
    }
    f_in.close();
    #pragma omp parallel
    {
      uint64_t tid = omp_get_thread_num();
      uint64_t block_num_offset = start_read_bin/num_reads_per_block;
      uint64_t block_num = tid;
      bool done = false;
      while(!done) {
        uint64_t start_read_num = block_num*num_reads_per_block;
        uint64_t end_read_num = (block_num + 1)*num_reads_per_block;
        if(start_read_num >= num_reads_bin)
          break;
        if(end_read_num >= num_reads_bin) {
          done = true;
          end_read_num = num_reads_bin;
        }
        uint32_t num_reads_block = (uin32_t)(end_read_num - start_read_num);
        std::string outfile_name = file_name + "." + std::to_string(block_num_offset+block_num);

        if(mode == "id") {
          compress_id_block(outfile_name.c_str(), str_array + start_read_num, num_reads_block);
        }
        else {
          // store lengths in array for quality compression
          for(uint64_t i = start_read_num; i < end_read_num; i++)
            read_lengths_array[i] = str_array[i].size();
          bsc::BSC_str_array_compress(outfile_name.c_str(), str_array + start_read_num, num_reads_block, read_lengths_array + start_read_num);
        }
        block_num += num_thr;
      }
    } // omp parallel
  }
  if(mode == "quality")
    delete[] read_lengths_array;
}

} // namespace spring
