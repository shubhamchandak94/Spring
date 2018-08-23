#include <cstdio>
#include <fstream>
#include <iostream>

#include "util.h"
#include "pe_encode.h"

namespace spring {

void pe_encode(const std::string &temp_dir, const compression_params &cp) {
  // Read some parameters
  uint32_t numreads = cp.num_reads;
  uint32_t num_reads_by_2 = numreads/2;

  std::string basedir = temp_dir;
  std::string file_order = basedir + "/read_order.bin";
  uint32_t *order_array = new uint32_t [numreads];
  // stores index mapping position in reordered file to
  // position in original file
  // later stores index generated in this step - maps position
  // in reordered file to position in intended decompressed file

  uint32_t *inverse_order_array = new uint32_t [numreads];
  // stores index mapping position in original file to
  // position in reordered file

  std::ifstream fin_order(file_order, std::ios::binary);
  uint32_t order;
  for (uint32_t i = 0; i < numreads; i++) {
    fin_order.read((char *)&order, sizeof(uint32_t));
    order_array[i] = order;
    inverse_order_array[order] = i;
  }
  fin_order.close();

  // First fill positions in new order array corresponding to reads
  // in file 1. These reads are decompressed in the same order as
  // in the reordered file
  uint32_t pos_in_file_1 = 0;
  for (uint32_t i = 0; i < numreads; i++) {
    if(order_array[i] < numreads_by_2)
      order_array[i] = pos_in_file_1++;
  }

  // Now fill positions in new order array corresponding to reads in
  // file 2. These are automatically decided by the pairing.
  for (uint32_t i = 0; i < numreads; i++) {
    if(order_array[i] >=numreads_by_2) {
      uint32_t pos_in_original = order_array[i];
      uint32_t pos_of_pair_in_original = pos_in_original - numreads_by_2;
      uint32_t pos_of_pair_in_reordered = inverse_order_array[pos_of_pair_in_original];
      uint32_t new_order_of_pair = order_array[pos_of_pair_in_reordered];
      order_array[i] = new_order_of_pair + numreads_by_2;
    }
  }

  // Write to tmp file and replace
  std::ifstream fout_order(file_order+".tmp", std::ios::binary);
  for (uint32_t i = 0; i < numreads; i++) {
    fout_order.write((char *)&order_array[i], sizeof(uint32_t));
  }
  fout_order.close();

  remove(file_order.c_str());
  rename((file_order + ".tmp").c_str(), file_order.c_str());

  delete[] order_array;
  delete[] inverse_order_array;
  return;
}

} // namespace spring
