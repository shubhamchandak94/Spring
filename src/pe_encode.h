#ifndef SPRING_PE_ENCODE_H_
#define SPRING_PE_ENCODE_H_

#include <string>
#include "algorithms/SPRING/pe_encode.h"

namespace spring {

struct pe_encode_global {
  std::string infile_order;
  std::string infilenumreads;

  std::string outfile_order_paired;
  std::string outfile_paired_flag_first;

  uint32_t numreads, numreads_by_2;
};

void populate_arrays(uint32_t* read_order, uint32_t* read_inverse_order,
                     pe_encode_global& peg);
// populate arrays:
// read_order = pos in reordered file to pos in original file
// read_inverse_order = pos in original file to pos in reordered file

void write_order_paired(uint32_t* read_order, uint32_t* read_inverse_order,
                        pe_encode_global& peg);
// write to all output files
// order_paired - store relative position of paired read once per pair (store
// for the read occuring first in the reordered file)
// For each pair, paired_flag_first stores 1 is 1st read comes first.

void packbits(pe_encode_global& peg);
// pack flag files into 1 bit per flag

void generate_order_preserve(uint32_t* read_order, pe_encode_global& peg);
// generate order file for half the reads

void pe_encode_main(std::string& working_dir, bool preserve_order);

}  // namespace spring

#endif  // SPRING_PE_ENCODE_H_
