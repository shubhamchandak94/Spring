#include "pe_encode.h"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>

namespace spring {

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

void pe_encode_main(std::string& temp_dir, bool preserve_order) {
  pe_encode_global* peg_ptr = new pe_encode_global;
  pe_encode_global& peg = *peg_ptr;
  std::string basedir = temp_dir;
  peg.infilenumreads = basedir + "/numreads.bin";
  peg.infile_order = basedir + "/read_order.bin";
  peg.outfile_order_paired = basedir + "/read_order_paired.bin";
  peg.outfile_paired_flag_first = basedir + "/read_paired_flag_first.bin";

  std::ifstream f_numreads(peg.infilenumreads, std::ios::binary);
  f_numreads.seekg(4);
  f_numreads.read((char*)&peg.numreads, sizeof(uint32_t));
  f_numreads.close();
  peg.numreads_by_2 = peg.numreads / 2;

  uint32_t* read_order = new uint32_t[peg.numreads];
  uint32_t* read_inverse_order = new uint32_t[peg.numreads];
  populate_arrays(read_order, read_inverse_order, peg);
  write_order_paired(read_order, read_inverse_order, peg);
  packbits(peg);

  if (preserve_order == true) generate_order_preserve(read_order, peg);

  delete[] read_order;
  delete[] read_inverse_order;
  delete peg_ptr;
}

void populate_arrays(uint32_t* read_order, uint32_t* read_inverse_order,
                     pe_encode_global& peg) {
  // read file read_order
  std::ifstream f_order(peg.infile_order, std::ios::binary);
  for (uint32_t i = 0; i < peg.numreads; i++) {
    f_order.read((char*)&read_order[i], sizeof(uint32_t));
  }
  f_order.close();

  // now fill read_inverse_order
  for (uint32_t i = 0; i < peg.numreads; i++) {
    read_inverse_order[read_order[i]] = i;
  }
  return;
}

void write_order_paired(uint32_t* read_order, uint32_t* read_inverse_order,
                        pe_encode_global& peg) {
  std::ofstream f_flag_first(peg.outfile_paired_flag_first);
  std::ofstream f_order_paired(peg.outfile_order_paired, std::ios::binary);
  for (uint32_t i = 0; i < peg.numreads; i++) {
    if (read_order[i] < peg.numreads_by_2)  // first read of pair
    {
      if (read_inverse_order[read_order[i] + peg.numreads_by_2] >
          i)  // pair not already seen
      {
        uint32_t temp =
            (read_inverse_order[read_order[i] + peg.numreads_by_2] - i);
        f_order_paired.write((char*)&temp, sizeof(uint32_t));
        f_flag_first << '1';
      }
    } else {
      if (read_inverse_order[read_order[i] - peg.numreads_by_2] >
          i)  // pair not already seen
      {
        f_flag_first << '0';
        uint32_t temp =
            (read_inverse_order[read_order[i] - peg.numreads_by_2] - i);
        f_order_paired.write((char*)&temp, sizeof(uint32_t));
      }
    }
  }
  f_flag_first.close();
  f_order_paired.close();
}

void packbits(pe_encode_global& peg) {
  // flag_first
  std::ifstream in_flag_first(peg.outfile_paired_flag_first);
  std::ofstream f_flag_first(peg.outfile_paired_flag_first + ".tmp",
                             std::ios::binary);
  std::ofstream f_flag_first_tail(peg.outfile_paired_flag_first + ".tail");

  uint8_t chartoint[128];
  chartoint[(uint8_t)'0'] = 0;
  chartoint[(uint8_t)'1'] = 1;
  in_flag_first.close();
  in_flag_first.open(peg.outfile_paired_flag_first);
  char chararray[8];
  uint8_t packedchar;
  for (uint64_t i = 0; i < peg.numreads_by_2 / 8; i++) {
    in_flag_first.read(chararray, 8);

    packedchar = 128 * chartoint[(uint8_t)chararray[7]] +
                 64 * chartoint[(uint8_t)chararray[6]] +
                 32 * chartoint[(uint8_t)chararray[5]] +
                 16 * chartoint[(uint8_t)chararray[4]] +
                 8 * chartoint[(uint8_t)chararray[3]] +
                 4 * chartoint[(uint8_t)chararray[2]] +
                 2 * chartoint[(uint8_t)chararray[1]] +
                 1 * chartoint[(uint8_t)chararray[0]];
    f_flag_first.write((char*)&packedchar, sizeof(uint8_t));
  }
  f_flag_first.close();
  in_flag_first.read(chararray, peg.numreads_by_2 % 8);
  for (unsigned int i = 0; i < peg.numreads_by_2 % 8; i++)
    f_flag_first_tail << chararray[i];
  f_flag_first_tail.close();
  in_flag_first.close();
  remove((peg.outfile_paired_flag_first).c_str());
  rename((peg.outfile_paired_flag_first + ".tmp").c_str(),
         (peg.outfile_paired_flag_first).c_str());
  return;
}

void generate_order_preserve(uint32_t* read_order, pe_encode_global& peg) {
  std::ofstream fout_order(peg.infile_order, std::ios::binary);
  for (uint32_t i = 0; i < peg.numreads; i++) {
    if (read_order[i] < peg.numreads_by_2) {
      fout_order.write((char*)&read_order[i], sizeof(uint32_t));
    }
  }
  fout_order.close();
  return;
}

}  // namespace spring
