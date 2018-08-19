#include "encoder.h"
#include <omp.h>
#include <algorithm>
#include <array>
#include <bitset>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <list>
#include <numeric>
#include <string>
#include <string>
#include <vector>
#include "bcm/bcm.h"
namespace spring {

std::string buildcontig(std::list<contig_reads> &current_contig,
                        const uint32_t &list_size) {
  static const char longtochar[5] = {'A', 'C', 'G', 'T', 'N'};
  static const long chartolong[128] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 4, 0,
    0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  };
  if (list_size == 1) return (current_contig.front()).read;
  auto current_contig_it = current_contig.begin();
  int64_t currentpos = 0, currentsize = 0, to_insert;
  std::vector<std::array<long, 4>> count;
  for (; current_contig_it != current_contig.end(); ++current_contig_it) {
    if (current_contig_it == current_contig.begin())  // first read
      to_insert = (*current_contig_it).read_length;
    else {
      currentpos = (*current_contig_it).pos;
      if (currentpos + (*current_contig_it).read_length > currentsize)
        to_insert = currentpos + (*current_contig_it).read_length - currentsize;
      else
        to_insert = 0;
    }
    count.insert(count.end(), to_insert, {0, 0, 0, 0});
    currentsize = currentsize + to_insert;
    for (long i = 0; i < (*current_contig_it).read_length; i++)
      count[currentpos + i]
           [chartolong[(uint8_t)(*current_contig_it).read[i]]] += 1;
  }
  std::string ref(count.size(), 'A');
  for (size_t i = 0; i < count.size(); i++) {
    long max = 0, indmax = 0;
    for (long j = 0; j < 4; j++)
      if (count[i][j] > max) {
        max = count[i][j];
        indmax = j;
      }
    ref[i] = longtochar[indmax];
  }
  return ref;
}

void writecontig(const std::string &ref, std::list<contig_reads> &current_contig,
                 std::ofstream &f_seq, std::ofstream &f_pos,
                 std::ofstream &f_noise, std::ofstream &f_noisepos,
                 std::ofstream &f_order, std::ofstream &f_RC,
                 std::ofstream &f_readlength, const encoder_global &eg, uint64_t &abs_pos) {
  f_seq << ref;
  uint16_t pos_var;
  long prevj = 0;
  auto current_contig_it = current_contig.begin();
  long currentpos;
  uint64_t abs_current_pos;
  for (; current_contig_it != current_contig.end(); ++current_contig_it) {
    currentpos = (*current_contig_it).pos;
    prevj = 0;
    for (long j = 0; j < (*current_contig_it).read_length; j++)
      if ((*current_contig_it).read[j] != ref[currentpos + j]) {
        f_noise << eg.enc_noise[(uint8_t)ref[currentpos + j]]
                               [(uint8_t)(*current_contig_it).read[j]];
        pos_var = j - prevj;
	      f_noisepos.write((char*)&pos_var,sizeof(uint16_t));
        prevj = j;
      }
    f_noise << "\n";
    abs_current_pos = abs_pos + currentpos;
    f_pos.write((char*)&abs_current_pos,sizeof(uint64_t));
    f_order.write((char *)&((*current_contig_it).order), sizeof(uint32_t));
    f_readlength.write((char *)&((*current_contig_it).read_length),
                       sizeof(uint16_t));
    f_RC << (*current_contig_it).RC;
  }
  abs_pos += ref.size();
  return;
}

void pack_compress_seq(const encoder_global &eg, uint64_t *file_len_seq_thr) {
#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    // seq
    std::ifstream in_seq(eg.outfile_seq + '.' + std::to_string(tid));
    std::ofstream f_seq(eg.outfile_seq + '.' + std::to_string(tid) + ".tmp",
                        std::ios::binary);
    std::ofstream f_seq_tail(eg.outfile_seq + '.' + std::to_string(tid) +
                             ".tail");
    uint64_t file_len = 0;
    char c;
    while (in_seq >> std::noskipws >> c) file_len++;
    file_len_seq_thr[tid] = file_len;
    uint8_t basetoint[128];
    basetoint[(uint8_t)'A'] = 0;
    basetoint[(uint8_t)'C'] = 1;
    basetoint[(uint8_t)'G'] = 2;
    basetoint[(uint8_t)'T'] = 3;

    in_seq.close();
    in_seq.open(eg.outfile_seq + '.' + std::to_string(tid));
    char dnabase[8];
    uint8_t dnabin;
    for (uint64_t i = 0; i < file_len / 4; i++) {
      in_seq.read(dnabase, 4);

      dnabin = 64 * basetoint[(uint8_t)dnabase[3]] +
               16 * basetoint[(uint8_t)dnabase[2]] +
               4 * basetoint[(uint8_t)dnabase[1]] +
               basetoint[(uint8_t)dnabase[0]];
      f_seq.write((char *)&dnabin, sizeof(uint8_t));
    }
    f_seq.close();
    in_seq.read(dnabase, file_len % 4);
    for (unsigned int i = 0; i < file_len % 4; i++) f_seq_tail << dnabase[i];
    f_seq_tail.close();
    in_seq.close();
    bcm::bcm_compress((eg.outfile_seq + '.' + std::to_string(tid) + ".tmp").c_str(),
                      (eg.outfile_seq + '.' + std::to_string(tid) + ".bcm").c_str());
    remove((eg.outfile_seq + '.' + std::to_string(tid)).c_str());
    remove((eg.outfile_seq + '.' + std::to_string(tid) + ".tmp").c_str());
  }

  return;
}

void getDataParams(encoder_global &eg, const compression_params &cp) {
  uint32_t numreads_clean, numreads_total;
  numreads_clean = cp.num_reads_clean[0] + cp.num_reads_clean[1];
  numreads_total = cp.num_reads;

  std::ifstream myfile_s(eg.infile + ".singleton", std::ifstream::in);
  eg.numreads_s = 0;
  std::string line;
  while (std::getline(myfile_s, line)) ++eg.numreads_s;
  myfile_s.close();
  eg.numreads = numreads_clean - eg.numreads_s;
  eg.numreads_N = numreads_total - numreads_clean;

  std::cout << "Maximum Read length: " << eg.max_readlen << std::endl;
  std::cout << "Number of non-singleton reads: " << eg.numreads << std::endl;
  std::cout << "Number of singleton reads: " << eg.numreads_s << std::endl;
  std::cout << "Number of reads with N: " << eg.numreads_N << std::endl;
}

void correct_order(uint32_t *order_s, const encoder_global &eg) {
  uint32_t numreads_total = eg.numreads + eg.numreads_s + eg.numreads_N;
  bool *read_flag_N = new bool[numreads_total]();
  // bool array indicating N reads
  for (uint32_t i = 0; i < eg.numreads_N; i++) {
    read_flag_N[order_s[eg.numreads_s + i]] = true;
  }

  uint32_t *cumulative_N_reads = new uint32_t[eg.numreads + eg.numreads_s];
  // number of reads occuring before pos in clean reads
  uint32_t pos_in_clean = 0, num_N_reads_till_now = 0;
  for (uint32_t i = 0; i < numreads_total; i++) {
    if (read_flag_N[i] == true)
      num_N_reads_till_now++;
    else
      cumulative_N_reads[pos_in_clean++] = num_N_reads_till_now;
  }

  // First correct the order for singletons
  for (uint32_t i = 0; i < eg.numreads_s; i++)
    order_s[i] += cumulative_N_reads[order_s[i]];

  // Now correct for clean reads (this is stored on file)
  for (int tid = 0; tid < eg.num_thr; tid++) {
    std::ifstream fin_order(eg.infile_order + '.' + std::to_string(tid),
                            std::ios::binary);
    std::ofstream fout_order(
        eg.infile_order + '.' + std::to_string(tid) + ".tmp", std::ios::binary);
    uint32_t pos;
    fin_order.read((char *)&pos, sizeof(uint32_t));
    while (!fin_order.eof()) {
      pos += cumulative_N_reads[pos];
      fout_order.write((char *)&pos, sizeof(uint32_t));
      fin_order.read((char *)&pos, sizeof(uint32_t));
    }
    fin_order.close();
    fout_order.close();
    remove((eg.infile_order + '.' + std::to_string(tid)).c_str());
    rename((eg.infile_order + '.' + std::to_string(tid) + ".tmp").c_str(),
           (eg.infile_order + '.' + std::to_string(tid)).c_str());
  }
  remove(eg.infile_order_N.c_str());
  delete[] read_flag_N;
  delete[] cumulative_N_reads;
  return;
}

}  // namespace spring
