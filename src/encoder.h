#ifndef SPRING_ENCODER_H_
#define SPRING_ENCODER_H_

#include <omp.h>
#include <algorithm>
#include <bitset>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include "bitset_util.h"
#include "params.h"
#include "util.h"

namespace spring {

template <size_t bitset_size>
struct encoder_global_b {
  std::bitset<bitset_size> **basemask;
  int max_readlen;
  // bitset for A,G,C,T,N at each position
  // used in stringtobitset, and bitsettostring
  std::bitset<bitset_size>
      mask63;  // bitset with 63 bits set to 1 (used in bitsettostring for
               // conversion to ullong)
  encoder_global_b(int max_readlen_param) {
    max_readlen = max_readlen_param;
    basemask = new std::bitset<bitset_size> *[max_readlen_param];
    for (int i = 0; i < max_readlen_param; i++)
      basemask[i] = new std::bitset<bitset_size>[128];
  }
  ~encoder_global_b() {
    for (int i = 0; i < max_readlen; i++) delete[] basemask[i];
    delete[] basemask;
  }
};

struct encoder_global {
  uint32_t numreads, numreads_s, numreads_N;
  int numdict_s = NUM_DICT_ENCODER;

  int max_readlen, num_thr;

  std::string basedir;
  std::string infile;
  std::string infile_flag;
  std::string infile_pos;
  std::string infile_seq;
  std::string infile_RC;
  std::string infile_readlength;
  std::string infile_N;
  std::string outfile_unaligned;
  std::string outfile_seq;
  std::string outfile_pos;
  std::string outfile_noise;
  std::string outfile_noisepos;
  std::string infile_order;
  std::string infile_order_N;

  char enc_noise[128][128];
};

struct contig_reads {
  std::string read;
  int64_t pos;
  char RC;
  uint32_t order;
  uint16_t read_length;
};

std::string buildcontig(std::list<contig_reads> &current_contig,
                        const uint32_t &list_size);

void writecontig(const std::string &ref,
                 std::list<contig_reads> &current_contig, std::ofstream &f_seq,
                 std::ofstream &f_pos, std::ofstream &f_noise,
                 std::ofstream &f_noisepos, std::ofstream &f_order,
                 std::ofstream &f_RC, std::ofstream &f_readlength,
                 const encoder_global &eg, uint64_t &abs_pos);

void pack_compress_seq(const encoder_global &eg, uint64_t *file_len_seq_thr);

void getDataParams(encoder_global &eg, const compression_params &cp);

void correct_order(uint32_t *order_s, const encoder_global &eg);

template <size_t bitset_size>
std::string bitsettostring(std::bitset<bitset_size> b, const uint16_t readlen,
                           const encoder_global_b<bitset_size> &egb) {
  // destroys bitset b
  static const char revinttochar[8] = {'A', 'N', 'G', 0, 'C', 0, 'T', 0};
  std::string s;
  s.resize(readlen);
  unsigned long long ull;
  for (int i = 0; i < 3 * readlen / 63 + 1; i++) {
    ull = (b & egb.mask63).to_ullong();
    b >>= 63;
    for (int j = 21 * i; j < 21 * i + 21 && j < readlen; j++) {
      s[j] = revinttochar[ull % 8];
      ull /= 8;
    }
  }
  return s;
}

template <size_t bitset_size>
void encode(std::bitset<bitset_size> *read, bbhashdict *dict, uint32_t *order_s,
            uint16_t *read_lengths_s, const encoder_global &eg,
            const encoder_global_b<bitset_size> &egb) {
  static const int thresh_s = THRESH_ENCODER;
  static const int maxsearch = MAX_SEARCH_ENCODER;
  omp_lock_t *read_lock = new omp_lock_t[eg.numreads_s + eg.numreads_N];
  omp_lock_t *dict_lock = new omp_lock_t[eg.numreads_s + eg.numreads_N];
  for (uint64_t j = 0; j < eg.numreads_s + eg.numreads_N; j++) {
    omp_init_lock(&read_lock[j]);
    omp_init_lock(&dict_lock[j]);
  }
  bool *remainingreads = new bool[eg.numreads_s + eg.numreads_N];
  std::fill(remainingreads, remainingreads + eg.numreads_s + eg.numreads_N, 1);

  std::bitset<bitset_size> *mask1 = new std::bitset<bitset_size>[eg.numdict_s];
  generateindexmasks<bitset_size>(mask1, dict, eg.numdict_s, 3);
  std::bitset<bitset_size> **mask =
      new std::bitset<bitset_size> *[eg.max_readlen];
  for (int i = 0; i < eg.max_readlen; i++)
    mask[i] = new std::bitset<bitset_size>[eg.max_readlen];
  generatemasks<bitset_size>(mask, eg.max_readlen, 3);
  std::cout << "Encoding reads\n";
#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    std::ifstream f(eg.infile + '.' + std::to_string(tid));
    std::ifstream in_flag(eg.infile_flag + '.' + std::to_string(tid));
    std::ifstream in_pos(eg.infile_pos + '.' + std::to_string(tid),
                         std::ios::binary);
    std::ifstream in_order(eg.infile_order + '.' + std::to_string(tid),
                           std::ios::binary);
    std::ifstream in_RC(eg.infile_RC + '.' + std::to_string(tid));
    std::ifstream in_readlength(
        eg.infile_readlength + '.' + std::to_string(tid), std::ios::binary);
    std::ofstream f_seq(eg.outfile_seq + '.' + std::to_string(tid));
    std::ofstream f_pos(eg.outfile_pos + '.' + std::to_string(tid));
    std::ofstream f_noise(eg.outfile_noise + '.' + std::to_string(tid));
    std::ofstream f_noisepos(eg.outfile_noisepos + '.' + std::to_string(tid));
    std::ofstream f_order(eg.infile_order + '.' + std::to_string(tid) + ".tmp",
                          std::ios::binary);
    std::ofstream f_RC(eg.infile_RC + '.' + std::to_string(tid) + ".tmp");
    std::ofstream f_readlength(
        eg.infile_readlength + '.' + std::to_string(tid) + ".tmp",
        std::ios::binary);
    int64_t dictidx[2];  // to store the start and end index (end not inclusive)
                         // in the dict read_id array
    uint64_t startposidx;  // index in startpos
    uint64_t ull;
    uint64_t abs_pos = 0;  // absolute position in reference (total length
    // of all contigs till now)
    bool flag = 0;
    // flag to check if match was found or not
    std::string current, ref;
    std::bitset<bitset_size> forward_bitset, reverse_bitset, b;
    char c = '0', rc = 'd';
    std::list<contig_reads> current_contig;
    int64_t p;
    uint16_t rl;
    uint32_t ord, list_size = 0;  // list_size variable introduced because
                                  // list::size() was running very slowly
                                  // on UIUC machine
    std::list<uint32_t> *deleted_rids = new std::list<uint32_t>[eg.numdict_s];
    bool done = false;
    while (!done) {
      if (!(in_flag >> c)) done = true;
      if (!done) {
        std::getline(f, current);
        rc = in_RC.get();
        in_pos.read((char *)&p, sizeof(int64_t));
        in_order.read((char *)&ord, sizeof(uint32_t));
        in_readlength.read((char *)&rl, sizeof(uint16_t));
      }
      if (c == '0' || done || list_size > 10000000)  // limit on list size so
                                                     // that memory doesn't get
                                                     // too large
      {
        if (list_size != 0) {
          // sort contig according to pos
          current_contig.sort([](const contig_reads &a, const contig_reads &b) {
            return a.pos < b.pos;
          });
          // make first pos zero and shift all pos values accordingly
          auto current_contig_it = current_contig.begin();
          int64_t first_pos = (*current_contig_it).pos;
          for (; current_contig_it != current_contig.end(); ++current_contig_it)
            (*current_contig_it).pos -= first_pos;

          ref = buildcontig(current_contig, list_size);
          // try to align the singleton reads to ref
          // first create bitsets from first readlen positions of ref
          forward_bitset.reset();
          reverse_bitset.reset();
          if ((int64_t)ref.size() >= eg.max_readlen) {
            stringtobitset(ref.substr(0, eg.max_readlen), eg.max_readlen,
                           forward_bitset, egb.basemask);
            stringtobitset(reverse_complement(ref.substr(0, eg.max_readlen),
                                              eg.max_readlen),
                           eg.max_readlen, reverse_bitset, egb.basemask);
            for (long j = 0; j < (int64_t)ref.size() - eg.max_readlen + 1;
                 j++) {
              // search for singleton reads
              for (int rev = 0; rev < 2; rev++) {
                for (int l = 0; l < eg.numdict_s; l++) {
                  if (!rev)
                    b = forward_bitset & mask1[l];
                  else
                    b = reverse_bitset & mask1[l];
                  ull = (b >> 3 * dict[l].start).to_ullong();
                  startposidx = dict[l].bphf->lookup(ull);
                  if (startposidx >= dict[l].numkeys)  // not found
                    continue;
                  // check if any other thread is modifying same dictpos
                  if (!omp_test_lock(&dict_lock[startposidx])) continue;
                  dict[l].findpos(dictidx, startposidx);
                  if (dict[l].empty_bin[startposidx])  // bin is empty
                  {
                    omp_unset_lock(&dict_lock[startposidx]);
                    continue;
                  }
                  uint64_t ull1 =
                      ((read[dict[l].read_id[dictidx[0]]] & mask1[l]) >>
                       3 * dict[l].start)
                          .to_ullong();
                  if (ull ==
                      ull1)  // checking if ull is actually the key for this bin
                  {
                    for (int64_t i = dictidx[1] - 1;
                         i >= dictidx[0] && i >= dictidx[1] - maxsearch; i--) {
                      auto rid = dict[l].read_id[i];
                      int hamming;
                      if (!rev)
                        hamming =
                            ((forward_bitset ^ read[rid]) &
                             mask[0][eg.max_readlen - read_lengths_s[rid]])
                                .count();
                      else
                        hamming =
                            ((reverse_bitset ^ read[rid]) &
                             mask[0][eg.max_readlen - read_lengths_s[rid]])
                                .count();
                      if (hamming <= thresh_s) {
                        omp_set_lock(&read_lock[rid]);
                        if (remainingreads[rid]) {
                          remainingreads[rid] = 0;
                          flag = 1;
                        }
                        omp_unset_lock(&read_lock[rid]);
                      }
                      if (flag == 1)  // match found
                      {
                        flag = 0;
                        list_size++;
                        char rc = rev ? 'r' : 'd';
                        long pos =
                            rev ? (j + eg.max_readlen - read_lengths_s[rid])
                                : j;
                        std::string read_string =
                            rev ? reverse_complement(
                                      bitsettostring<bitset_size>(
                                          read[rid], read_lengths_s[rid], egb),
                                      read_lengths_s[rid])
                                : bitsettostring<bitset_size>(
                                      read[rid], read_lengths_s[rid], egb);
                        current_contig.push_back({read_string, pos, rc,
                                                  order_s[rid],
                                                  read_lengths_s[rid]});
                        for (int l1 = 0; l1 < eg.numdict_s; l1++) {
                          if (read_lengths_s[rid] > dict[l1].end)
                            deleted_rids[l1].push_back(rid);
                        }
                      }
                    }
                  }
                  omp_unset_lock(&dict_lock[startposidx]);
                  // delete from dictionaries
                  for (int l1 = 0; l1 < eg.numdict_s; l1++)
                    for (auto it = deleted_rids[l1].begin();
                         it != deleted_rids[l1].end();) {
                      b = read[*it] & mask1[l1];
                      ull = (b >> 3 * dict[l1].start).to_ullong();
                      startposidx = dict[l1].bphf->lookup(ull);
                      if (!omp_test_lock(&dict_lock[startposidx])) {
                        ++it;
                        continue;
                      }
                      dict[l1].findpos(dictidx, startposidx);
                      dict[l1].remove(dictidx, startposidx, *it);
                      it = deleted_rids[l1].erase(it);
                      omp_unset_lock(&dict_lock[startposidx]);
                    }
                }
              }
              if (j !=
                  (int64_t)ref.size() -
                      eg.max_readlen)  // not at last position,shift bitsets
              {
                forward_bitset >>= 3;
                forward_bitset = forward_bitset & mask[0][0];
                forward_bitset |=
                    egb.basemask[eg.max_readlen - 1]
                                [(uint8_t)ref[j + eg.max_readlen]];
                reverse_bitset <<= 3;
                reverse_bitset = reverse_bitset & mask[0][0];
                reverse_bitset |= egb.basemask[0][(
                    uint8_t)chartorevchar[(uint8_t)ref[j + eg.max_readlen]]];
              }

            }  // end for
          }    // end if
          // sort contig according to pos
          current_contig.sort([](const contig_reads &a, const contig_reads &b) {
            return a.pos < b.pos;
          });
          writecontig(ref, current_contig, f_seq, f_pos, f_noise, f_noisepos,
                      f_order, f_RC, f_readlength, eg, abs_pos);
        }
        if (!done) {
          current_contig = {{current, p, rc, ord, rl}};
          list_size = 1;
        }
      } else if (c == '1')  // read found during rightward search
      {
        current_contig.push_back({current, p, rc, ord, rl});
        list_size++;
      }
    }
    f.close();
    in_flag.close();
    in_pos.close();
    in_order.close();
    in_RC.close();
    in_readlength.close();
    f_seq.close();
    f_pos.close();
    f_noise.close();
    f_noisepos.close();
    f_order.close();
    f_readlength.close();
    f_RC.close();
    delete[] deleted_rids;
  } // end omp parallel

  // Combine files produced by the threads
  std::ofstream f_order(eg.infile_order);
  std::ofstream f_readlength(eg.infile_readlength);
  std::ofstream f_noisepos(eg.outfile_noisepos);
  std::ofstream f_noise(eg.outfile_noise);
  std::ofstream f_RC(eg.infile_RC);

  for (int tid = 0; tid < eg.num_thr; tid++) {
    std::ifstream in_order(eg.infile_order + '.' + std::to_string(tid) +
                           ".tmp");
    std::ifstream in_readlength(eg.infile_readlength + '.' +
                                std::to_string(tid) + ".tmp");
    std::ifstream in_RC(eg.infile_RC + '.' + std::to_string(tid) + ".tmp");
    std::ifstream in_noisepos(eg.outfile_noisepos + '.' + std::to_string(tid));
    std::ifstream in_noise(eg.outfile_noise + '.' + std::to_string(tid));
    f_order << in_order.rdbuf();
    f_order.clear();  // clear error flag in case in_order is empty
    f_noisepos << in_noisepos.rdbuf();
    f_noisepos.clear();  // clear error flag in case in_noisepos is empty
    f_noise << in_noise.rdbuf();
    f_noise.clear();  // clear error flag in case in_noise is empty
    f_readlength << in_readlength.rdbuf();
    f_readlength.clear();  // clear error flag in case in_readlength is empty
    f_RC << in_RC.rdbuf();
    f_RC.clear();  // clear error flag in case in_RC is empty

    remove((eg.infile_order + '.' + std::to_string(tid)).c_str());
    remove((eg.infile_order + '.' + std::to_string(tid) + ".tmp").c_str());
    remove((eg.infile_readlength + '.' + std::to_string(tid)).c_str());
    remove((eg.infile_readlength + '.' + std::to_string(tid) + ".tmp").c_str());
    remove((eg.outfile_noisepos + '.' + std::to_string(tid)).c_str());
    remove((eg.outfile_noise + '.' + std::to_string(tid)).c_str());
    remove((eg.infile_RC + '.' + std::to_string(tid) + ".tmp").c_str());
    remove((eg.infile_RC + '.' + std::to_string(tid)).c_str());
    remove((eg.infile_flag + '.' + std::to_string(tid)).c_str());
    remove((eg.infile_pos + '.' + std::to_string(tid)).c_str());
    remove((eg.infile + '.' + std::to_string(tid)).c_str());
  }
  f_order.close();
  f_readlength.close();
  // write remaining singleton reads now
  std::ofstream f_unaligned(eg.outfile_unaligned);
  f_order.open(eg.infile_order, std::ios::binary | std::ofstream::app);
  f_readlength.open(eg.infile_readlength,
                    std::ios::binary | std::ofstream::app);
  uint32_t matched_s = eg.numreads_s;
  for (uint32_t i = 0; i < eg.numreads_s; i++)
    if (remainingreads[i] == 1) {
      matched_s--;
      f_order.write((char *)&order_s[i], sizeof(uint32_t));
      f_readlength.write((char *)&read_lengths_s[i], sizeof(uint16_t));
      f_unaligned << bitsettostring<bitset_size>(read[i], read_lengths_s[i],
                                                 egb);
    }
  uint32_t matched_N = eg.numreads_N;
  for (uint32_t i = eg.numreads_s; i < eg.numreads_s + eg.numreads_N; i++)
    if (remainingreads[i] == 1) {
      matched_N--;
      f_unaligned << bitsettostring<bitset_size>(read[i], read_lengths_s[i],
                                                 egb);
      f_order.write((char *)&order_s[i], sizeof(uint32_t));
      f_readlength.write((char *)&read_lengths_s[i], sizeof(uint16_t));
    }
  f_order.close();
  f_readlength.close();
  f_unaligned.close();
  delete[] remainingreads;
  delete[] dict_lock;
  delete[] read_lock;
  for (int i = 0; i < eg.max_readlen; i++) delete[] mask[i];
  delete[] mask;
  delete[] mask1;

  // pack read_Seq and convert read_pos into 8 byte non-diff (absolute)
  // positions
  uint64_t *file_len_seq_thr = new uint64_t[eg.num_thr];
  uint64_t abs_pos = 0;
  uint64_t abs_pos_thr;
  pack_compress_seq(eg, file_len_seq_thr);
  std::ofstream fout_pos(eg.outfile_pos, std::ios::binary);
  for (int tid = 0; tid < eg.num_thr; tid++) {
    std::ifstream fin_pos(eg.outfile_pos + '.' + std::to_string(tid),
                          std::ios::binary);
    fin_pos.read((char *)&abs_pos_thr, sizeof(uint64_t));
    while (!fin_pos.eof()) {
      abs_pos_thr += abs_pos;
      fout_pos.write((char *)&abs_pos_thr, sizeof(uint64_t));
      fin_pos.read((char *)&abs_pos_thr, sizeof(uint64_t));
    }
    fin_pos.close();
    remove((eg.outfile_pos + '.' + std::to_string(tid)).c_str());
    abs_pos += file_len_seq_thr[tid];
  }
  fout_pos.close();
  delete[] file_len_seq_thr;

  std::cout << "Encoding done:\n";
  std::cout << matched_s << " singleton reads were aligned\n";
  std::cout << matched_N << " reads with N were aligned\n";
  return;
}

template <size_t bitset_size>
void setglobalarrays(encoder_global &eg, encoder_global_b<bitset_size> &egb) {
  for (int i = 0; i < 63; i++) egb.mask63[i] = 1;
  for (int i = 0; i < eg.max_readlen; i++) {
    egb.basemask[i][(uint8_t)'A'][3 * i] = 0;
    egb.basemask[i][(uint8_t)'A'][3 * i + 1] = 0;
    egb.basemask[i][(uint8_t)'A'][3 * i + 2] = 0;
    egb.basemask[i][(uint8_t)'C'][3 * i] = 0;
    egb.basemask[i][(uint8_t)'C'][3 * i + 1] = 0;
    egb.basemask[i][(uint8_t)'C'][3 * i + 2] = 1;
    egb.basemask[i][(uint8_t)'G'][3 * i] = 0;
    egb.basemask[i][(uint8_t)'G'][3 * i + 1] = 1;
    egb.basemask[i][(uint8_t)'G'][3 * i + 2] = 0;
    egb.basemask[i][(uint8_t)'T'][3 * i] = 0;
    egb.basemask[i][(uint8_t)'T'][3 * i + 1] = 1;
    egb.basemask[i][(uint8_t)'T'][3 * i + 2] = 1;
    egb.basemask[i][(uint8_t)'N'][3 * i] = 1;
    egb.basemask[i][(uint8_t)'N'][3 * i + 1] = 0;
    egb.basemask[i][(uint8_t)'N'][3 * i + 2] = 0;
  }

  // enc_noise uses substitution statistics from Minoche et al.
  eg.enc_noise[(uint8_t)'A'][(uint8_t)'C'] = '0';
  eg.enc_noise[(uint8_t)'A'][(uint8_t)'G'] = '1';
  eg.enc_noise[(uint8_t)'A'][(uint8_t)'T'] = '2';
  eg.enc_noise[(uint8_t)'A'][(uint8_t)'N'] = '3';
  eg.enc_noise[(uint8_t)'C'][(uint8_t)'A'] = '0';
  eg.enc_noise[(uint8_t)'C'][(uint8_t)'G'] = '1';
  eg.enc_noise[(uint8_t)'C'][(uint8_t)'T'] = '2';
  eg.enc_noise[(uint8_t)'C'][(uint8_t)'N'] = '3';
  eg.enc_noise[(uint8_t)'G'][(uint8_t)'T'] = '0';
  eg.enc_noise[(uint8_t)'G'][(uint8_t)'A'] = '1';
  eg.enc_noise[(uint8_t)'G'][(uint8_t)'C'] = '2';
  eg.enc_noise[(uint8_t)'G'][(uint8_t)'N'] = '3';
  eg.enc_noise[(uint8_t)'T'][(uint8_t)'G'] = '0';
  eg.enc_noise[(uint8_t)'T'][(uint8_t)'C'] = '1';
  eg.enc_noise[(uint8_t)'T'][(uint8_t)'A'] = '2';
  eg.enc_noise[(uint8_t)'T'][(uint8_t)'N'] = '3';
  eg.enc_noise[(uint8_t)'N'][(uint8_t)'A'] = '0';
  eg.enc_noise[(uint8_t)'N'][(uint8_t)'G'] = '1';
  eg.enc_noise[(uint8_t)'N'][(uint8_t)'C'] = '2';
  eg.enc_noise[(uint8_t)'N'][(uint8_t)'T'] = '3';
  return;
}

template <size_t bitset_size>
void readsingletons(std::bitset<bitset_size> *read, uint32_t *order_s,
                    uint16_t *read_lengths_s, const encoder_global &eg,
                    const encoder_global_b<bitset_size> &egb) {
  // not parallelized right now since these are very small number of reads
  std::ifstream f(eg.infile + ".singleton", std::ifstream::in);
  std::string s;
  for (uint32_t i = 0; i < eg.numreads_s; i++) {
    std::getline(f, s);
    read_lengths_s[i] = s.length();
    stringtobitset<bitset_size>(s, read_lengths_s[i], read[i], egb.basemask);
  }
  f.close();
  remove((eg.infile + ".singleton").c_str());
  f.open(eg.infile_N);
  for (uint32_t i = eg.numreads_s; i < eg.numreads_s + eg.numreads_N; i++) {
    std::getline(f, s);
    read_lengths_s[i] = s.length();
    stringtobitset<bitset_size>(s, read_lengths_s[i], read[i], egb.basemask);
  }
  std::ifstream f_order_s(eg.infile_order + ".singleton", std::ios::binary);
  for (uint32_t i = 0; i < eg.numreads_s; i++)
    f_order_s.read((char *)&order_s[i], sizeof(uint32_t));
  f_order_s.close();
  remove((eg.infile_order + ".singleton").c_str());
  std::ifstream f_order_N(eg.infile_order_N, std::ios::binary);
  for (uint32_t i = eg.numreads_s; i < eg.numreads_s + eg.numreads_N; i++)
    f_order_N.read((char *)&order_s[i], sizeof(uint32_t));
  f_order_N.close();
}

template <size_t bitset_size>
void encoder_main(const std::string &temp_dir, const compression_params &cp) {
  encoder_global_b<bitset_size> *egb_ptr =
      new encoder_global_b<bitset_size>(cp.max_readlen);
  encoder_global *eg_ptr = new encoder_global;
  encoder_global_b<bitset_size> &egb = *egb_ptr;
  encoder_global &eg = *eg_ptr;

  eg.basedir = temp_dir;
  eg.infile = eg.basedir + "/temp.dna";
  eg.infile_pos = eg.basedir + "/temppos.txt";
  eg.infile_flag = eg.basedir + "/tempflag.txt";
  eg.infile_order = eg.basedir + "/read_order.bin";
  eg.infile_order_N = eg.basedir + "/read_order_N.bin";
  eg.infile_RC = eg.basedir + "/read_rev.txt";
  eg.infile_readlength = eg.basedir + "/read_lengths.bin";
  eg.infile_N = eg.basedir + "/input_N.dna";
  eg.outfile_seq = eg.basedir + "/read_seq.bin";
  eg.outfile_pos = eg.basedir + "/read_pos.bin";
  eg.outfile_noise = eg.basedir + "/read_noise.txt";
  eg.outfile_noisepos = eg.basedir + "/read_noisepos.bin";
  eg.outfile_unaligned = eg.basedir + "/read_unaligned.txt";

  eg.max_readlen = cp.max_readlen;
  eg.num_thr = cp.num_thr;

  omp_set_num_threads(eg.num_thr);
  getDataParams(eg, cp);  // populate numreads
  setglobalarrays<bitset_size>(eg, egb);
  std::bitset<bitset_size> *read =
      new std::bitset<bitset_size>[eg.numreads_s + eg.numreads_N];
  uint32_t *order_s = new uint32_t[eg.numreads_s + eg.numreads_N];
  uint16_t *read_lengths_s = new uint16_t[eg.numreads_s + eg.numreads_N];
  readsingletons<bitset_size>(read, order_s, read_lengths_s, eg, egb);
  correct_order(order_s, eg);

  bbhashdict *dict = new bbhashdict[eg.numdict_s];
  if (eg.max_readlen > 50) {
    dict[0].start = 0;
    dict[0].end = 20;
    dict[1].start = 21;
    dict[1].end = 41;
  } else {
    dict[0].start = 0;
    dict[0].end = 20 * eg.max_readlen / 50;
    dict[1].start = 20 * eg.max_readlen / 50 + 1;
    dict[1].end = 41 * eg.max_readlen / 50;
  }

  constructdictionary<bitset_size>(read, dict, read_lengths_s, eg.numdict_s,
                                   eg.numreads_s + eg.numreads_N, 3, eg.basedir,
                                   eg.num_thr);
  encode<bitset_size>(read, dict, order_s, read_lengths_s, eg, egb);
  remove(eg.infile_N.c_str());
  delete[] read;
  delete[] dict;
  delete[] order_s;
  delete[] read_lengths_s;
  delete eg_ptr;
  delete egb_ptr;
}

}  // namespace spring

#endif  // SPRING_ENCODER_H_
