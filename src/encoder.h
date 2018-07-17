#ifndef SPRING_ENCODER_H_
#define SPRING_ENCODER_H_

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
#include "algorithms/SPRING/util.h"

namespace spring {

template <size_t bitset_size>
struct encoder_global_b {
  std::bitset<bitset_size>
      basemask[MAX_READ_LEN][128];  // bitset for A,G,C,T,N at each position
  // used in stringtobitset, and bitsettostring
  std::bitset<bitset_size>
      mask63;  // bitset with 63 bits set to 1 (used in bitsettostring for
               // conversion to ullong)
};

struct encoder_global {
  uint32_t numreads, numreads_s, numreads_N;
  int numdict_s = 2;
  int thresh_s = 24;
  int maxsearch = 1000;
  int max_readlen, num_thr;

  std::string basedir;
  std::string infile;
  std::string infile_flag;
  std::string infile_pos;
  std::string infile_seq;
  std::string infile_RC;
  std::string infile_readlength;
  std::string infile_N;
  std::string outfile_N;
  std::string outfile_seq;
  std::string outfile_pos;
  std::string outfile_noise;
  std::string outfile_noisepos;
  std::string outfile_singleton;
  std::string infile_order;
  std::string infile_order_N;
  std::string infilenumreads;

  char longtochar[5] = {'A', 'C', 'G', 'T', 'N'};
  long chartolong[128];
  char enc_noise[128][128];

  // Some global arrays (some initialized in setglobalarrays())
  char revinttochar[8] = {'A', 'N', 'G', '#',
                          'C', '#', 'T', '#'};  // used in bitsettostring
  char chartorevchar[128];  // A-T etc for reverse complement
};

struct contig_reads {
  std::string read;
  int64_t pos;
  char RC;
  uint32_t order;
  uint8_t read_length;
};

std::string buildcontig(std::list<contig_reads> &current_contig,
                        uint32_t list_size, encoder_global &eg);

void writecontig(std::string &ref, std::list<contig_reads> &current_contig,
                 std::ofstream &f_seq, std::ofstream &f_pos,
                 std::ofstream &f_noise, std::ofstream &f_noisepos,
                 std::ofstream &f_order, std::ofstream &f_RC,
                 std::ofstream &f_readlength, encoder_global &eg);

void packbits(encoder_global &eg);

void getDataParams(encoder_global &eg);

void correct_order(uint32_t *order_s, encoder_global &eg);

template <size_t bitset_size>
std::string bitsettostring(std::bitset<bitset_size> b, uint8_t readlen,
                           encoder_global &eg,
                           encoder_global_b<bitset_size> &egb) {
  char s[MAX_READ_LEN + 1];
  s[readlen] = '\0';
  unsigned long long ull;
  for (int i = 0; i < 3 * readlen / 63 + 1; i++) {
    ull = (b & egb.mask63).to_ullong();
    b >>= 63;
    for (int j = 21 * i; j < 21 * i + 21 && j < readlen; j++) {
      s[j] = eg.revinttochar[ull % 8];
      ull /= 8;
    }
  }
  std::string s1 = s;
  return s1;
}

template <size_t bitset_size>
void encode(std::bitset<bitset_size> *read, bbhashdict *dict, uint32_t *order_s,
            uint8_t *read_lengths_s, encoder_global &eg,
            encoder_global_b<bitset_size> &egb) {
  omp_lock_t *read_lock = new omp_lock_t[eg.numreads_s + eg.numreads_N];
  omp_lock_t *dict_lock = new omp_lock_t[eg.numreads_s + eg.numreads_N];
  for (uint64_t j = 0; j < eg.numreads_s + eg.numreads_N; j++) {
    omp_init_lock(&read_lock[j]);
    omp_init_lock(&dict_lock[j]);
  }
  bool *remainingreads = new bool[eg.numreads_s + eg.numreads_N];
  std::fill(remainingreads, remainingreads + eg.numreads_s + eg.numreads_N, 1);

  std::bitset<bitset_size> mask1[eg.numdict_s];
  generateindexmasks<bitset_size>(mask1, dict, eg.numdict_s, 3);
  std::bitset<bitset_size> mask[MAX_READ_LEN][MAX_READ_LEN];
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
    bool flag = 0;
    // flag to check if match was found or not
    std::string current, ref;
    std::bitset<bitset_size> forward_bitset, reverse_bitset, b;
    char c, rc;
    std::list<contig_reads> current_contig;
    int64_t p;
    uint8_t rl;
    uint32_t ord, list_size = 0;  // list_size variable introduced because
                                  // list::size() was running very slowly
                                  // on UIUC machine
    std::list<uint32_t> deleted_rids[eg.numdict_s];
    bool done = false;
    while (!done) {
      if (!(in_flag >> c)) done = true;
      if (!done) {
        std::getline(f, current);
        rc = in_RC.get();
        in_pos.read((char *)&p, sizeof(int64_t));
        in_order.read((char *)&ord, sizeof(uint32_t));
        in_readlength.read((char *)&rl, sizeof(uint8_t));
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

          ref = buildcontig(current_contig, list_size, eg);
          // try to align the singleton reads to ref
          // first create bitsets from first readlen positions of ref
          forward_bitset.reset();
          reverse_bitset.reset();
          if ((int64_t)ref.size() >= eg.max_readlen) {
            stringtobitset(ref.substr(0, eg.max_readlen), eg.max_readlen,
                           forward_bitset, egb.basemask);
            stringtobitset(reverse_complement(ref.substr(0, eg.max_readlen),
                                              eg.max_readlen, eg.chartorevchar),
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
                         i >= dictidx[0] && i >= dictidx[1] - eg.maxsearch;
                         i--) {
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
                      if (hamming <= eg.thresh_s) {
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
                                          read[rid], read_lengths_s[rid], eg,
                                          egb),
                                      read_lengths_s[rid], eg.chartorevchar)
                                : bitsettostring<bitset_size>(
                                      read[rid], read_lengths_s[rid], eg, egb);
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
                reverse_bitset |= egb.basemask[0][(uint8_t)eg.chartorevchar[(
                    uint8_t)ref[j + eg.max_readlen]]];
              }

            }  // end for
          }    // end if
          // sort contig according to pos
          current_contig.sort([](const contig_reads &a, const contig_reads &b) {
            return a.pos < b.pos;
          });
          writecontig(ref, current_contig, f_seq, f_pos, f_noise, f_noisepos,
                      f_order, f_RC, f_readlength, eg);
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
  }

  // Combine files produced by the threads
  std::ofstream f_order(eg.infile_order);
  std::ofstream f_readlength(eg.infile_readlength);
  for (int tid = 0; tid < eg.num_thr; tid++) {
    std::ifstream in_order(eg.infile_order + '.' + std::to_string(tid) +
                           ".tmp");
    std::ifstream in_readlength(eg.infile_readlength + '.' +
                                std::to_string(tid) + ".tmp");
    f_order << in_order.rdbuf();
    f_order.clear();  // clear error flag in case in_order is empty
    f_readlength << in_readlength.rdbuf();
    f_readlength.clear();  // clear error flag in case in_readlength is empty

    remove((eg.infile_order + '.' + std::to_string(tid)).c_str());
    remove((eg.infile_order + '.' + std::to_string(tid) + ".tmp").c_str());
    remove((eg.infile_readlength + '.' + std::to_string(tid)).c_str());
    remove((eg.infile_readlength + '.' + std::to_string(tid) + ".tmp").c_str());
    rename((eg.infile_RC + '.' + std::to_string(tid) + ".tmp").c_str(),
           (eg.infile_RC + '.' + std::to_string(tid)).c_str());
    remove((eg.infile_RC + '.' + std::to_string(tid) + ".tmp").c_str());
    remove((eg.infile_flag + '.' + std::to_string(tid)).c_str());
    remove((eg.infile_pos + '.' + std::to_string(tid)).c_str());
    remove((eg.infile + '.' + std::to_string(tid)).c_str());
  }
  f_order.close();
  f_readlength.close();
  // write remaining singleton reads now
  std::ofstream f_singleton(eg.outfile_singleton);
  f_order.open(eg.infile_order, std::ios::binary | std::ofstream::app);
  f_readlength.open(eg.infile_readlength,
                    std::ios::binary | std::ofstream::app);
  std::ofstream f_N(eg.outfile_N);
  uint32_t matched_s = eg.numreads_s;
  for (uint32_t i = 0; i < eg.numreads_s; i++)
    if (remainingreads[i] == 1) {
      matched_s--;
      f_order.write((char *)&order_s[i], sizeof(uint32_t));
      f_readlength.write((char *)&read_lengths_s[i], sizeof(uint8_t));
      f_singleton << bitsettostring<bitset_size>(read[i], read_lengths_s[i], eg,
                                                 egb);
    }
  uint32_t matched_N = eg.numreads_N;
  for (uint32_t i = eg.numreads_s; i < eg.numreads_s + eg.numreads_N; i++)
    if (remainingreads[i] == 1) {
      matched_N--;
      f_N << bitsettostring<bitset_size>(read[i], read_lengths_s[i], eg, egb);
      f_order.write((char *)&order_s[i], sizeof(uint32_t));
      f_readlength.write((char *)&read_lengths_s[i], sizeof(uint8_t));
    }
  f_order.close();
  f_readlength.close();
  f_N.close();
  f_singleton.close();
  delete[] remainingreads;
  delete[] dict_lock;
  delete[] read_lock;
  packbits(eg);
  std::cout << "Encoding done:\n";
  std::cout << matched_s << " singleton reads were aligned\n";
  std::cout << matched_N << " reads with N were aligned\n";
  return;
}

template <size_t bitset_size>
void setglobalarrays(encoder_global &eg, encoder_global_b<bitset_size> &egb) {
  eg.chartorevchar[(uint8_t)'A'] = 'T';
  eg.chartorevchar[(uint8_t)'C'] = 'G';
  eg.chartorevchar[(uint8_t)'G'] = 'C';
  eg.chartorevchar[(uint8_t)'T'] = 'A';
  eg.chartorevchar[(uint8_t)'N'] = 'N';

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
  eg.chartolong[(uint8_t)'A'] = 0;
  eg.chartolong[(uint8_t)'C'] = 1;
  eg.chartolong[(uint8_t)'G'] = 2;
  eg.chartolong[(uint8_t)'T'] = 3;
  eg.chartolong[(uint8_t)'N'] = 4;
  return;
}

template <size_t bitset_size>
void readsingletons(std::bitset<bitset_size> *read, uint32_t *order_s,
                    uint8_t *read_lengths_s, encoder_global &eg,
                    encoder_global_b<bitset_size> &egb) {
  // not parallelized right now since these are very small number of reads
  std::ifstream f(eg.infile + ".singleton", std::ifstream::in);
  std::string s;
  for (uint32_t i = 0; i < eg.numreads_s; i++) {
    std::getline(f, s);
    read_lengths_s[i] = s.length();
    stringtobitset<bitset_size>(s, read_lengths_s[i], read[i], egb.basemask);
  }
  f.close();
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
  std::ifstream f_order_N(eg.infile_order_N, std::ios::binary);
  for (uint32_t i = eg.numreads_s; i < eg.numreads_s + eg.numreads_N; i++)
    f_order_N.read((char *)&order_s[i], sizeof(uint32_t));
  f_order_N.close();
}

template <size_t bitset_size>
void encoder_main(const std::string &working_dir, int max_readlen, int num_thr) {
  encoder_global_b<bitset_size> *egb_ptr = new encoder_global_b<bitset_size>;
  encoder_global *eg_ptr = new encoder_global;
  encoder_global_b<bitset_size> &egb = *egb_ptr;
  encoder_global &eg = *eg_ptr;

  eg.basedir = working_dir;
  eg.infile = eg.basedir + "/temp.dna";
  eg.infile_pos = eg.basedir + "/temppos.txt";
  eg.infile_flag = eg.basedir + "/tempflag.txt";
  eg.infile_order = eg.basedir + "/read_order.bin";
  eg.infile_order_N = eg.basedir + "/read_order_N.bin";
  eg.infile_RC = eg.basedir + "/read_rev.txt";
  eg.infile_readlength = eg.basedir + "/read_lengths.bin";
  eg.infile_N = eg.basedir + "/input_N.dna";
  eg.outfile_N = eg.basedir + "/unaligned_N.txt";
  eg.outfile_seq = eg.basedir + "/read_seq.txt";
  eg.outfile_pos = eg.basedir + "/read_pos.txt";
  eg.outfile_noise = eg.basedir + "/read_noise.txt";
  eg.outfile_noisepos = eg.basedir + "/read_noisepos.txt";
  eg.outfile_singleton = eg.basedir + "/unaligned_singleton.txt";
  eg.infilenumreads = eg.basedir + "/numreads.bin";

  eg.max_readlen = max_readlen;
  eg.num_thr = num_thr;

  omp_set_num_threads(eg.num_thr);
  getDataParams(eg);  // populate numreads
  setglobalarrays<bitset_size>(eg, egb);
  std::bitset<bitset_size> *read =
      new std::bitset<bitset_size>[eg.numreads_s + eg.numreads_N];
  uint32_t *order_s = new uint32_t[eg.numreads_s + eg.numreads_N];
  uint8_t *read_lengths_s = new uint8_t[eg.numreads_s + eg.numreads_N];
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