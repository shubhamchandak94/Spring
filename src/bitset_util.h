/*
* Copyright 2018 University of Illinois Board of Trustees and Stanford University. All Rights Reserved.
* Licensed under the “Non-exclusive Research Use License for SPRING Software” license (the "License");
* You may not use this file except in compliance with the License.
* The License is included in the distribution as license.pdf file.
 
* Software distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and limitations under the License.
*/

#ifndef SPRING_BITSET_UTIL_H_
#define SPRING_BITSET_UTIL_H_

#include <omp.h>
#include <algorithm>
#include <bitset>
#include <fstream>
#include <string>
#include "BooPHF.h"
#include "params.h"
namespace spring {

typedef boomphf::SingleHashFunctor<u_int64_t> hasher_t;
typedef boomphf::mphf<u_int64_t, hasher_t> boophf_t;

class bbhashdict {
 public:
  boophf_t *bphf;
  int start;
  int end;
  uint32_t numkeys;
  uint32_t dict_numreads;  // number of reads in this dict (for variable length)
  uint32_t *startpos;
  uint32_t *read_id;
  bool *empty_bin;
  void findpos(int64_t *dictidx, const uint64_t &startposidx);
  void remove(int64_t *dictidx, const uint64_t &startposidx,
              const int64_t current);
  bbhashdict() {
    bphf = NULL;
    startpos = NULL;
    read_id = NULL;
  }
  ~bbhashdict() {
    delete[] startpos;
    delete[] read_id;
    delete[] empty_bin;
    delete bphf;
  }
};

template <size_t bitset_size>
void stringtobitset(const std::string &s, const uint16_t readlen,
                    std::bitset<bitset_size> &b,
                    std::bitset<bitset_size> **basemask) {
  for (int i = 0; i < readlen; i++) b |= basemask[i][(uint8_t)s[i]];
}

template <size_t bitset_size>
void generateindexmasks(std::bitset<bitset_size> *mask1, bbhashdict *dict,
                        int numdict, int bpb) {
  for (int j = 0; j < numdict; j++) mask1[j].reset();
  for (int j = 0; j < numdict; j++)
    for (int i = bpb * dict[j].start; i < bpb * (dict[j].end + 1); i++)
      mask1[j][i] = 1;
  return;
}

template <size_t bitset_size>
void constructdictionary(std::bitset<bitset_size> *read, bbhashdict *dict,
                         uint16_t *read_lengths, const int numdict,
                         const uint32_t &numreads, const int bpb,
                         const std::string &basedir, const int &num_thr) {
  std::bitset<bitset_size> *mask = new std::bitset<bitset_size>[numdict];
  generateindexmasks<bitset_size>(mask, dict, numdict, bpb);
  for (int j = 0; j < numdict; j++) {
    uint64_t *ull = new uint64_t[numreads];
#pragma omp parallel
    {
      std::bitset<bitset_size> b;
      int tid = omp_get_thread_num();
      uint64_t i, stop;
      i = uint64_t(tid) * numreads / omp_get_num_threads();
      stop = uint64_t(tid + 1) * numreads / omp_get_num_threads();
      if (tid == omp_get_num_threads() - 1) stop = numreads;
      // compute keys and and store in ull
      for (; i < stop; i++) {
        b = read[i] & mask[j];
        ull[i] = (b >> bpb * dict[j].start).to_ullong();
      }
    }  // parallel end

    // remove keys corresponding to reads shorter than dict_end[j]
    dict[j].dict_numreads = 0;
    for (uint32_t i = 0; i < numreads; i++) {
      if (read_lengths[i] > dict[j].end) {
        ull[dict[j].dict_numreads] = ull[i];
        dict[j].dict_numreads++;
      }
    }

// write ull to file
#pragma omp parallel
    {
      int tid = omp_get_thread_num();
      std::ofstream foutkey(
          basedir + std::string("/keys.bin.") + std::to_string(tid),
          std::ios::binary);
      uint64_t i, stop;
      i = uint64_t(tid) * dict[j].dict_numreads / omp_get_num_threads();
      stop = uint64_t(tid + 1) * dict[j].dict_numreads / omp_get_num_threads();
      if (tid == omp_get_num_threads() - 1) stop = dict[j].dict_numreads;
      for (; i < stop; i++) foutkey.write((char *)&ull[i], sizeof(uint64_t));
      foutkey.close();
    }  // parallel end

    // deduplicating ull
    std::sort(ull, ull + dict[j].dict_numreads);
    uint32_t k = 0;
    for (uint32_t i = 1; i < dict[j].dict_numreads; i++)
      if (ull[i] != ull[k]) ull[++k] = ull[i];
    dict[j].numkeys = k + 1;
    // construct mphf
    auto data_iterator =
        boomphf::range(static_cast<const u_int64_t *>(ull),
                       static_cast<const u_int64_t *>(ull + dict[j].numkeys));
    double gammaFactor = 5.0;  // balance between speed and memory
    dict[j].bphf = new boomphf::mphf<u_int64_t, hasher_t>(
        dict[j].numkeys, data_iterator, num_thr, gammaFactor, true, false);

    delete[] ull;

// compute hashes for all reads
#pragma omp parallel
    {
      int tid = omp_get_thread_num();
      std::ifstream finkey(
          basedir + std::string("/keys.bin.") + std::to_string(tid),
          std::ios::binary);
      std::ofstream fouthash(basedir + std::string("/hash.bin.") +
                                 std::to_string(tid) + '.' + std::to_string(j),
                             std::ios::binary);
      uint64_t currentkey, currenthash;
      uint64_t i, stop;
      i = uint64_t(tid) * dict[j].dict_numreads / omp_get_num_threads();
      stop = uint64_t(tid + 1) * dict[j].dict_numreads / omp_get_num_threads();
      if (tid == omp_get_num_threads() - 1) stop = dict[j].dict_numreads;
      for (; i < stop; i++) {
        finkey.read((char *)&currentkey, sizeof(uint64_t));
        currenthash = (dict[j].bphf)->lookup(currentkey);
        fouthash.write((char *)&currenthash, sizeof(uint64_t));
      }
      finkey.close();
      remove(
          (basedir + std::string("/keys.bin.") + std::to_string(tid)).c_str());
      fouthash.close();
    }  // parallel end
  }

  // for rest of the function, use numdict threads to parallelize
  omp_set_num_threads(std::min(numdict, num_thr));
#pragma omp parallel
  {
#pragma omp for
    for (int j = 0; j < numdict; j++) {
      // fill startpos by first storing numbers and then doing cumulative sum
      dict[j].startpos =
          new uint32_t[dict[j].numkeys +
                       1]();  // 1 extra to store end pos of last key
      uint64_t currenthash;
      for (int tid = 0; tid < num_thr; tid++) {
        std::ifstream finhash(basedir + std::string("/hash.bin.") +
                                  std::to_string(tid) + '.' + std::to_string(j),
                              std::ios::binary);
        finhash.read((char *)&currenthash, sizeof(uint64_t));
        while (!finhash.eof()) {
          dict[j].startpos[currenthash + 1]++;
          finhash.read((char *)&currenthash, sizeof(uint64_t));
        }
        finhash.close();
      }

      dict[j].empty_bin = new bool[dict[j].numkeys]();
      for (uint32_t i = 1; i < dict[j].numkeys; i++)
        dict[j].startpos[i] = dict[j].startpos[i] + dict[j].startpos[i - 1];

      // insert elements in the dict array
      dict[j].read_id = new uint32_t[dict[j].dict_numreads];
      uint32_t i = 0;
      for (int tid = 0; tid < num_thr; tid++) {
        std::ifstream finhash(basedir + std::string("/hash.bin.") +
                                  std::to_string(tid) + '.' + std::to_string(j),
                              std::ios::binary);
        finhash.read((char *)&currenthash, sizeof(uint64_t));
        while (!finhash.eof()) {
          while (read_lengths[i] <= dict[j].end) i++;
          dict[j].read_id[dict[j].startpos[currenthash]++] = i;
          i++;
          finhash.read((char *)&currenthash, sizeof(uint64_t));
        }
        finhash.close();
        remove((basedir + std::string("/hash.bin.") + std::to_string(tid) +
                '.' + std::to_string(j))
                   .c_str());
      }

      // correcting startpos array modified during insertion
      for (int64_t keynum = dict[j].numkeys; keynum >= 1; keynum--)
        dict[j].startpos[keynum] = dict[j].startpos[keynum - 1];
      dict[j].startpos[0] = 0;
    }  // for end
  }    // parallel end
  omp_set_num_threads(num_thr);
  delete[] mask;
  return;
}

template <size_t bitset_size>
void generatemasks(std::bitset<bitset_size> **mask, const int max_readlen,
                   const int bpb) {
  // mask for zeroing the end bits (needed while reordering to compute Hamming
  // distance between shifted reads)
  for (int i = 0; i < max_readlen; i++) {
    for (int j = 0; j < max_readlen; j++) {
      mask[i][j].reset();
      for (int k = bpb * i; k < bpb * max_readlen - bpb * j; k++)
        mask[i][j][k] = 1;
    }
  }
  return;
}

template <size_t bitset_size>
void chartobitset(char *s, const int readlen, std::bitset<bitset_size> &b,
                  std::bitset<bitset_size> **basemask) {
  b.reset();
  for (int i = 0; i < readlen; i++) b |= basemask[i][(uint8_t)s[i]];
  return;
}

}  // namespace spring

#endif  // SPRING_BITSET_UTIL_H_
