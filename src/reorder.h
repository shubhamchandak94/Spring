#ifndef SPRING_REORDER_H_
#define SPRING_REORDER_H_

#include <omp.h>
#include <algorithm>
#include <atomic>
#include <bitset>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "util.h"

namespace spring {

template <size_t bitset_size>
struct reorder_global {
  uint32_t num_locks =
      0x1000000;  // limits on number of locks (power of 2 for fast mod)

  uint32_t numreads = 0;

  int maxshift, numdict = 2, maxsearch = 1000, num_thr, max_readlen;
  uint thresh = 4;
  std::string basedir;
  std::string infile;
  std::string infilenumreads;
  std::string outfile;
  std::string outfileRC;
  std::string outfileflag;
  std::string outfilepos;
  std::string outfileorder;
  std::string outfilereadlength;

  // Some global arrays (some initialized in setglobalarrays())
  char revinttochar[4] = {'A', 'G', 'C', 'T'};  // used in bitsettostring
  char inttochar[4] = {'A', 'C', 'G', 'T'};
  char chartorevchar[128];  // A-T etc for reverse complement
  int chartoint[128];       // A-0,C-1 etc. used in updaterefcount

  std::bitset<bitset_size>
      basemask[MAX_READ_LEN][128];  // bitset for A,G,C,T at each position
  // used in stringtobitset, chartobitset and bitsettostring
  std::bitset<bitset_size>
      mask64;  // bitset with 64 bits set to 1 (used in bitsettostring for
               // conversion to ullong)
};

template <size_t bitset_size>
void bitsettostring(std::bitset<bitset_size> b, char *s, uint8_t readlen,
                    reorder_global<bitset_size> &rg) {
  unsigned long long ull;
  for (int i = 0; i < 2 * readlen / 64 + 1; i++) {
    ull = (b & rg.mask64).to_ullong();
    b >>= 64;
    for (int j = 32 * i; j < 32 * i + 32 && j < readlen; j++) {
      s[j] = rg.revinttochar[ull % 4];
      ull /= 4;
    }
  }
  s[readlen] = '\0';
  return;
}

template <size_t bitset_size>
void setglobalarrays(reorder_global<bitset_size> &rg) {
  rg.chartorevchar[(uint8_t)'A'] = 'T';
  rg.chartorevchar[(uint8_t)'C'] = 'G';
  rg.chartorevchar[(uint8_t)'G'] = 'C';
  rg.chartorevchar[(uint8_t)'T'] = 'A';
  rg.chartoint[(uint8_t)'A'] = 0;
  rg.chartoint[(uint8_t)'C'] = 1;
  rg.chartoint[(uint8_t)'G'] = 2;
  rg.chartoint[(uint8_t)'T'] = 3;

  for (int i = 0; i < 64; i++) rg.mask64[i] = 1;
  for (int i = 0; i < rg.max_readlen; i++) {
    rg.basemask[i][(uint8_t)'A'][2 * i] = 0;
    rg.basemask[i][(uint8_t)'A'][2 * i + 1] = 0;
    rg.basemask[i][(uint8_t)'C'][2 * i] = 0;
    rg.basemask[i][(uint8_t)'C'][2 * i + 1] = 1;
    rg.basemask[i][(uint8_t)'G'][2 * i] = 1;
    rg.basemask[i][(uint8_t)'G'][2 * i + 1] = 0;
    rg.basemask[i][(uint8_t)'T'][2 * i] = 1;
    rg.basemask[i][(uint8_t)'T'][2 * i + 1] = 1;
  }
  return;
}

template <size_t bitset_size>
void updaterefcount(std::bitset<bitset_size> &cur,
                    std::bitset<bitset_size> &ref,
                    std::bitset<bitset_size> &revref, int count[][MAX_READ_LEN],
                    bool resetcount, bool rev, int shift, uint8_t cur_readlen,
                    int &ref_len, reorder_global<bitset_size> &rg)
// for var length, shift represents shift of start positions, if read length is
// small, may not need to shift actually
{
  char s[MAX_READ_LEN + 1], s1[MAX_READ_LEN + 1], *current;
  bitsettostring<bitset_size>(cur, s, cur_readlen, rg);
  if (rev == false)
    current = s;
  else {
    reverse_complement(s, s1, cur_readlen, rg.chartorevchar);
    current = s1;
  }

  if (resetcount == true)  // resetcount - unmatched read so start over
  {
    std::memset(count, 0, sizeof(count[0][0]) * 4 * MAX_READ_LEN);
    for (int i = 0; i < cur_readlen; i++) {
      count[rg.chartoint[(uint8_t)current[i]]][i] = 1;
    }
    ref_len = cur_readlen;
  } else {
    if (!rev) {
      // shift count
      for (int i = 0; i < ref_len - shift; i++) {
        for (int j = 0; j < 4; j++) count[j][i] = count[j][i + shift];
        if (i < cur_readlen) count[rg.chartoint[(uint8_t)current[i]]][i] += 1;
      }

      // for the new positions set count to 1
      for (int i = ref_len - shift; i < cur_readlen; i++) {
        for (int j = 0; j < 4; j++) count[j][i] = 0;
        count[rg.chartoint[(uint8_t)current[i]]][i] = 1;
      }
      ref_len = std::max<int>(ref_len - shift, cur_readlen);
    } else  // reverse case is quite complicated in the variable length case
    {
      if (cur_readlen - shift >= ref_len) {
        for (int i = cur_readlen - shift - ref_len; i < cur_readlen - shift;
             i++) {
          for (int j = 0; j < 4; j++)
            count[j][i] = count[j][i - (cur_readlen - shift - ref_len)];
          count[rg.chartoint[(uint8_t)current[i]]][i] += 1;
        }
        for (int i = 0; i < cur_readlen - shift - ref_len; i++) {
          for (int j = 0; j < 4; j++) count[j][i] = 0;
          count[rg.chartoint[(uint8_t)current[i]]][i] = 1;
        }
        for (int i = cur_readlen - shift; i < cur_readlen; i++) {
          for (int j = 0; j < 4; j++) count[j][i] = 0;
          count[rg.chartoint[(uint8_t)current[i]]][i] = 1;
        }
        ref_len = cur_readlen;
      } else if (ref_len + shift <= rg.max_readlen) {
        for (int i = ref_len - cur_readlen + shift; i < ref_len; i++)
          count[rg.chartoint[(
              uint8_t)current[i - (ref_len - cur_readlen + shift)]]][i] += 1;
        for (int i = ref_len; i < ref_len + shift; i++) {
          for (int j = 0; j < 4; j++) count[j][i] = 0;
          count[rg.chartoint[(
              uint8_t)current[i - (ref_len - cur_readlen + shift)]]][i] = 1;
        }
        ref_len = ref_len + shift;
      } else {
        for (int i = 0; i < rg.max_readlen - shift; i++) {
          for (int j = 0; j < 4; j++)
            count[j][i] = count[j][i + (ref_len + shift - rg.max_readlen)];
        }
        for (int i = rg.max_readlen - cur_readlen; i < rg.max_readlen - shift;
             i++)
          count[rg.chartoint[(
              uint8_t)current[i - (rg.max_readlen - cur_readlen)]]][i] += 1;
        for (int i = rg.max_readlen - shift; i < rg.max_readlen; i++) {
          for (int j = 0; j < 4; j++) count[j][i] = 0;
          count[rg.chartoint[(
              uint8_t)current[i - (rg.max_readlen - cur_readlen)]]][i] = 1;
        }
        ref_len = rg.max_readlen;
      }
    }
    // find max of each position to get ref

    for (int i = 0; i < ref_len; i++) {
      int max = 0, indmax = 0;
      for (int j = 0; j < 4; j++)
        if (count[j][i] > max) {
          max = count[j][i];
          indmax = j;
        }
      current[i] = rg.inttochar[indmax];
    }
  }
  chartobitset<bitset_size>(current, ref_len, ref, rg.basemask);
  char revcurrent[MAX_READ_LEN];
  reverse_complement(current, revcurrent, ref_len, rg.chartorevchar);
  chartobitset<bitset_size>(revcurrent, ref_len, revref, rg.basemask);

  return;
}

template <size_t bitset_size>
void readDnaFile(std::bitset<bitset_size> *read, uint8_t *read_lengths,
                 reorder_global<bitset_size> &rg) {
#pragma omp parallel
  {
    uint32_t tid = omp_get_thread_num();
    std::ifstream f(rg.infile, std::ifstream::in);
    std::string s;
    uint32_t i = 0;
    while (std::getline(f, s)) {
      if (i % rg.num_thr == tid) {
        read_lengths[i] = (uint8_t)s.length();
        stringtobitset<bitset_size>(s, read_lengths[i], read[i], rg.basemask);
        i++;
      } else {
        i++;
        continue;
      }
    }
    f.close();
  }
  remove(rg.infile.c_str());
  return;
}

template <size_t bitset_size>
bool search_match(std::bitset<bitset_size> &ref,
                  std::bitset<bitset_size> *mask1, omp_lock_t *dict_lock,
                  omp_lock_t *read_lock,
                  std::bitset<bitset_size> mask[MAX_READ_LEN][MAX_READ_LEN],
                  uint8_t *read_lengths, bool *remainingreads,
                  std::bitset<bitset_size> *read, bbhashdict *dict, uint32_t &k,
                  bool rev, int shift, int &ref_len,
                  reorder_global<bitset_size> &rg) {
  std::bitset<bitset_size> b;
  uint64_t ull;
  int64_t dictidx[2];    // to store the start and end index (end not inclusive)
                         // in the dict read_id array
  uint64_t startposidx;  // index in startpos
  bool flag = 0;
  for (int l = 0; l < rg.numdict; l++) {
    // check if dictionary index is within bound
    if (!rev) {
      if (dict[l].end + shift >= ref_len) continue;
    } else {
      if (dict[l].end >= ref_len + shift || dict[l].start <= shift) continue;
    }
    b = ref & mask1[l];
    ull = (b >> 2 * dict[l].start).to_ullong();
    startposidx = dict[l].bphf->lookup(ull);
    if (startposidx >= dict[l].numkeys)  // not found
      continue;
    // check if any other thread is modifying same dictpos
    omp_set_lock(&dict_lock[startposidx & 0xFFFFFF]);
    dict[l].findpos(dictidx, startposidx);
    if (dict[l].empty_bin[startposidx])  // bin is empty
    {
      omp_unset_lock(&dict_lock[startposidx & 0xFFFFFF]);
      continue;
    }
    uint64_t ull1 =
        ((read[dict[l].read_id[dictidx[0]]] & mask1[l]) >> 2 * dict[l].start)
            .to_ullong();
    if (ull == ull1)  // checking if ull is actually the key for this bin
    {
      for (int64_t i = dictidx[1] - 1;
           i >= dictidx[0] && i >= dictidx[1] - rg.maxsearch; i--) {
        auto rid = dict[l].read_id[i];
        size_t hamming;
        if (!rev)
          hamming = ((ref ^ read[rid]) &
                     mask[0][rg.max_readlen -
                             std::min<int>(ref_len - shift, read_lengths[rid])])
                        .count();
        else
          hamming =
              ((ref ^ read[rid]) &
               mask[shift][rg.max_readlen -
                           std::min<int>(ref_len + shift, read_lengths[rid])])
                  .count();
        if (hamming <= rg.thresh) {
          omp_set_lock(&read_lock[rid & 0xFFFFFF]);
          if (remainingreads[rid]) {
            remainingreads[rid] = 0;
            k = rid;
            flag = 1;
          }
          omp_unset_lock(&read_lock[rid & 0xFFFFFF]);
          if (flag == 1) break;
        }
      }
    }
    omp_unset_lock(&dict_lock[startposidx & 0xFFFFFF]);
    if (flag == 1) break;
  }
  return flag;
}

template <size_t bitset_size>
void reorder(std::bitset<bitset_size> *read, bbhashdict *dict,
             uint8_t *read_lengths, reorder_global<bitset_size> &rg) {
  omp_lock_t *dict_lock = new omp_lock_t[rg.num_locks];
  omp_lock_t *read_lock = new omp_lock_t[rg.num_locks];
  for (uint j = 0; j < rg.num_locks; j++) {
    omp_init_lock(&dict_lock[j]);
    omp_init_lock(&read_lock[j]);
  }
  std::bitset<bitset_size> mask[MAX_READ_LEN][MAX_READ_LEN];
  generatemasks<bitset_size>(mask, rg.max_readlen, 2);
  std::bitset<bitset_size> mask1[rg.numdict];
  generateindexmasks<bitset_size>(mask1, dict, rg.numdict, 2);
  bool *remainingreads = new bool[rg.numreads];
  std::fill(remainingreads, remainingreads + rg.numreads, 1);

  // we go through remainingreads array from behind as that speeds up deletion
  // from bin arrays

  uint32_t firstread = 0, unmatched[rg.num_thr];
#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    std::string tid_str = std::to_string(tid);
    std::ofstream foutRC(rg.outfileRC + '.' + tid_str, std::ofstream::out);
    std::ofstream foutflag(rg.outfileflag + '.' + tid_str, std::ofstream::out);
    std::ofstream foutpos(rg.outfilepos + '.' + tid_str,
                          std::ofstream::out | std::ios::binary);
    std::ofstream foutorder(rg.outfileorder + '.' + tid_str,
                            std::ofstream::out | std::ios::binary);
    std::ofstream foutorder_s(rg.outfileorder + ".singleton." + tid_str,
                              std::ofstream::out | std::ios::binary);
    std::ofstream foutlength(rg.outfilereadlength + '.' + tid_str,
                             std::ofstream::out | std::ios::binary);

    unmatched[tid] = 0;
    std::bitset<bitset_size> ref, revref, b;

    int64_t first_rid;
    // first_rid represents first read of contig, used for left searching
    int count[4][MAX_READ_LEN];
    int64_t dictidx[2];  // to store the start and end index (end not inclusive)
                         // in the dict read_id array
    uint64_t startposidx;  // index in startpos
    bool flag = 0, done = 0, prev_unmatched = false, left_search_start = false,
         left_search = false;
    // left_search_start - true when left search is starting (to avoid double
    // deletion of first read)
    // left_search - true during left search - needed so that we know when to
    // pick arbitrary read
    int64_t current, prev;
    uint64_t ull;
    int ref_len;
    int64_t ref_pos, cur_read_pos;
    // starting positions of the ref and the current read in the contig, can be
    // negative during left search or due to RC
    // useful for sorting according to starting position in the encoding stage.

    int64_t remainingpos =
        rg.numreads -
        1;  // used for searching next unmatched read when no match is found
#pragma omp critical
    {  // doing initial setup and first read
      current = firstread;
      firstread +=
          rg.numreads / omp_get_num_threads();  // spread out first read equally
      remainingreads[current] = 0;
      unmatched[tid]++;
    }
#pragma omp barrier
    updaterefcount<bitset_size>(read[current], ref, revref, count, true, false,
                                0, read_lengths[current], ref_len, rg);
    cur_read_pos = 0;
    ref_pos = 0;
    first_rid = current;
    prev_unmatched = true;
    prev = current;
    while (!done) {
      // delete the read from the corresponding dictionary bins (unless we are
      // starting left search)
      if (!left_search_start) {
        for (int l = 0; l < rg.numdict; l++) {
          if (read_lengths[current] <= dict[l].end) continue;
          b = read[current] & mask1[l];
          ull = (b >> 2 * dict[l].start).to_ullong();
          startposidx = dict[l].bphf->lookup(ull);
          // check if any other thread is modifying same dictpos
          omp_set_lock(&dict_lock[startposidx & 0xFFFFFF]);
          dict[l].findpos(dictidx, startposidx);
          dict[l].remove(dictidx, startposidx, current);
          omp_unset_lock(&dict_lock[startposidx & 0xFFFFFF]);
        }
      } else {
        left_search_start = false;
      }
      flag = 0;
      uint32_t k;
      for (int shift = 0; shift < rg.maxshift; shift++) {
        // find forward match
        flag = search_match<bitset_size>(ref, mask1, dict_lock, read_lock, mask,
                                         read_lengths, remainingreads, read,
                                         dict, k, false, shift, ref_len, rg);
        if (flag == 1) {
          current = k;
          int ref_len_old = ref_len;
          updaterefcount<bitset_size>(read[current], ref, revref, count, false,
                                      false, shift, read_lengths[current],
                                      ref_len, rg);
          if (!left_search) {
            cur_read_pos = ref_pos + shift;
            ref_pos = cur_read_pos;
          } else {
            cur_read_pos =
                ref_pos + ref_len_old - shift - read_lengths[current];
            ref_pos = ref_pos + ref_len_old - shift - ref_len;
          }
          if (prev_unmatched == true)  // prev read not singleton, write it now
          {
            foutRC << 'd';
            foutorder.write((char *)&prev, sizeof(uint32_t));
            foutflag << 0;  // for unmatched
            int64_t zero = 0;
            foutpos.write((char *)&zero, sizeof(int64_t));
            foutlength.write((char *)&read_lengths[prev], sizeof(uint8_t));
          }
          foutRC << (left_search ? 'r' : 'd');
          foutorder.write((char *)&current, sizeof(uint32_t));
          foutflag << 1;  // for matched
          foutpos.write((char *)&cur_read_pos, sizeof(int64_t));
          foutlength.write((char *)&read_lengths[current], sizeof(uint8_t));

          prev_unmatched = false;
          break;
        }

        // find reverse match
        flag = search_match<bitset_size>(
            revref, mask1, dict_lock, read_lock, mask, read_lengths,
            remainingreads, read, dict, k, true, shift, ref_len, rg);
        if (flag == 1) {
          current = k;
          int ref_len_old = ref_len;
          updaterefcount<bitset_size>(read[current], ref, revref, count, false,
                                      true, shift, read_lengths[current],
                                      ref_len, rg);
          if (!left_search) {
            cur_read_pos =
                ref_pos + ref_len_old + shift - read_lengths[current];
            ref_pos = ref_pos + ref_len_old + shift - ref_len;
          } else {
            cur_read_pos = ref_pos - shift;
            ref_pos = cur_read_pos;
          }
          if (prev_unmatched == true)  // prev read not singleton, write it now
          {
            foutRC << 'd';
            foutorder.write((char *)&prev, sizeof(uint32_t));
            foutflag << 0;  // for unmatched
            int64_t zero = 0;
            foutpos.write((char *)&zero, sizeof(int64_t));
            foutlength.write((char *)&read_lengths[prev], sizeof(uint8_t));
          }
          foutRC << (left_search ? 'd' : 'r');
          foutorder.write((char *)&current, sizeof(uint32_t));
          foutflag << 1;  // for matched
          foutpos.write((char *)&cur_read_pos, sizeof(int64_t));
          foutlength.write((char *)&read_lengths[current], sizeof(uint8_t));

          prev_unmatched = false;
          break;
        }

        revref <<= 2;
        ref >>= 2;
      }
      if (flag == 0)  // no match found
      {
        if (!left_search)  // start left search
        {
          left_search = true;
          left_search_start = true;
          // update ref and count with RC of first_read
          updaterefcount<bitset_size>(read[first_rid], ref, revref, count, true,
                                      true, 0, read_lengths[first_rid], ref_len,
                                      rg);
          ref_pos = 0;
          cur_read_pos = 0;
        } else  // left search done, now pick arbitrary read and start new
                // contig
        {
          left_search = false;
          for (int64_t j = remainingpos; j >= 0; j--) {
            if (remainingreads[j] == 1) {
              omp_set_lock(&read_lock[j & 0xffffff]);
              if (remainingreads[j])  // checking again inside critical block
              {
                current = j;
                remainingpos = j - 1;
                remainingreads[j] = 0;
                flag = 1;
                unmatched[tid]++;
              }
              omp_unset_lock(&read_lock[j & 0xffffff]);
              if (flag == 1) break;
            }
          }
          if (flag == 0) {
            if (prev_unmatched ==
                true)  // last read was singleton, write it now
            {
              foutorder_s.write((char *)&prev, sizeof(uint32_t));
            }
            done = 1;  // no reads left
          } else {
            updaterefcount<bitset_size>(read[current], ref, revref, count, true,
                                        false, 0, read_lengths[current],
                                        ref_len, rg);
            ref_pos = 0;
            cur_read_pos = 0;
            if (prev_unmatched == true)  // prev read singleton, write it now
            {
              foutorder_s.write((char *)&prev, sizeof(uint32_t));
            }
            prev_unmatched = true;
            first_rid = current;
            prev = current;
          }
        }
      }
    }  // while (!done) end

    foutRC.close();
    foutorder.close();
    foutflag.close();
    foutpos.close();
    foutorder_s.close();
    foutlength.close();
  }  // parallel end

  delete[] remainingreads;
  delete[] dict_lock;
  delete[] read_lock;
  std::cout << "Reordering done, "
            << std::accumulate(unmatched, unmatched + rg.num_thr, 0)
            << " were unmatched\n";
  return;
}

template <size_t bitset_size>
void writetofile(std::bitset<bitset_size> *read, uint8_t *read_lengths,
                 reorder_global<bitset_size> &rg) {
// convert bitset to string for all num_thr files in parallel
#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    std::string tid_str = std::to_string(tid);
    std::ofstream fout(rg.outfile + '.' + tid_str, std::ofstream::out);
    std::ofstream fout_s(rg.outfile + ".singleton." + tid_str,
                         std::ofstream::out);
    std::ifstream finRC(rg.outfileRC + '.' + tid_str, std::ifstream::in);
    std::ifstream finorder(rg.outfileorder + '.' + tid_str,
                           std::ifstream::in | std::ios::binary);
    std::ifstream finorder_s(rg.outfileorder + ".singleton." + tid_str,
                             std::ifstream::in | std::ios::binary);
    char s[MAX_READ_LEN + 1], s1[MAX_READ_LEN + 1];
    uint32_t current;
    char c;
    while (finRC >> std::noskipws >> c)  // read character by character
    {
      finorder.read((char *)&current, sizeof(uint32_t));
      bitsettostring<bitset_size>(read[current], s, read_lengths[current], rg);
      if (c == 'd') {
        fout << s << "\n";
      } else {
        reverse_complement(s, s1, read_lengths[current], rg.chartorevchar);
        fout << s1 << "\n";
      }
    }
    finorder_s.read((char *)&current, sizeof(uint32_t));
    while (!finorder_s.eof()) {
      bitsettostring<bitset_size>(read[current], s, read_lengths[current], rg);
      fout_s << s << "\n";
      finorder_s.read((char *)&current, sizeof(uint32_t));
    }

    fout.close();
    fout_s.close();
    finRC.close();
    finorder.close();
    finorder_s.close();
  }

  // Now combine the num_thr order files
  std::ofstream fout_s(rg.outfile + ".singleton", std::ofstream::out);
  std::ofstream foutorder_s(rg.outfileorder + ".singleton",
                            std::ofstream::out | std::ios::binary);
  for (int tid = 0; tid < rg.num_thr; tid++) {
    std::string tid_str = std::to_string(tid);
    std::ifstream fin_s(rg.outfile + ".singleton." + tid_str,
                        std::ifstream::in);
    std::ifstream finorder_s(rg.outfileorder + ".singleton." + tid_str,
                             std::ifstream::in | std::ios::binary);

    fout_s << fin_s.rdbuf();  // write entire file
    foutorder_s << finorder_s.rdbuf();

    // clear error flags which occur on rdbuffing empty file
    fout_s.clear();
    foutorder_s.clear();

    fin_s.close();
    finorder_s.close();

    remove((rg.outfile + ".singleton." + tid_str).c_str());
    remove((rg.outfileorder + ".singleton." + tid_str).c_str());
  }
  fout_s.close();
  foutorder_s.close();
  return;
}

template <size_t bitset_size>
void reorder_main(const std::string &temp_dir, int max_readlen, int num_thr) {
  reorder_global<bitset_size> *rg_pointer = new reorder_global<bitset_size>;
  reorder_global<bitset_size> &rg = *rg_pointer;
  rg.basedir = temp_dir;
  rg.infile = rg.basedir + "/input_clean.dna";
  rg.outfile = rg.basedir + "/temp.dna";
  rg.outfileRC = rg.basedir + "/read_rev.txt";
  rg.outfileflag = rg.basedir + "/tempflag.txt";
  rg.outfilepos = rg.basedir + "/temppos.txt";
  rg.outfileorder = rg.basedir + "/read_order.bin";
  rg.outfilereadlength = rg.basedir + "/read_lengths.bin";
  rg.infilenumreads = rg.basedir + "/numreads.bin";

  rg.max_readlen = max_readlen;
  rg.num_thr = num_thr;
  rg.maxshift = rg.max_readlen / 2;
  bbhashdict *dict = new bbhashdict[rg.numdict];
  dict[0].start = rg.max_readlen > 100
                      ? rg.max_readlen / 2 - 32
                      : rg.max_readlen / 2 - rg.max_readlen * 32 / 100;
  dict[0].end = rg.max_readlen / 2 - 1;
  dict[1].start = rg.max_readlen / 2;
  dict[1].end = rg.max_readlen > 100
                    ? rg.max_readlen / 2 - 1 + 32
                    : rg.max_readlen / 2 - 1 + rg.max_readlen * 32 / 100;

  std::ifstream f_numreads(rg.infilenumreads, std::ios::binary);
  f_numreads.read((char *)&rg.numreads, sizeof(uint32_t));

  omp_set_num_threads(rg.num_thr);
  setglobalarrays(rg);
  std::bitset<bitset_size> *read = new std::bitset<bitset_size>[rg.numreads];
  uint8_t *read_lengths = new uint8_t[rg.numreads];
  std::cout << "Reading file: " << rg.infile << std::endl;
  readDnaFile<bitset_size>(read, read_lengths, rg);

  std::cout << "Constructing dictionaries\n";
  constructdictionary<bitset_size>(read, dict, read_lengths, rg.numdict,
                                   rg.numreads, 2, rg.basedir, rg.num_thr);
  std::cout << "Reordering reads\n";
  reorder<bitset_size>(read, dict, read_lengths, rg);
  std::cout << "Writing to file\n";
  writetofile<bitset_size>(read, read_lengths, rg);
  delete[] read;
  delete[] dict;
  delete[] read_lengths;
  delete rg_pointer;
  std::cout << "Done!\n";
}

}  // namespace spring

#endif  // SPRING_REORDER_H_
