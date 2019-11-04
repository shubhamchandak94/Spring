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

#include "util.h"
#include <algorithm>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <fstream>
#include <stdexcept>
#include <string>

#include "omp.h"

namespace spring {

uint32_t read_fastq_block(std::istream *fin, std::string *id_array,
                          std::string *read_array, std::string *quality_array,
                          const uint32_t &num_reads) {
  uint32_t num_done = 0;
  std::string comment;
  for (; num_done < num_reads; num_done++) {
    if (!std::getline(*fin, id_array[num_done])) break;
    remove_CR_from_end(id_array[num_done]);
    if (!std::getline(*fin, read_array[num_done]))
      throw std::runtime_error(
          "Invalid FASTQ file. Number of lines not multiple of 4");
    remove_CR_from_end(read_array[num_done]);
    if (!std::getline(*fin, comment))
      throw std::runtime_error(
          "Invalid FASTQ file. Number of lines not multiple of 4");
    if (!std::getline(*fin, quality_array[num_done]))
      throw std::runtime_error(
          "Invalid FASTQ file. Number of lines not multiple of 4");
    remove_CR_from_end(quality_array[num_done]);
  }
  return num_done;
}

void write_fastq_block(std::ofstream &fout, std::string *id_array,
                       std::string *read_array, std::string *quality_array,
                       const uint32_t &num_reads, const bool preserve_quality,
                       const int &num_thr, const bool &gzip_flag) {
  if (!gzip_flag) {
    for (uint32_t i = 0; i < num_reads; i++) {
      fout << id_array[i] << "\n";
      fout << read_array[i] << "\n";
      if (preserve_quality) {
        fout << "+\n";
        fout << quality_array[i] << "\n";
      }
    }
  } else {
    if (num_reads == 0) return;
    std::string *gzip_compressed = new std::string[num_thr];

    uint64_t *start_read_num = new uint64_t[num_thr];
    uint64_t *end_read_num = new uint64_t[num_thr];
    uint64_t num_reads_per_thread =
        1 + ((num_reads - 1) / num_thr);  // ceiling function
    for (uint32_t i = 0; i < num_thr; i++) {
      if (i == 0)
        start_read_num[i] = 0;
      else
        start_read_num[i] = end_read_num[i - 1];
      if (start_read_num[i] > num_reads) start_read_num[i] = num_reads;
      end_read_num[i] = start_read_num[i] + num_reads_per_thread;
      if (end_read_num[i] > num_reads) end_read_num[i] = num_reads;
    }
#pragma omp parallel num_threads(num_thr)
    {
      int tid = omp_get_thread_num();
      boost::iostreams::filtering_ostream out;
      out.push(boost::iostreams::gzip_compressor());
      out.push(boost::iostreams::back_inserter(gzip_compressed[tid]));

      for (uint64_t i = start_read_num[tid]; i < end_read_num[tid]; i++) {
        out << id_array[i] << "\n";
        out << read_array[i] << "\n";
        if (preserve_quality) {
          out << "+\n";
          out << quality_array[i] << "\n";
        }
      }
      boost::iostreams::close(out);

    }  // end omp parallel
    for (uint32_t i = 0; i < num_thr; i++)
      fout.write(&(gzip_compressed[i][0]), gzip_compressed[i].size());
    delete[] gzip_compressed;
    delete[] start_read_num;
    delete[] end_read_num;
  }
}

void reverse_complement(char *s, char *s1, const int readlen) {
  for (int j = 0; j < readlen; j++)
    s1[j] = chartorevchar[(uint8_t)s[readlen - j - 1]];
  s1[readlen] = '\0';
  return;
}

std::string reverse_complement(const std::string &s, const int readlen) {
  std::string s1;
  s1.resize(readlen);
  for (int j = 0; j < readlen; j++)
    s1[j] = chartorevchar[(uint8_t)s[readlen - j - 1]];
  return s1;
}

void remove_CR_from_end(std::string &str) {
  if (str.size())
    if (str[str.size() - 1] == '\r') str.resize(str.size() - 1);
}

bool is_permutation(uint32_t *order_array, const uint32_t &numreads) {
    bool *seen = new bool[numreads]();
    for (uint32_t i = 0; i < numreads; i++) {
        if (seen[order_array[i]] == true) {
            delete[] seen;
            return false;
        }
        seen[order_array[i]] = true;
    }
    delete[] seen;
    return true;
}

}  // namespace spring
