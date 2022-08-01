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
#include <boost/filesystem.hpp>
#include <fstream>
#include <stdexcept>
#include <string>

#include "id_compression/include/sam_block.h"
#include "omp.h"
#include "qvz/include/qvz.h"

namespace spring {

uint32_t read_fastq_block(std::istream *fin, std::string *id_array,
                          std::string *read_array, std::string *quality_array,
                          const uint32_t &num_reads, const bool &fasta_flag) {
  uint32_t num_done = 0;
  std::string comment;
  for (; num_done < num_reads; num_done++) {
    if (!std::getline(*fin, id_array[num_done])) break;
    remove_CR_from_end(id_array[num_done]);
    if (!std::getline(*fin, read_array[num_done]))
      throw std::runtime_error(
          "Invalid FASTQ(A) file. Number of lines not multiple of 4(2)");
    remove_CR_from_end(read_array[num_done]);
    if (fasta_flag)
      continue;
    if (!std::getline(*fin, comment))
      throw std::runtime_error(
          "Invalid FASTQ(A) file. Number of lines not multiple of 4(2)");
    if (!std::getline(*fin, quality_array[num_done]))
      throw std::runtime_error(
          "Invalid FASTQ(A) file. Number of lines not multiple of 4(2)");
    remove_CR_from_end(quality_array[num_done]);
  }
  return num_done;
}

void write_fastq_block(std::ofstream &fout, std::string *id_array,
                       std::string *read_array, std::string *quality_array,
                       const uint32_t &num_reads, const bool preserve_quality,
                       const int &num_thr, const bool &gzip_flag,
                       const int &gzip_level) {
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
    for (uint32_t i = 0; i < (uint32_t)num_thr; i++) {
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
      out.push(boost::iostreams::gzip_compressor(boost::iostreams::gzip_params(gzip_level)));
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
    for (uint32_t i = 0; i < (uint32_t)num_thr; i++)
      fout.write(&(gzip_compressed[i][0]), gzip_compressed[i].size());
    delete[] gzip_compressed;
    delete[] start_read_num;
    delete[] end_read_num;
  }
}

void compress_id_block(const char *outfile_name, std::string *id_array,
                       const uint32_t &num_ids) {
  struct id_comp::compressor_info_t comp_info;
  comp_info.numreads = num_ids;
  comp_info.mode = COMPRESSION;
  comp_info.id_array = id_array;
  comp_info.fcomp = fopen(outfile_name, "w");
  if (!comp_info.fcomp) {
    perror(outfile_name);
    throw std::runtime_error("ID compression: File output error");
  }
  id_comp::compress((void *)&comp_info);
  fclose(comp_info.fcomp);
}

void decompress_id_block(const char *infile_name, std::string *id_array,
                         const uint32_t &num_ids) {
  struct id_comp::compressor_info_t comp_info;
  comp_info.numreads = num_ids;
  comp_info.mode = DECOMPRESSION;
  comp_info.id_array = id_array;
  comp_info.fcomp = fopen(infile_name, "r");
  if (!comp_info.fcomp) {
    perror(infile_name);
    throw std::runtime_error("ID compression: File input error");
  }
  id_comp::decompress((void *)&comp_info);
  fclose(comp_info.fcomp);
}

void quantize_quality(std::string *quality_array, const uint32_t &num_lines,
                      char *quantization_table) {
  for (uint32_t i = 0; i < num_lines; i++)
    for (uint32_t j = 0; j < quality_array[i].size(); j++)
      quality_array[i][j] = quantization_table[(uint8_t)quality_array[i][j]];
  return;
}

void quantize_quality_qvz(std::string *quality_array, const uint32_t &num_lines,
                          uint32_t *str_len_array, double qv_ratio) {
  struct qvz::qv_options_t opts;
  opts.verbose = 0;
  opts.stats = 0;
  opts.clusters = 1;
  opts.uncompressed = 0;
  opts.ratio = qv_ratio;
  opts.distortion = DISTORTION_MSE;
  opts.mode = MODE_FIXED;
  size_t max_readlen =
      *(std::max_element(str_len_array, str_len_array + num_lines));
  qvz::encode(&opts, max_readlen, num_lines, quality_array, str_len_array);
}

void generate_illumina_binning_table(char *illumina_binning_table) {
  for (uint8_t i = 0; i <= 33 + 1; i++) illumina_binning_table[i] = 33 + 0;
  for (uint8_t i = 33 + 2; i <= 33 + 9; i++) illumina_binning_table[i] = 33 + 6;
  for (uint8_t i = 33 + 10; i <= 33 + 19; i++)
    illumina_binning_table[i] = 33 + 15;
  for (uint8_t i = 33 + 20; i <= 33 + 24; i++)
    illumina_binning_table[i] = 33 + 22;
  for (uint8_t i = 33 + 25; i <= 33 + 29; i++)
    illumina_binning_table[i] = 33 + 27;
  for (uint8_t i = 33 + 30; i <= 33 + 34; i++)
    illumina_binning_table[i] = 33 + 33;
  for (uint8_t i = 33 + 35; i <= 33 + 39; i++)
    illumina_binning_table[i] = 33 + 37;
  for (uint8_t i = 33 + 40; i <= 127; i++) illumina_binning_table[i] = 33 + 40;
}

void generate_binary_binning_table(char *binary_binning_table,
                                   const unsigned int thr,
                                   const unsigned int high,
                                   const unsigned int low) {
  for (uint8_t i = 0; i < 33 + thr; i++) binary_binning_table[i] = 33 + low;
  for (uint8_t i = 33 + thr; i <= 127; i++) binary_binning_table[i] = 33 + high;
}

// ID patterns
// code 0: no pattern found
// code 1: */1 and */2 where * are same in both
// code 2: * and * where * are same in both
// code 3: * 1:# and * 2:# where * and # are common to both and * contains no
// space (used in new versions)
uint8_t find_id_pattern(const std::string &id_1, const std::string &id_2) {
  if (id_1.length() != id_2.length()) return 0;
  if (id_1 == id_2) return 2;
  size_t len = id_1.length();
  size_t i;
  if (id_1[len - 1] == '1' && id_2[len - 1] == '2') {
    // compare rest
    for (i = 0; i < len - 1; i++)
      if (id_1[i] != id_2[i]) break;
    if (i == len - 1) return 1;
  }
  for (i = 0; i < len; i++) {
    if (id_1[i] != id_2[i]) break;
    if (id_1[i] == ' ') {
      if (i < len - 1 && id_1[i + 1] == '1' && id_2[i + 1] == '2')
        i++;
      else
        break;
    }
  }
  if (i == len) return 3;
  return 0;
}

bool check_id_pattern(const std::string &id_1, const std::string &id_2,
                      const uint8_t paired_id_code) {
  if (id_1.length() != id_2.length()) return false;
  size_t len = id_1.length();
  size_t i;
  switch (paired_id_code) {
    case 1:
      if (id_1[len - 1] == '1' && id_2[len - 1] == '2') {
        // compare rest
        for (i = 0; i < len - 1; i++)
          if (id_1[i] != id_2[i]) break;
        if (i == len - 1) return true;
      }
      break;
    case 2:
      if (id_1 == id_2) return true;
      break;
    case 3:
      for (i = 0; i < len; i++) {
        if (id_1[i] != id_2[i]) break;
        if (id_1[i] == ' ') {
          if (i < len - 1 && id_1[i + 1] == '1' && id_2[i + 1] == '2')
            i++;
          else
            break;
        }
      }
      if (i == len) return true;
      break;
    default:
      throw std::runtime_error("Invalid paired id code.");
  }
  return false;
}

void modify_id(std::string &id, const uint8_t paired_id_code) {
  if (paired_id_code == 2)
    return;
  else if (paired_id_code == 1) {
    id.back() = '2';
    return;
  } else if (paired_id_code == 3) {
    int i = 0;
    while (id[i] != ' ') i++;
    id[i + 1] = '2';
    return;
  }
}

void write_dna_in_bits(const std::string &read, std::ofstream &fout) {
  uint8_t dna2int[128];
  dna2int[(uint8_t)'A'] = 0;
  dna2int[(uint8_t)'C'] = 2; // chosen to align with the bitset representation
  dna2int[(uint8_t)'G'] = 1;
  dna2int[(uint8_t)'T'] = 3;
  uint8_t bitarray[128];
  uint8_t pos_in_bitarray = 0;
  uint16_t readlen = read.size();
  fout.write((char *)&readlen, sizeof(uint16_t));
  for (int i = 0; i < readlen / 4; i++) {
    bitarray[pos_in_bitarray] = 0;
    for (int j = 0; j < 4; j++)
      bitarray[pos_in_bitarray] |= (dna2int[(uint8_t)read[4 * i + j]]<<(2*j));
    pos_in_bitarray++;
  }
  if (readlen % 4 != 0) {
    int i = readlen / 4;
    bitarray[pos_in_bitarray] = 0;
    for (int j = 0; j < readlen % 4; j++)
      bitarray[pos_in_bitarray] |= (dna2int[(uint8_t)read[4 * i + j]]<<(2*j));
    pos_in_bitarray++;
  }
  fout.write((char *)&bitarray[0], pos_in_bitarray);
  return;
}

void read_dna_from_bits(std::string &read, std::ifstream &fin) {
  uint16_t readlen;
  uint8_t bitarray[128];
  const char int2dna[4] = {'A','G','C','T'};
  fin.read((char *)&readlen, sizeof(uint16_t));
  read.resize(readlen);
  uint16_t num_bytes_to_read = ((uint32_t)readlen+4-1)/4;
  fin.read((char*)&bitarray[0],num_bytes_to_read);
  uint8_t pos_in_bitarray = 0;
  for (int i = 0; i < readlen / 4; i++) {
    for (int j = 0; j < 4; j++) {
      read[4 * i + j] = int2dna[bitarray[pos_in_bitarray] & 3];
      bitarray[pos_in_bitarray]>>=2;
    }
    pos_in_bitarray++;
  }
  if (readlen % 4 != 0) {
    int i = readlen / 4;
    for (int j = 0; j < readlen % 4; j++) {
      read[4 * i + j] = int2dna[bitarray[pos_in_bitarray] & 3];
      bitarray[pos_in_bitarray]>>=2;
    }
    pos_in_bitarray++;
  }
}

void write_dnaN_in_bits(const std::string &read, std::ofstream &fout) {
  uint8_t dna2int[128];
  dna2int[(uint8_t)'A'] = 0;
  dna2int[(uint8_t)'C'] = 2; // chosen to align with the bitset representation
  dna2int[(uint8_t)'G'] = 1;
  dna2int[(uint8_t)'T'] = 3;
  dna2int[(uint8_t)'N'] = 4;
  uint8_t bitarray[256];
  uint8_t pos_in_bitarray = 0;
  uint16_t readlen = read.size();
  fout.write((char *)&readlen, sizeof(uint16_t));
  for (int i = 0; i < readlen / 2; i++) {
    bitarray[pos_in_bitarray] = 0;
    for (int j = 0; j < 2; j++)
      bitarray[pos_in_bitarray] |= (dna2int[(uint8_t)read[2 * i + j]]<<(4*j));
    pos_in_bitarray++;
  }
  if (readlen % 2 != 0) {
    int i = readlen / 2;
    bitarray[pos_in_bitarray] = 0;
    for (int j = 0; j < readlen % 2; j++)
      bitarray[pos_in_bitarray] |= (dna2int[(uint8_t)read[2 * i + j]]<<(4*j));
    pos_in_bitarray++;
  }
  fout.write((char *)&bitarray[0], pos_in_bitarray);
  return;
}

void read_dnaN_from_bits(std::string &read, std::ifstream &fin) {
  uint16_t readlen;
  uint8_t bitarray[256];
  const char int2dna[5] = {'A','G','C','T','N'};
  fin.read((char *)&readlen, sizeof(uint16_t));
  read.resize(readlen);
  uint16_t num_bytes_to_read = ((uint32_t)readlen+2-1)/2;
  fin.read((char*)&bitarray[0],num_bytes_to_read);
  uint8_t pos_in_bitarray = 0;
  for (int i = 0; i < readlen / 2; i++) {
    for (int j = 0; j < 2; j++) {
      read[2 * i + j] = int2dna[bitarray[pos_in_bitarray] & 15];
      bitarray[pos_in_bitarray]>>=4;
    }
    pos_in_bitarray++;
  }
  if (readlen % 2 != 0) {
    int i = readlen / 2;
    for (int j = 0; j < readlen % 2; j++) {
      read[2 * i + j] = int2dna[bitarray[pos_in_bitarray] & 15];
      bitarray[pos_in_bitarray]>>=4;
    }
    pos_in_bitarray++;
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

size_t get_directory_size(const std::string &temp_dir) {
  namespace fs = boost::filesystem;
  size_t size = 0;
  fs::path p{temp_dir};
  fs::directory_iterator itr{p};
  for (; itr != fs::directory_iterator{}; ++itr) {
    size += fs::file_size(itr->path());
  }
  return size;
}

// below functions based on code at https://github.com/sean-/postgresql-varint/blob/trunk/src/varint.c
// also on https://github.com/shubhamchandak94/CDTC/blob/master/src/util.cpp

// Used to pack some temporary files efficiently
uint64_t zigzag_encode64(const int64_t n) {
  return (n << 1) ^ (n >> 63);
}

int64_t zigzag_decode64(const uint64_t n) {
  return (n >> 1) ^ -((int64_t)(n & 1));
}

void write_var_int64(const int64_t val, std::ofstream &fout) {
  uint64_t uval = zigzag_encode64(val);
  uint8_t byte;
  while (uval > 127) {
    byte = (uint8_t)(uval & 0x7f) | 0x80;
    fout.write((char*)&byte, sizeof(uint8_t));
    uval >>= 7;
  }
  byte = (uint8_t)(uval & 0x7f);
  fout.write((char*)&byte, sizeof(uint8_t));
}

int64_t read_var_int64(std::ifstream &fin) {
  uint64_t uval = 0;
  uint8_t byte;
  uint8_t shift = 0;
  do {
    fin.read((char*)&byte, sizeof(uint8_t));
    uval |= ((byte & 0x7f) << shift);
    shift += 7;
  } while(byte & 0x80);
  return zigzag_decode64(uval);
}

}  // namespace spring
