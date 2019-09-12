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

#ifndef SPRING_UTIL_H_
#define SPRING_UTIL_H_

#include <fstream>
#include <string>

namespace spring {

static const char chartorevchar[128] = {
    0, 0,   0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0,
    0, 0,   0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0,
    0, 0,   0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 'T',
    0, 'G', 0, 0, 0, 'C', 0, 0, 0, 0, 0, 0, 'N', 0, 0, 0, 0, 0, 'A', 0, 0, 0,
    0, 0,   0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0,
    0, 0,   0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0};
struct compression_params {
  bool paired_end;
  bool preserve_order;
  bool preserve_quality;
  bool preserve_id;
  bool long_flag;
  bool qvz_flag;
  bool ill_bin_flag;
  bool bin_thr_flag;
  double qvz_ratio;
  unsigned int bin_thr_thr;
  unsigned int bin_thr_high;
  unsigned int bin_thr_low;
  uint32_t num_reads;
  uint32_t num_reads_clean[2];
  uint32_t max_readlen;
  uint8_t paired_id_code;
  bool paired_id_match;
  int num_reads_per_block;
  int num_reads_per_block_long;
  int num_thr;
};

uint32_t read_fastq_block(std::istream *fin, std::string *id_array,
                          std::string *read_array, std::string *quality_array,
                          const uint32_t &num_reads);

void write_fastq_block(std::ofstream &fout, std::string *id_array,
                       std::string *read_array, std::string *quality_array,
                       const uint32_t &num_reads,
                       const bool preserve_quality, const int &num_thr,
                       const bool &gzip_flag);

void compress_id_block(const char *outfile_name, std::string *id_array,
                       const uint32_t &num_ids);

void decompress_id_block(const char *infile_name, std::string *id_array,
                         const uint32_t &num_ids);

void quantize_quality(std::string *quality_array, const uint32_t &num_lines,
                      char *quantization_table);

void quantize_quality_qvz(std::string *quality_array, const uint32_t &num_lines,
                          uint32_t *str_len_array, double qv_ratio);

void generate_illumina_binning_table(char *illumina_binning_table);

void generate_binary_binning_table(char *binary_binning_table,
                                   const unsigned int thr,
                                   const unsigned int high,
                                   const unsigned int low);

uint8_t find_id_pattern(const std::string &id_1, const std::string &id_2);

bool check_id_pattern(const std::string &id_1, const std::string &id_2,
                      const uint8_t paired_id_code);

void modify_id(std::string &id, const uint8_t paired_id_code);

void write_dna_in_bits(const std::string &read, std::ofstream &fout);

void read_dna_from_bits(std::string &read, std::ifstream &fin);
void reverse_complement(char *s, char *s1, const int readlen);

std::string reverse_complement(const std::string &s, const int readlen);

void remove_CR_from_end(std::string &str);

}  // namespace spring

#endif  // SPRING_UTIL_H_
