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

#include "decompress.h"
#include <omp.h>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include "libbsc/bsc.h"
#include "util.h"

namespace spring {

void set_dec_noise_array(char **dec_noise);

void decompress_short(const std::string &temp_dir, const std::string &outfile_1,
                      const std::string &outfile_2,
                      const compression_params &cp, const int &num_thr,
                      const uint64_t &start_num, const uint64_t &end_num,
                      const bool &gzip_flag) {
  std::string basedir = temp_dir;

  std::string file_seq = basedir + "/read_seq.bin";
  std::string file_flag = basedir + "/read_flag.txt";
  std::string file_pos = basedir + "/read_pos.bin";
  std::string file_pos_pair = basedir + "/read_pos_pair.bin";
  std::string file_RC = basedir + "/read_rev.txt";
  std::string file_RC_pair = basedir + "/read_rev_pair.txt";
  std::string file_readlength = basedir + "/read_lengths.bin";
  std::string file_unaligned = basedir + "/read_unaligned.txt";
  std::string file_noise = basedir + "/read_noise.txt";
  std::string file_noisepos = basedir + "/read_noisepos.bin";
  std::string infilequality[2];
  std::string infileid[2];

  infilequality[0] = basedir + "/quality_1";
  infilequality[1] = basedir + "/quality_2";
  infileid[0] = basedir + "/id_1";
  infileid[1] = basedir + "/id_2";

  uint32_t num_reads = cp.num_reads;
  uint8_t paired_id_code = cp.paired_id_code;
  bool paired_id_match = cp.paired_id_match;
  uint32_t num_reads_per_block = cp.num_reads_per_block;
  bool paired_end = cp.paired_end;
  bool preserve_id = cp.preserve_id;
  bool preserve_quality = cp.preserve_quality;
  bool preserve_order = cp.preserve_order;

  std::string outfile[2] = {outfile_1, outfile_2};
  std::ofstream fout[2];

  for (int j = 0; j < 2; j++) {
    if (j == 1 && !paired_end) continue;
    if (gzip_flag)
      fout[j].open(outfile[j], std::ios::binary);
    else
      fout[j].open(outfile[j]);
  }

  // Check that we were able to open the output files
  if (!fout[0].is_open()) throw std::runtime_error("Error opening output file");
  if (paired_end)
    if (!fout[1].is_open())
      throw std::runtime_error("Error opening output file");

  uint64_t num_reads_per_step = (uint64_t)num_thr * num_reads_per_block;

  // allocate less if the total number of reads is small
  if (paired_end) {
    if (num_reads_per_step > num_reads / 2) num_reads_per_step = num_reads / 2;
  } else {
    if (num_reads_per_step > num_reads) num_reads_per_step = num_reads;
  }

  std::string *read_array_1 = new std::string[num_reads_per_step];
  std::string *read_array_2 = NULL;
  if (paired_end) read_array_2 = new std::string[num_reads_per_step];
  std::string *id_array = new std::string[num_reads_per_step];
  std::string *quality_array = NULL;
  if (preserve_quality) quality_array = new std::string[num_reads_per_step];
  uint32_t *read_lengths_array_1 = new uint32_t[num_reads_per_step];
  uint32_t *read_lengths_array_2 = NULL;
  if (paired_end) read_lengths_array_2 = new uint32_t[num_reads_per_step];
  char **dec_noise;
  dec_noise = new char *[128];
  for (int i = 0; i < 128; i++) dec_noise[i] = new char[128];
  set_dec_noise_array(dec_noise);

  omp_set_num_threads(num_thr);

  // Decompress read_seq and store in a string
  std::string seq;
  int num_thr_e = cp.num_thr;  // number of encoding threads
  decompress_unpack_seq(file_seq, num_thr_e, num_thr);
  for (int tid_e = 0; tid_e < num_thr_e; tid_e++) {
    uint64_t prev_len = seq.size();
    uint64_t file_len;
    std::ifstream in_seq(file_seq + '.' + std::to_string(tid_e));
    in_seq.seekg(0, in_seq.end);
    file_len = in_seq.tellg();
    in_seq.seekg(0);
    seq.resize(prev_len + file_len);
    in_seq.read(&seq[prev_len], file_len);
    in_seq.close();
  }

  bool done = false;
  uint32_t num_blocks_done = start_num / num_reads_per_block;
  uint32_t num_reads_done =
      num_blocks_done *
      num_reads_per_block;  // denotes number of pairs done for PE
  while (!done) {
    uint32_t num_reads_cur_step = num_reads_per_step;
    if (paired_end) {
      if (num_reads_done + num_reads_cur_step >= num_reads / 2) {
        num_reads_cur_step = num_reads / 2 - num_reads_done;
      }
    } else {
      if (num_reads_done + num_reads_cur_step >= num_reads) {
        num_reads_cur_step = num_reads - num_reads_done;
      }
    }
    if (num_reads_cur_step == 0) break;
    for (int j = 0; j < 2; j++) {
      if (j == 1 && !paired_end) continue;
#pragma omp parallel
      {
        uint64_t tid = omp_get_thread_num();
        if (tid * num_reads_per_block < num_reads_cur_step) {
          uint32_t num_reads_thr = std::min((uint64_t)num_reads_cur_step,
                                            (tid + 1) * num_reads_per_block) -
                                   tid * num_reads_per_block;

          if (j == 0) {
            // Read decompression done when j = 0 (even for PE)
            uint32_t block_num = num_blocks_done + tid;

            // Decompress files with bsc
            std::string outfile_bsc =
                file_flag + '.' + std::to_string(block_num);
            std::string infile_bsc = outfile_bsc + ".bsc";
            bsc::BSC_decompress(infile_bsc.c_str(), outfile_bsc.c_str());

            outfile_bsc = file_pos + '.' + std::to_string(block_num);
            infile_bsc = outfile_bsc + ".bsc";
            bsc::BSC_decompress(infile_bsc.c_str(), outfile_bsc.c_str());

            outfile_bsc = file_noise + '.' + std::to_string(block_num);
            infile_bsc = outfile_bsc + ".bsc";
            bsc::BSC_decompress(infile_bsc.c_str(), outfile_bsc.c_str());

            outfile_bsc = file_noisepos + '.' + std::to_string(block_num);
            infile_bsc = outfile_bsc + ".bsc";
            bsc::BSC_decompress(infile_bsc.c_str(), outfile_bsc.c_str());

            outfile_bsc = file_unaligned + '.' + std::to_string(block_num);
            infile_bsc = outfile_bsc + ".bsc";
            bsc::BSC_decompress(infile_bsc.c_str(), outfile_bsc.c_str());

            outfile_bsc = file_readlength + '.' + std::to_string(block_num);
            infile_bsc = outfile_bsc + ".bsc";
            bsc::BSC_decompress(infile_bsc.c_str(), outfile_bsc.c_str());

            outfile_bsc = file_RC + '.' + std::to_string(block_num);
            infile_bsc = outfile_bsc + ".bsc";
            bsc::BSC_decompress(infile_bsc.c_str(), outfile_bsc.c_str());

            if (paired_end) {
              outfile_bsc = file_pos_pair + '.' + std::to_string(block_num);
              infile_bsc = outfile_bsc + ".bsc";
              bsc::BSC_decompress(infile_bsc.c_str(), outfile_bsc.c_str());

              outfile_bsc = file_RC_pair + '.' + std::to_string(block_num);
              infile_bsc = outfile_bsc + ".bsc";
              bsc::BSC_decompress(infile_bsc.c_str(), outfile_bsc.c_str());
            }
            // Open files
            std::ifstream f_flag(file_flag + '.' + std::to_string(block_num));
            std::ifstream f_noise(file_noise + '.' + std::to_string(block_num));
            std::ifstream f_noisepos(
                file_noisepos + '.' + std::to_string(block_num),
                std::ios::binary);
            std::ifstream f_pos(file_pos + '.' + std::to_string(block_num),
                                std::ios::binary);
            std::ifstream f_RC(file_RC + '.' + std::to_string(block_num));
            std::ifstream f_unaligned(file_unaligned + '.' +
                                      std::to_string(block_num));
            std::ifstream f_readlength(
                file_readlength + '.' + std::to_string(block_num),
                std::ios::binary);
            std::ifstream f_pos_pair;
            std::ifstream f_RC_pair;
            if (paired_end) {
              f_pos_pair.open(file_pos_pair + '.' + std::to_string(block_num),
                              std::ios::binary);
              f_RC_pair.open(file_RC_pair + '.' + std::to_string(block_num));
            }

            char flag;
            uint64_t pos_1, pos_2, prevpos;
            bool singleton_1, singleton_2;
            char RC_1, RC_2;
            uint16_t rl;
            uint16_t diffpos_16;
            bool first_read_of_block = true;
            for (uint32_t i = tid * num_reads_per_block;
                 i < tid * num_reads_per_block + num_reads_thr; i++) {
              f_flag >> flag;
              f_readlength.read((char *)&rl, sizeof(uint16_t));
              read_lengths_array_1[i] = rl;
              singleton_1 = (flag == '2') || (flag == '4');
              if (!singleton_1) {
                if (preserve_order)
                  f_pos.read((char *)&pos_1, sizeof(uint64_t));
                else {
                  if (first_read_of_block) {
                    // in order non-preserving mode, if first read (or read 1 in
                    // first pair) is a singleton, then all the rest are.
                    first_read_of_block = false;
                    f_pos.read((char *)&pos_1, sizeof(uint64_t));
                    prevpos = pos_1;
                  } else {
                    f_pos.read((char *)&diffpos_16, sizeof(uint16_t));
                    if (diffpos_16 == 65535)
                      f_pos.read((char *)&pos_1, sizeof(uint64_t));
                    else
                      pos_1 = prevpos + diffpos_16;
                    prevpos = pos_1;
                  }
                }
                f_RC >> RC_1;
                std::string read = seq.substr(pos_1, read_lengths_array_1[i]);
                std::string noise;
                uint16_t noisepos, prevnoisepos = 0;
                std::getline(f_noise, noise);
                for (uint16_t k = 0; k < noise.size(); k++) {
                  f_noisepos.read((char *)&noisepos, sizeof(uint16_t));
                  noisepos += prevnoisepos;
                  read[noisepos] =
                      dec_noise[(uint8_t)read[noisepos]][(uint8_t)noise[k]];
                  prevnoisepos = noisepos;
                }
                if (RC_1 == 'd')
                  read_array_1[i] = read;
                else
                  read_array_1[i] =
                      reverse_complement(read, read_lengths_array_1[i]);
              } else {
                read_array_1[i].resize(read_lengths_array_1[i]);
                f_unaligned.read(&read_array_1[i][0], read_lengths_array_1[i]);
              }

              if (paired_end) {
                int16_t pos_pair_16;
                singleton_2 = (flag == '2') || (flag == '3');
                f_readlength.read((char *)&rl, sizeof(uint16_t));
                read_lengths_array_2[i] = rl;
                if (!singleton_2) {
                  if (flag == '1' || flag == '4') {
                    // reads 1 and 2 encoded independently
                    f_pos.read((char *)&pos_2, sizeof(uint64_t));
                    f_RC >> RC_2;
                  } else {
                    // read 2 pos and RC encoded in terms of read 1
                    char RC_relative;
                    f_pos_pair.read((char *)&pos_pair_16, sizeof(int16_t));
                    pos_2 = pos_1 + pos_pair_16;
                    f_RC_pair >> RC_relative;
                    if (RC_relative == '0')
                      RC_2 = (RC_1 == 'd') ? 'r' : 'd';
                    else
                      RC_2 = (RC_1 == 'd') ? 'd' : 'r';
                  }
                  std::string read = seq.substr(pos_2, read_lengths_array_2[i]);
                  std::string noise;
                  uint16_t noisepos, prevnoisepos = 0;
                  std::getline(f_noise, noise);
                  for (uint16_t k = 0; k < noise.size(); k++) {
                    f_noisepos.read((char *)&noisepos, sizeof(uint16_t));
                    noisepos += prevnoisepos;
                    read[noisepos] =
                        dec_noise[(uint8_t)read[noisepos]][(uint8_t)noise[k]];
                    prevnoisepos = noisepos;
                  }
                  if (RC_2 == 'd')
                    read_array_2[i] = read;
                  else
                    read_array_2[i] =
                        reverse_complement(read, read_lengths_array_2[i]);
                } else {
                  read_array_2[i].resize(read_lengths_array_2[i]);
                  f_unaligned.read(&read_array_2[i][0],
                                   read_lengths_array_2[i]);
                }
              }
            }

            // Close files
            f_flag.close();
            f_noise.close();
            f_noisepos.close();
            f_pos.close();
            f_RC.close();
            f_unaligned.close();
            f_readlength.close();
            if (paired_end) {
              f_pos_pair.close();
              f_RC_pair.close();
            }
          }
          // Decompress ids and quality
          uint32_t *read_lengths_array;
          std::string infile_name;
          if (j == 0)
            read_lengths_array = read_lengths_array_1;
          else
            read_lengths_array = read_lengths_array_2;
          if (preserve_quality) {
            // Decompress qualities
            infile_name =
                infilequality[j] + "." + std::to_string(num_blocks_done + tid);
            bsc::BSC_str_array_decompress(
                infile_name.c_str(), quality_array + tid * num_reads_per_block,
                num_reads_thr, read_lengths_array + tid * num_reads_per_block);
          }
          if (!preserve_id) {
            // Fill id array with fake ids
            for (uint32_t i = tid * num_reads_per_block;
                 i < tid * num_reads_per_block + num_reads_thr; i++)
              id_array[i] = "@" + std::to_string(num_reads_done + i + 1) + "/" +
                            std::to_string(j + 1);
          } else {
            if (j == 1 && paired_id_match) {
              // id match found, so modify id array appropriately
              for (uint32_t i = tid * num_reads_per_block;
                   i < tid * num_reads_per_block + num_reads_thr; i++)
                modify_id(id_array[i], paired_id_code);
            } else {
              // Decompress ids
              infile_name =
                  infileid[j] + "." + std::to_string(num_blocks_done + tid);
              decompress_id_block(infile_name.c_str(),
                                  id_array + tid * num_reads_per_block,
                                  num_reads_thr);
            }
          }
        }
      }  // end omp parallel
      std::string *read_array;
      if (j == 0)
        read_array = read_array_1;
      else
        read_array = read_array_2;
      uint32_t num_reads_cur_step_output = num_reads_cur_step;
      if (num_reads_done + num_reads_cur_step_output >= end_num) {
        num_reads_cur_step_output = end_num - num_reads_done;
        done = true;
      }

      if (num_blocks_done == start_num / num_reads_per_block) {
        // first blocks
        uint32_t shift = start_num % num_reads_per_block;
        write_fastq_block(fout[j], id_array + shift, read_array + shift,
                          quality_array + shift,
                          num_reads_cur_step_output - shift, preserve_quality,
                          num_thr, gzip_flag);
      } else {
        write_fastq_block(fout[j], id_array, read_array, quality_array,
                          num_reads_cur_step_output, preserve_quality, num_thr,
                          gzip_flag);
      }
    }
    num_reads_done += num_reads_cur_step;
    num_blocks_done += num_thr;
  }

  fout[0].close();
  if (paired_end) fout[1].close();

  delete[] read_array_1;
  if (paired_end) delete[] read_array_2;
  delete[] id_array;
  if (preserve_quality) delete[] quality_array;
  delete[] read_lengths_array_1;
  if (paired_end) delete[] read_lengths_array_2;
  for (int i = 0; i < 128; i++) delete[] dec_noise[i];
  delete[] dec_noise;
}

void decompress_long(const std::string &temp_dir, const std::string &outfile_1,
                     const std::string &outfile_2, const compression_params &cp,
                     const int &num_thr, const uint64_t &start_num,
                     const uint64_t &end_num, const bool &gzip_flag) {
  std::string infileread[2];
  std::string infilequality[2];
  std::string infileid[2];
  std::string infilereadlength[2];
  std::string basedir = temp_dir;
  infileread[0] = basedir + "/read_1";
  infileread[1] = basedir + "/read_2";
  infilequality[0] = basedir + "/quality_1";
  infilequality[1] = basedir + "/quality_2";
  infileid[0] = basedir + "/id_1";
  infileid[1] = basedir + "/id_2";
  infilereadlength[0] = basedir + "/readlength_1";
  infilereadlength[1] = basedir + "/readlength_2";

  uint32_t num_reads = cp.num_reads;
  uint8_t paired_id_code = cp.paired_id_code;
  bool paired_id_match = cp.paired_id_match;
  uint32_t num_reads_per_block = cp.num_reads_per_block_long;
  bool paired_end = cp.paired_end;
  bool preserve_id = cp.preserve_id;
  bool preserve_quality = cp.preserve_quality;

  std::string outfile[2] = {outfile_1, outfile_2};
  std::ofstream fout[2];

  for (int j = 0; j < 2; j++) {
    if (j == 1 && !paired_end) continue;
    if (gzip_flag)
      fout[j].open(outfile[j], std::ios::binary);
    else
      fout[j].open(outfile[j]);
  }

  // Check that we were able to open the output files
  if (!fout[0].is_open()) throw std::runtime_error("Error opening output file");
  if (paired_end)
    if (!fout[1].is_open())
      throw std::runtime_error("Error opening output file");

  uint64_t num_reads_per_step = (uint64_t)num_thr * num_reads_per_block;

  // allocate less if the total number of reads is small
  if (paired_end) {
    if (num_reads_per_step > num_reads / 2) num_reads_per_step = num_reads / 2;
  } else {
    if (num_reads_per_step > num_reads) num_reads_per_step = num_reads;
  }

  std::string *read_array = new std::string[num_reads_per_step];
  std::string *id_array = new std::string[num_reads_per_step];
  std::string *quality_array = NULL;
  if (preserve_quality) quality_array = new std::string[num_reads_per_step];
  uint32_t *read_lengths_array = new uint32_t[num_reads_per_step];

  omp_set_num_threads(num_thr);

  bool done = false;

  uint32_t num_blocks_done = start_num / num_reads_per_block;
  uint32_t num_reads_done =
      num_blocks_done *
      num_reads_per_block;  // denotes number of pairs done for PE
  while (!done) {
    uint32_t num_reads_cur_step = num_reads_per_step;
    if (paired_end) {
      if (num_reads_done + num_reads_cur_step >= num_reads / 2) {
        num_reads_cur_step = num_reads / 2 - num_reads_done;
      }
    } else {
      if (num_reads_done + num_reads_cur_step >= num_reads) {
        num_reads_cur_step = num_reads - num_reads_done;
      }
    }
    if (num_reads_cur_step == 0) break;
    for (int j = 0; j < 2; j++) {
      if (j == 1 && !paired_end) continue;
#pragma omp parallel
      {
        uint64_t tid = omp_get_thread_num();
        if (tid * num_reads_per_block < num_reads_cur_step) {
          uint32_t num_reads_thr = std::min((uint64_t)num_reads_cur_step,
                                            (tid + 1) * num_reads_per_block) -
                                   tid * num_reads_per_block;

          // Decompress read lengths file and read into array
          std::string infile_name = infilereadlength[j] + "." +
                                    std::to_string(num_blocks_done + tid) +
                                    ".bsc";
          std::string outfile_name =
              infilereadlength[j] + "." + std::to_string(num_blocks_done + tid);
          bsc::BSC_decompress(infile_name.c_str(), outfile_name.c_str());
          std::ifstream fin_readlength(outfile_name, std::ios::binary);
          for (uint32_t i = tid * num_reads_per_block;
               i < tid * num_reads_per_block + num_reads_thr; i++)
            fin_readlength.read((char *)&read_lengths_array[i],
                                sizeof(uint32_t));
          fin_readlength.close();

          // Decompress reads
          infile_name =
              infileread[j] + "." + std::to_string(num_blocks_done + tid);
          bsc::BSC_str_array_decompress(
              infile_name.c_str(), read_array + tid * num_reads_per_block,
              num_reads_thr, read_lengths_array + tid * num_reads_per_block);

          if (preserve_quality) {
            // Decompress qualities
            infile_name =
                infilequality[j] + "." + std::to_string(num_blocks_done + tid);
            bsc::BSC_str_array_decompress(
                infile_name.c_str(), quality_array + tid * num_reads_per_block,
                num_reads_thr, read_lengths_array + tid * num_reads_per_block);
          }
          if (!preserve_id) {
            // Fill id array with fake ids
            for (uint32_t i = tid * num_reads_per_block;
                 i < tid * num_reads_per_block + num_reads_thr; i++)
              id_array[i] = "@" + std::to_string(num_reads_done + i + 1) + "/" +
                            std::to_string(j + 1);
          } else {
            if (j == 1 && paired_id_match) {
              // id match found, so modify id array appropriately
              for (uint32_t i = tid * num_reads_per_block;
                   i < tid * num_reads_per_block + num_reads_thr; i++)
                modify_id(id_array[i], paired_id_code);
            } else {
              // Decompress ids
              infile_name =
                  infileid[j] + "." + std::to_string(num_blocks_done + tid);
              decompress_id_block(infile_name.c_str(),
                                  id_array + tid * num_reads_per_block,
                                  num_reads_thr);
            }
          }
        }
      }  // end omp parallel
      uint32_t num_reads_cur_step_output = num_reads_cur_step;
      if (num_reads_done + num_reads_cur_step_output >= end_num) {
        num_reads_cur_step_output = end_num - num_reads_done;
        done = true;
      }
      if (num_blocks_done == start_num / num_reads_per_block) {
        // first blocks
        uint32_t shift = start_num % num_reads_per_block;
        write_fastq_block(fout[j], id_array + shift, read_array + shift,
                          quality_array + shift,
                          num_reads_cur_step_output - shift, preserve_quality,
                          num_thr, gzip_flag);
      } else {
        write_fastq_block(fout[j], id_array, read_array, quality_array,
                          num_reads_cur_step_output, preserve_quality, num_thr,
                          gzip_flag);
      }
    }
    num_reads_done += num_reads_cur_step;
    num_blocks_done += num_thr;
  }

  fout[0].close();
  if (paired_end) fout[1].close();

  delete[] read_array;
  delete[] id_array;
  if (preserve_quality) delete[] quality_array;
  delete[] read_lengths_array;
}

void decompress_unpack_seq(const std::string &infile_seq, const int &num_thr_e,
                           const int &num_thr) {
#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    for (int tid_e = tid * num_thr_e / num_thr;
         tid_e < (tid + 1) * num_thr_e / num_thr; tid_e++) {
      bsc::BSC_decompress(
          (infile_seq + '.' + std::to_string(tid_e) + ".bsc").c_str(),
          (infile_seq + '.' + std::to_string(tid_e)).c_str());

      std::ofstream f_seq(infile_seq + '.' + std::to_string(tid_e) + ".tmp");
      std::ifstream in_seq(infile_seq + '.' + std::to_string(tid_e),
                           std::ios::binary);
      std::ifstream in_seq_tail(infile_seq + '.' + std::to_string(tid_e) +
                                ".tail");
      char inttobase[4];
      inttobase[0] = 'A';
      inttobase[1] = 'C';
      inttobase[2] = 'G';
      inttobase[3] = 'T';

      uint8_t dnabin;
      in_seq.read((char *)&dnabin, sizeof(uint8_t));
      while (!in_seq.eof()) {
        f_seq << inttobase[dnabin % 4];
        dnabin /= 4;
        f_seq << inttobase[dnabin % 4];
        dnabin /= 4;
        f_seq << inttobase[dnabin % 4];
        dnabin /= 4;
        f_seq << inttobase[dnabin % 4];
        in_seq.read((char *)&dnabin, sizeof(uint8_t));
      }
      in_seq.close();

      f_seq << in_seq_tail.rdbuf();
      in_seq_tail.close();
      f_seq.close();

      remove((infile_seq + '.' + std::to_string(tid_e)).c_str());
      remove((infile_seq + '.' + std::to_string(tid_e) + ".tail").c_str());
      rename((infile_seq + '.' + std::to_string(tid_e) + ".tmp").c_str(),
             (infile_seq + '.' + std::to_string(tid_e)).c_str());
    }  // for end
  }    // parallel end
}

void set_dec_noise_array(char **dec_noise) {
  dec_noise[(uint8_t)'A'][(uint8_t)'0'] = 'C';
  dec_noise[(uint8_t)'A'][(uint8_t)'1'] = 'G';
  dec_noise[(uint8_t)'A'][(uint8_t)'2'] = 'T';
  dec_noise[(uint8_t)'A'][(uint8_t)'3'] = 'N';
  dec_noise[(uint8_t)'C'][(uint8_t)'0'] = 'A';
  dec_noise[(uint8_t)'C'][(uint8_t)'1'] = 'G';
  dec_noise[(uint8_t)'C'][(uint8_t)'2'] = 'T';
  dec_noise[(uint8_t)'C'][(uint8_t)'3'] = 'N';
  dec_noise[(uint8_t)'G'][(uint8_t)'0'] = 'T';
  dec_noise[(uint8_t)'G'][(uint8_t)'1'] = 'A';
  dec_noise[(uint8_t)'G'][(uint8_t)'2'] = 'C';
  dec_noise[(uint8_t)'G'][(uint8_t)'3'] = 'N';
  dec_noise[(uint8_t)'T'][(uint8_t)'0'] = 'G';
  dec_noise[(uint8_t)'T'][(uint8_t)'1'] = 'C';
  dec_noise[(uint8_t)'T'][(uint8_t)'2'] = 'A';
  dec_noise[(uint8_t)'T'][(uint8_t)'3'] = 'N';
  dec_noise[(uint8_t)'N'][(uint8_t)'0'] = 'A';
  dec_noise[(uint8_t)'N'][(uint8_t)'1'] = 'G';
  dec_noise[(uint8_t)'N'][(uint8_t)'2'] = 'C';
  dec_noise[(uint8_t)'N'][(uint8_t)'3'] = 'T';
}

}  // namespace spring
