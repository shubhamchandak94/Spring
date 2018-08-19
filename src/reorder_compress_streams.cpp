#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <cstdio>
#include <cstdlib> // abs
#include <omp.h>

#include "reorder_compress_streams.h"
#include "params.h"
#include "bcm/bcm.h"
#include "util.h"

namespace spring {

void reorder_compress_streams (const std::string &temp_dir, const compression_params &cp) {
  if (cp.preserve_order == false)
    throw std::runtime_error("Not implemented");
  std::string basedir = temp_dir;
  std::string file_flag = basedir + "/read_flag.txt";
  // possible flags for PE (for SE):
  // 0: both reads aligned and distance b/w pair < 32767 (SE read aligned)
  // 1: both reads aligned and distance b/w pair >= 32767
  // 2: both reads unaligned (SE read unaligned)
  // 3: first read aligned, second not
  // 4: first read unaligned, second aligned
  std::string file_pos = basedir + "/read_pos.bin";
  // For order preserving mode PE (SE):
  // Flag 0: Store position of first read (Store pos of read)
  // Flag 1: Store position of first and second read
  // Flag 2: Nothing (Nothing)
  // Flag 3: Store position of first read
  // Flag 4: Store position of second read
  std::string file_pos_pair = basedir + "/read_pos_pair.bin";
  // For PE and flag 0, store the distance b/w PE reads with 2 bytes (signed)
  std::string file_RC = basedir + "/read_rev.txt";
  // PE (SE)
  // Flag 0: Orientation of first read (Orientation of SE read)
  // Flag 1: Orientation of first and second read
  // Flag 2: Nothing (Nothing)
  // Flag 3: Orientation of first read
  // Flag 4: Orientation of second read
  std::string file_RC_pair = basedir + "/read_rev_pair.txt";
  // For PE mode and flag 0, store 0 if PE reads have opposite orientation,
  // otherwise store 1
  std::string file_readlength = basedir + "/read_lengths.bin";
  // store read length with 2 bytes for all reads (for PE, store for both
  // reads by interleaving)
  std::string file_unaligned = basedir + "/read_unaligned.txt";
  // store unaligned reads without any newlines
  std::string file_noise = basedir + "/read_noise.txt";
  // store noise separated by newlines (interleave for PE if flag is 0 or 1)
  std::string file_noisepos = basedir + "/read_noisepos.bin";
  // store noisepos with 2 bytes, no newlines, otherwise similar to noise
  std::string file_order = basedir + "/read_order.bin";

  // load some params
  uint32_t num_reads = cp.num_reads, num_reads_aligned = 0, num_reads_unaligned;
  uint32_t num_reads_by_2 = num_reads/2;
  int num_thr = cp.num_thr;
  bool paired_end = cp.paired_end;

  char *RC_arr = new char [num_reads];
  uint16_t *read_length_arr = new uint16_t [num_reads];
  bool *flag_arr = new bool [num_reads];
  uint64_t *pos_in_noise_arr = new uint64_t [num_reads];
  uint64_t *pos_arr = new uint64_t [num_reads];
  uint16_t *noise_len_arr = new uint16_t [num_reads];

  // read streams for aligned reads
  std::ifstream f_order(file_order, std::ios::binary);
  std::ifstream f_RC(file_RC);
  std::ifstream f_readlength(file_readlength, std::ios::binary);
  std::ifstream f_noise(file_noise);
  std::ifstream f_noisepos(file_noisepos, std::ios::binary);
  std::ifstream f_pos(file_pos, std::ios::binary);
  f_noisepos.seekg(0, f_noisepos.end);
  uint64_t noise_array_size = f_noisepos.tellg()/2;
  f_noisepos.seekg(0, f_noisepos.beg);
  // divide by 2 because we have 2 bytes per noise
  char *noise_arr = new char [noise_array_size];
  uint16_t *noisepos_arr = new uint16_t [noise_array_size];
  char rc, noise_char;
  uint32_t order;
  uint64_t current_pos_noise_arr = 0;
  uint64_t current_pos_noisepos_arr = 0;
  uint64_t pos;
  uint16_t num_noise_in_curr_read;
  uint16_t read_length, noisepos;

  while(f_RC.get(rc)) {
    f_order.read((char*)&order, sizeof(uint32_t));
    f_readlength.read((char*)&read_length, sizeof(uint16_t));
    f_pos.read((char*)&pos, sizeof(uint64_t));
    RC_arr[order] = rc;
    read_length_arr[order] = read_length;
    flag_arr[order] = true; // aligned
    pos_arr[order] = pos;
    pos_in_noise_arr[order] = current_pos_noise_arr;
    num_noise_in_curr_read = 0;
    f_noise.get(noise_char);
    while(noise_char != '\n') {
      noise_arr[current_pos_noise_arr++] = noise_char;
      num_noise_in_curr_read++;
      f_noise.get(noise_char);
    }
    for(uint16_t i = 0; i < num_noise_in_curr_read; i++) {
      f_noisepos.read((char*)&noisepos, sizeof(uint16_t));
      noisepos_arr[current_pos_noisepos_arr] = noisepos;
      current_pos_noisepos_arr++;
    }
    noise_len_arr[order] = num_noise_in_curr_read;
    num_reads_aligned++;
  }
  f_noise.close();
  f_noisepos.close();
  f_RC.close();
  f_pos.close();

  // Now start with unaligned reads
  num_reads_unaligned = num_reads - num_reads_aligned;
  std::ifstream f_unaligned(file_unaligned);
  f_unaligned.seekg(0, f_unaligned.end);
  uint64_t unaligned_array_size = f_unaligned.tellg();
  f_unaligned.seekg(0, f_unaligned.beg);
  char *unaligned_arr = new char[unaligned_array_size];
  f_unaligned.read(unaligned_arr, unaligned_array_size);
  f_unaligned.close();
  uint64_t current_pos_in_unaligned_arr = 0;
  for(uint32_t i = 0; i < num_reads_unaligned; i++) {
    f_order.read((char*)&order, sizeof(uint32_t));
    f_readlength.read((char*)&read_length, sizeof(uint16_t));
    read_length_arr[order] = read_length;
    pos_arr[order] = current_pos_in_unaligned_arr;
    current_pos_in_unaligned_arr += read_length;
    flag_arr[order] = false; // unaligned
  }
  f_order.close();
  f_readlength.close();

  // delete old streams
  remove(file_noise.c_str());
  remove(file_noisepos.c_str());
  remove(file_RC.c_str());
  remove(file_order.c_str());
  remove(file_readlength.c_str());
  remove(file_unaligned.c_str());
  remove(file_pos.c_str());


  // Now generate new streams and compress chunks in parallel
  omp_set_num_threads(num_thr);
  uint32_t num_reads_per_chunk = NUM_READS_PER_CHUNK;
  // this is actually number of read pairs per chunk for PE
  #pragma omp parallel
  {
    uint64_t tid = omp_get_thread_num();
    uint64_t chunk_num = tid;
    bool done = false;
    while(!done) {
      uint64_t start_read_num = chunk_num*num_reads_per_chunk;
      uint64_t end_read_num = (chunk_num + 1)*num_reads_per_chunk;
      if(!paired_end) {
        if(start_read_num > num_reads)
          break;
        if(end_read_num > num_reads) {
          done = true;
          end_read_num = num_reads;
        }
      }
      else {
        if(start_read_num > num_reads_by_2)
          break;
        if(end_read_num > num_reads_by_2) {
          done = true;
          end_read_num = num_reads_by_2;
        }
      }
      // Open files
      std::ofstream f_flag(file_flag+'.'+std::to_string(chunk_num));
      std::ofstream f_noise(file_noise+'.'+std::to_string(chunk_num));
      std::ofstream f_noisepos(file_noisepos+'.'+std::to_string(chunk_num), std::ios::binary);
      std::ofstream f_pos(file_pos+'.'+std::to_string(chunk_num), std::ios::binary);
      std::ofstream f_RC(file_RC+'.'+std::to_string(chunk_num));
      std::ofstream f_unaligned(file_unaligned+'.'+std::to_string(chunk_num));
      std::ofstream f_readlength(file_readlength+'.'+std::to_string(chunk_num), std::ios::binary);
      std::ofstream f_pos_pair;
      std::ofstream f_RC_pair;
      if(paired_end) {
        f_pos_pair.open(file_pos_pair+'.'+std::to_string(chunk_num), std::ios::binary);
        f_RC_pair.open(file_RC_pair+'.'+std::to_string(chunk_num));
      }

      // Write streams
      for(uint64_t i = start_read_num; i < end_read_num; i++) {
        if(!paired_end) {
          f_readlength.write((char*)&read_length_arr[i], sizeof(uint16_t));
          if(flag_arr[i] == true) {
            f_flag << '0';
            f_RC << RC_arr[i];
            f_pos.write((char*)&pos_arr[i], sizeof(uint64_t));
            for(uint16_t j = 0; j < noise_len_arr[i]; j++) {
              f_noise << noise_arr[pos_in_noise_arr[i]+j];
              f_noisepos.write((char*)&noisepos_arr[pos_in_noise_arr[i]+j], sizeof(uint16_t));
            }
            f_noise << "\n";
          }
          else {
            f_flag << '2';
            f_unaligned.write(unaligned_arr + pos_arr[i], read_length_arr[i]);
          }
        }
        else {
          uint64_t i_p = num_reads_by_2 + i; //i_pair
          f_readlength.write((char*)&read_length_arr[i], sizeof(uint16_t));
          f_readlength.write((char*)&read_length_arr[i_p], sizeof(uint16_t));
          int64_t pos_pair = (int64_t)pos_arr[i_p] - (int64_t)pos_arr[i];
          int flag;
          if(flag_arr[i] && flag_arr[i_p] && abs(pos_pair) < 32767)
            flag = 0;
          else if(flag_arr[i] && flag_arr[i_p])
            flag = 1;
          else if(!flag_arr[i] && !flag_arr[i_p])
            flag = 2;
          else if(flag_arr[i] && !flag_arr[i_p])
            flag = 3;
          else if(!flag_arr[i] && flag_arr[i_p])
            flag = 4;
          f_flag << flag;
          if(flag == 0 && paired_end) {
            int16_t pos_pair_16 = (int16_t)pos_pair;
            f_pos_pair.write((char*)&pos_pair_16, sizeof(int16_t));
            if(RC_arr[i] != RC_arr[i_p])
              f_RC_pair << '0';
            else
              f_RC_pair << '1';
          }
          if(flag == 0 || flag == 1 || flag == 3) {
            // read 1 is aligned
            f_pos.write((char*)&pos_arr[i], sizeof(uint64_t));
            for(uint16_t j = 0; j < noise_len_arr[i]; j++) {
              f_noise << noise_arr[pos_in_noise_arr[i]+j];
              f_noisepos.write((char*)&noisepos_arr[pos_in_noise_arr[i]+j], sizeof(uint16_t));
            }
            f_noise << "\n";
            f_RC << RC_arr[i];
          }
          else {
            // read 1 is unaligned
            f_unaligned.write(unaligned_arr + pos_arr[i], read_length_arr[i]);
          }

          if(flag == 0 || flag == 1 || flag == 4) {
            // read 2 is aligned
            for(uint16_t j = 0; j < noise_len_arr[i_p]; j++) {
              f_noise << noise_arr[pos_in_noise_arr[i_p]+j];
              f_noisepos.write((char*)&noisepos_arr[pos_in_noise_arr[i_p]+j], sizeof(uint16_t));
            }
            f_noise << "\n";
            if(flag == 1 || flag == 4) {
              // read 2 is aligned but not paired properly
              f_pos.write((char*)&pos_arr[i_p], sizeof(uint64_t));
              f_RC << RC_arr[i_p];
            }
          }
          else {
            // read 2 is unaligned
            f_unaligned.write(unaligned_arr + pos_arr[i_p], read_length_arr[i_p]);
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
      if(paired_end) {
        f_pos_pair.close();
        f_RC_pair.close();
      }

      // Compress files with bcm and remove uncompressed files
      std::string infile_bcm = file_flag+'.'+std::to_string(chunk_num);
      std::string outfile_bcm = infile_bcm + ".bcm";
      bcm::bcm_compress(infile_bcm.c_str(), outfile_bcm.c_str());
      remove(infile_bcm.c_str());

      // TODO: Test impact of packing pos file into
      // minimum number of bits
      infile_bcm = file_pos+'.'+std::to_string(chunk_num);
      outfile_bcm = infile_bcm + ".bcm";
      bcm::bcm_compress(infile_bcm.c_str(), outfile_bcm.c_str());
      remove(infile_bcm.c_str());

      infile_bcm = file_noise+'.'+std::to_string(chunk_num);
      outfile_bcm = infile_bcm + ".bcm";
      bcm::bcm_compress(infile_bcm.c_str(), outfile_bcm.c_str());
      remove(infile_bcm.c_str());

      infile_bcm = file_noisepos+'.'+std::to_string(chunk_num);
      outfile_bcm = infile_bcm + ".bcm";
      bcm::bcm_compress(infile_bcm.c_str(), outfile_bcm.c_str());
      remove(infile_bcm.c_str());

      infile_bcm = file_unaligned+'.'+std::to_string(chunk_num);
      outfile_bcm = infile_bcm + ".bcm";
      bcm::bcm_compress(infile_bcm.c_str(), outfile_bcm.c_str());
      remove(infile_bcm.c_str());

      infile_bcm = file_readlength+'.'+std::to_string(chunk_num);
      outfile_bcm = infile_bcm + ".bcm";
      bcm::bcm_compress(infile_bcm.c_str(), outfile_bcm.c_str());
      remove(infile_bcm.c_str());

      infile_bcm = file_RC+'.'+std::to_string(chunk_num);
      outfile_bcm = infile_bcm + ".bcm";
      bcm::bcm_compress(infile_bcm.c_str(), outfile_bcm.c_str());
      remove(infile_bcm.c_str());

      if(paired_end) {
        infile_bcm = file_pos_pair+'.'+std::to_string(chunk_num);
        outfile_bcm = infile_bcm + ".bcm";
        bcm::bcm_compress(infile_bcm.c_str(), outfile_bcm.c_str());
        remove(infile_bcm.c_str());

        infile_bcm = file_RC_pair+'.'+std::to_string(chunk_num);
        outfile_bcm = infile_bcm + ".bcm";
        bcm::bcm_compress(infile_bcm.c_str(), outfile_bcm.c_str());
        remove(infile_bcm.c_str());
      }

      chunk_num += num_thr;
    }
  } // end omp parallel

  // deallocate
  delete[] RC_arr;
  delete[] read_length_arr;
  delete[] flag_arr;
  delete[] pos_in_noise_arr;
  delete[] pos_arr;
  delete[] noise_len_arr;
  delete[] noise_arr;
  delete[] noisepos_arr;

  return;
}

} // namespace spring
