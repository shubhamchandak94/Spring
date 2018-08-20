#include <iostream>
#include <stdexcept>
#include <fstream>
#include <string>
#include <omp.h>
#include "util.h"
#include "bcm/bcm.h"

namespace spring
{

void decompress_short(const std::string &temp_dir, const std::string &outfile_1,
const std::string &outfile_2, const compression_params &cp, const int &num_thr) {


}

void decompress_long(const std::string &temp_dir, const std::string &outfile_1,
const std::string &outfile_2, const compression_params &cp, const int &num_thr) {
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
  uint32_t num_reads_per_chunk = cp.num_reads_per_chunk_long;
  bool paired_end = cp.paired_end;
  bool preserve_id = cp.preserve_id;
  bool preserve_quality = cp.preserve_quality;

  std::ofstream fout[2];
  fout[0].open(outfile_1);
  if(paired_end)
    fout[1].open(outfile_2);

  // Check that we were able to open the output files
  if(!fout[0].is_open())
    throw std::runtime_error("Error opening output file");
  if(paired_end)
    if(!fout[1].is_open())
      throw std::runtime_error("Error opening output file");

  uint64_t num_reads_per_step = (uint64_t)num_thr*num_reads_per_chunk;

  // allocate less if the total number of reads is small
  if(paired_end) {
    if(num_reads_per_step > num_reads/2)
      num_reads_per_step = num_reads/2;
  }
  else {
    if(num_reads_per_step > num_reads)
      num_reads_per_step = num_reads;
  }

  std::string *read_array = new std::string[num_reads_per_step];
  std::string *id_array = new std::string[num_reads_per_step];
  std::string *quality_array;
  if(preserve_quality)
    quality_array = new std::string[num_reads_per_step];
  uint32_t *read_lengths_array = new uint32_t[num_reads_per_step];

  omp_set_num_threads(num_thr);

  bool done = false;
  uint32_t num_reads_done = 0; // denotes number of pairs done for PE
  uint32_t num_chunks_done = 0;
  while(!done) {
    uint32_t num_reads_cur_step = num_reads_per_step;
    if(paired_end) {
      if(num_reads_done + num_reads_cur_step >= num_reads/2) {
        num_reads_cur_step = num_reads/2 - num_reads_done;
        done = true;
      }
    }
    else {
      if(num_reads_done + num_reads_cur_step >= num_reads) {
        num_reads_cur_step = num_reads - num_reads_done;
        done = true;
      }
    }
    if(num_reads_cur_step == 0)
      break;
    for(int j = 0; j < 2; j++) {
      if(j == 1 && !paired_end)
        continue;
      #pragma omp parallel
      {
        uint64_t tid = omp_get_thread_num();
        if(tid*num_reads_per_chunk < num_reads_cur_step) {
          uint32_t num_reads_thr = std::min((uint64_t)num_reads_cur_step, (tid+1)*num_reads_per_chunk) - tid*num_reads_per_chunk;

          // Decompress read lengths file and read into array
          std::string infile_name = infilereadlength[j]+"."+std::to_string(num_chunks_done+tid) + ".bcm";
          std::string outfile_name = infilereadlength[j]+"."+std::to_string(num_chunks_done+tid);
          bcm::bcm_decompress(infile_name.c_str(), outfile_name.c_str());
          std::ifstream fin_readlength(outfile_name, std::ios::binary);
          for(uint32_t i = tid*num_reads_per_chunk; i < tid*num_reads_per_chunk + num_reads_thr; i++)
            fin_readlength.read((char*)&read_lengths_array[i], sizeof(uint32_t));
          fin_readlength.close();

          // Decompress reads
          infile_name = infileread[j] + "." + std::to_string(num_chunks_done+tid);
          bcm::bcm_str_array_decompress(infile_name.c_str(), read_array + tid*num_reads_per_chunk, num_reads_thr, read_lengths_array + tid*num_reads_per_chunk);

          if(preserve_quality) {
            // Decompress qualities
            infile_name = infilequality[j] + "." + std::to_string(num_chunks_done+tid);
            bcm::bcm_str_array_decompress(infile_name.c_str(), quality_array + tid*num_reads_per_chunk, num_reads_thr, read_lengths_array + tid*num_reads_per_chunk);
          }
          if(!preserve_id) {
            // Fill id array with fake ids
            for(uint32_t i = tid*num_reads_per_chunk; i < tid*num_reads_per_chunk + num_reads_thr; i++)
              id_array[i] = "@" + std::to_string(num_reads_done + i + 1) + "/" + std::to_string(j+1);
          }
          else {
            if(j == 1 && paired_id_match) {
              // id match found, so modify id array appropriately
              for(uint32_t i = tid*num_reads_per_chunk; i < tid*num_reads_per_chunk + num_reads_thr; i++)
                modify_id(id_array[i], paired_id_code);
            }
            else {
              // Decompress ids
              infile_name = infileid[j] + "." + std::to_string(num_chunks_done+tid);
              decompress_id_block(infile_name.c_str(), id_array + tid*num_reads_per_chunk, num_reads_thr);
            }
          }
        }
      } // end omp parallel
      write_fastq_block(fout[j], id_array, read_array, quality_array, num_reads_cur_step, preserve_quality);
    }
    num_reads_done += num_reads_cur_step;
    num_chunks_done += num_thr;
  }

  fout[0].close();
  if(paired_end)
    fout[1].close();

  delete[] read_array;
  delete[] id_array;
  if(preserve_quality)
    delete[] quality_array;
  delete[] read_lengths_array;
}

} // namespace spring
