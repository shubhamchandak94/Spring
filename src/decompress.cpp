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

  std::ofstream fout[2];
  fout[0].open(outfile_1);
  if(cp.paired_end)
    fout[1].open(outfile_2);
  uint32_t num_reads = cp.num_reads;
  uint8_t paired_id_code = cp.paired_id_code;
  bool paired_id_match = cp.paired_id_match;
  uint32_t num_reads_per_chunk = cp.num_reads_per_chunk_long;

  // Check that we were able to open the output files
  if(!fout[0].is_open())
    throw std::runtime_error("Error opening output file");
  if(cp.paired_end)
    if(!fout[1].is_open())
      throw std::runtime_error("Error opening output file");

}

} // namespace spring
