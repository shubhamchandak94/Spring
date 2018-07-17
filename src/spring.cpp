#include "algorithms/SPRING/spring.h"
#include <boost/filesystem.hpp>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include "algorithms/SPRING/encoder.h"
#include "algorithms/SPRING/pe_encode.h"
#include "algorithms/SPRING/preprocess.h"
#include "algorithms/SPRING/reorder.h"
#include "input/fastq/FastqFileReader.h"
#include "algorithms/SPRING/reorder_compress_quality_id.h"

namespace spring {

void generate_streams_SPRING(
    dsg::input::fastq::FastqFileReader *fastqFileReader1,
    dsg::input::fastq::FastqFileReader *fastqFileReader2, int num_thr,
    bool paired_end) {
  // generate random working directory in the current directory
  boost::filesystem::path p1 = boost::filesystem::unique_path();
  boost::filesystem::create_directory(p1);
  std::string working_dir = p1.native();
  std::cout << "Temporary directory: " << working_dir << "\n";
  int status = preprocess(fastqFileReader1, fastqFileReader2, working_dir,
                          paired_end, true, true);
  if (status != 0) throw std::runtime_error("Bad input file");
  std::ifstream f_meta(working_dir + "/read_meta.txt");
  std::string max_readlen_str;
  std::getline(f_meta, max_readlen_str);
  int max_readlen = std::stoi(max_readlen_str);

  call_reorder(working_dir, max_readlen, num_thr);
  call_encoder(working_dir, max_readlen, num_thr);
  if (paired_end == true) pe_encode_main(working_dir, false);
  fastqFileReader1->seekFromSet(0);
  if (paired_end == true)
       fastqFileReader2->seekFromSet(0);
  reorder_compress_quality_id(working_dir, max_readlen, num_thr,
  paired_end, false, true, true, fastqFileReader1, fastqFileReader2, "bsc",
  8.0);
  return;
}

void call_reorder(const std::string &working_dir, int max_readlen, int num_thr) {
  size_t bitset_size_reorder = (2 * max_readlen - 1) / 64 * 64 + 64;
  switch (bitset_size_reorder) {
    case 64:
      reorder_main<64>(working_dir, max_readlen, num_thr);
      break;
    case 128:
      reorder_main<128>(working_dir, max_readlen, num_thr);
      break;
    case 192:
      reorder_main<192>(working_dir, max_readlen, num_thr);
      break;
    case 256:
      reorder_main<256>(working_dir, max_readlen, num_thr);
      break;
    case 320:
      reorder_main<320>(working_dir, max_readlen, num_thr);
      break;
    case 384:
      reorder_main<384>(working_dir, max_readlen, num_thr);
      break;
    case 448:
      reorder_main<448>(working_dir, max_readlen, num_thr);
      break;
    case 512:
      reorder_main<512>(working_dir, max_readlen, num_thr);
      break;
    default:
      throw std::runtime_error("Wrong bitset size.");
  }
}

void call_encoder(const std::string &working_dir, int max_readlen, int num_thr) {
  size_t bitset_size_encoder = (3 * max_readlen - 1) / 64 * 64 + 64;
  switch (bitset_size_encoder) {
    case 64:
      encoder_main<64>(working_dir, max_readlen, num_thr);
      break;
    case 128:
      encoder_main<128>(working_dir, max_readlen, num_thr);
      break;
    case 192:
      encoder_main<192>(working_dir, max_readlen, num_thr);
      break;
    case 256:
      encoder_main<256>(working_dir, max_readlen, num_thr);
      break;
    case 320:
      encoder_main<320>(working_dir, max_readlen, num_thr);
      break;
    case 384:
      encoder_main<384>(working_dir, max_readlen, num_thr);
      break;
    case 448:
      encoder_main<448>(working_dir, max_readlen, num_thr);
      break;
    case 512:
      encoder_main<512>(working_dir, max_readlen, num_thr);
      break;
    case 576:
      encoder_main<576>(working_dir, max_readlen, num_thr);
      break;
    case 640:
      encoder_main<640>(working_dir, max_readlen, num_thr);
      break;
    case 704:
      encoder_main<704>(working_dir, max_readlen, num_thr);
      break;
    case 768:
      encoder_main<768>(working_dir, max_readlen, num_thr);
      break;
    default:
      throw std::runtime_error("Wrong bitset size.");
  }
}

}  // namespace spring
