#ifndef SPRING_SPRING_H_
#define SPRING_SPRING_H_

#include <string>
#include <util.h>

namespace spring {

void compress(std::string &temp_dir, std::vector<std::string>& infile_vec, std::vector<std::string>& outfile_vec, int &num_thr, bool &pairing_only_flag, bool &no_quality_flag, bool &no_ids_flag, bool &ill_bin_flag, bool &long_flag);

void decompress(std::string &temp_dir, std::vector<std::string>& infile_vec, std::vector<std::string>& outfile_vec, int &num_thr);

void call_reorder(const std::string &temp_dir, compression_params &cp);

void call_encoder(const std::string &temp_dir, compression_params &cp);

std::string random_string( size_t length );

} //namespace spring

#endif  // SPRING_SPRING_H_
