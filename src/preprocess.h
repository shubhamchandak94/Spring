#ifndef SPRING_PREPROCESS_H_
#define SPRING_PREPROCESS_H_

#include <string>
#include "util.h"

namespace spring {

void preprocess(std::string &infile_1, std::string &infile_2,
               std::string &temp_dir, bool &paired_end, bool &preserve_id,
               bool &preserve_quality, bool &preserve_order, bool &ill_bin_flag, std::string &quality_compressor, bool &long_flag, compression_params &p);

}  // namespace spring

#endif  // SPRING_PREPROCESS_H_
