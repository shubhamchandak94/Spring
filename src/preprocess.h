#ifndef SPRING_PREPROCESS_H_
#define SPRING_PREPROCESS_H_

#include <string>
#include "util.h"

namespace spring {

void preprocess(std::string &infile_1, std::string &infile_2,
               std::string &temp_dir, compression_params &cp);

}  // namespace spring

#endif  // SPRING_PREPROCESS_H_
