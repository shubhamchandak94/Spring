#ifndef SPRING_PE_ENCODE_H_
#define SPRING_PE_ENCODE_H_

#include <string>
#include "util.h"

namespace spring {

void pe_encode(const std::string &temp_dir, const compression_params &cp);

}  // namespace spring

#endif  // SPRING_PE_ENCODE_H_
