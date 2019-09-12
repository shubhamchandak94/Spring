#ifndef SPRING_CALL_TEMPLATE_FUNCTIONS_H_
#define SPRING_CALL_TEMPLATE_FUNCTIONS_H_

#include <string>
#include "util.h"

namespace spring {

void call_reorder(const std::string &temp_dir, compression_params &cp);

void call_encoder(const std::string &temp_dir, compression_params &cp);

}  // namespace spring

#endif  // SPRING_CALL_TEMPLATE_FUNCTIONS_H_
