#ifndef SPRING_REORDER_COMPRESS_STREAMS_H_
#define SPRING_REORDER_COMPRESS_STREAMS_H_

#include <string>
#include "util.h"

namespace spring {

void reorder_compress_streams(const std::string &temp_dir,
                              const compression_params &cp);

}  // namespace spring

#endif  // SPRING_REORDER_COMPRESS_STREAMS_H_
