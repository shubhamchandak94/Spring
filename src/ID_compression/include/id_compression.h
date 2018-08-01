#ifndef SPRING_ID_COMPRESSION_H
#define SPRING_ID_COMPRESSION_H

#include "ID_compression/include/Arithmetic_stream.h"
#include "ID_compression/include/sam_block.h"

namespace spring {
namespace id_comp {

uint8_t decompress_uint8t(Arithmetic_stream as, stream_model model);
int compress_uint8t(Arithmetic_stream as, stream_model model, uint8_t c);
int compress_rname(Arithmetic_stream as, rname_models models, char *rname);
int decompress_rname(Arithmetic_stream as, rname_models models, char *rname);

} // namespace id_comp
} // namespace spring
#endif
