#ifndef ID_COMPRESSION_H
#define ID_COMPRESSION_H

#include "Arithmetic_stream.h"
#include "sam_block.h"

uint8_t decompress_uint8t(Arithmetic_stream as, stream_model model);
int compress_uint8t(Arithmetic_stream as, stream_model model, uint8_t c);
int compress_rname(Arithmetic_stream as, rname_models models, char *rname);
int decompress_rname(Arithmetic_stream as, rname_models models, char *rname);
#endif
