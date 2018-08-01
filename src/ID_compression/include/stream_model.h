//
//  stream_stats.h
//  XC_s2fastqIO
//
//  Created by Mikel Hernaez on 11/5/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#ifndef SPRING_XC_s2fastqIO_stream_stats_h
#define SPRING_XC_s2fastqIO_stream_stats_h

#include "ID_compression/include/Arithmetic_stream.h"

namespace spring {
namespace id_comp {

typedef struct stream_model_t {
  uint32_t *counts;
  uint32_t alphabetCard;
  uint32_t step;
  uint32_t n;
  uint32_t rescale;

  // variables for the ppm0 implementation
  int32_t *alphabet;
  int32_t *alphaMap;
  uint8_t *alphaExist;
} * stream_model;

uint32_t update_model(stream_model stats, int32_t x);
void send_value_to_as(Arithmetic_stream as, stream_model model, int32_t x);
int read_value_from_as(Arithmetic_stream as, stream_model model);

void free_model(stream_model model);
void free_models_array(stream_model *model_ptr, uint32_t num_models);

} // namespace id_comp
} // namespace spring

#endif
