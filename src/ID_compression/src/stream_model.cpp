//
//  stream_model.c
//  XC_s2fastqIO
//
//  Created by Mikel Hernaez on 11/5/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#include "ID_compression/include/stream_model.h"
#include <stdbool.h>
#define DEBUG false

namespace spring {
namespace id_comp {

void free_models_array(stream_model *model_ptr, uint32_t num_models) {
  uint32_t i = 0;
  for (i = 0; i < num_models; i++) {
    free_model(model_ptr[i]);
  }
}
void free_model(stream_model model) {
  free(model->counts);
  if (model->alphabet) {
    free(model->alphabet);
  }
  if (model->alphaExist) {
    free(model->alphaExist);
  }
  if (model->alphaMap) {
    free(model->alphaMap);
  }
  free(model);
}

uint32_t update_model(stream_model model, int32_t x) {
  int32_t i = 0;

  // Update the statistics
  model->counts[x] += model->step, model->n += model->step;

  // Rescale if necessary
  // We rescale by dividing by two all the counts

  if (model->n >= model->rescale) {
    model->n = 0;
    for (i = 0; i < (int32_t)model->alphabetCard; i++) {
      model->counts[i] >>= 1, model->counts[i]++;
      model->n += model->counts[i];
    }
  }

  return 1;
}

void send_value_to_as(Arithmetic_stream as, stream_model model, int32_t x) {
  uint32_t i = 0;

  uint32_t cumCountX_1 = 0, cumCountX = 0;

  register uint32_t *count = model->counts;

  // Compute the cumulative counts of x and x-1
  assert(x < model->alphabetCard);

  for (i = 0; i < x; ++i) {
    // cumCountX_1 += model->counts[i];
    cumCountX_1 += *count, ++count;
  }
  // cumCountX = cumCountX_1 + model->counts[x];
  cumCountX = cumCountX_1 + *count;

  assert(cumCountX_1 < cumCountX);

  // Send value to the arithmetic encoder
  arithmetic_encoder_step(as, cumCountX_1, cumCountX, model->n);
  if (DEBUG) {
    printf("cumCountX: %d, cumCountX_1: %d\n", cumCountX, cumCountX_1);
    // printf("x: %d\n", x);
  }
}

int read_value_from_as(Arithmetic_stream as, stream_model model) {
  uint32_t i = 0;

  uint32_t cumCountX_1 = 0, cumCountX = 0, cumCount = 0, cumRange = 0;

  uint32_t x = 0;

  ////////////// ASSERT ///////////////////////
  //    for (int i = 0; i < model->alphabetCard; i++)
  //        foo += model->counts[i];

  //    assert(foo == model->n);
  ////////////////////////////////

  // Decode the symbol x
  cumRange = arithmetic_get_symbol_range(as, model->n);

  while (cumCount <= cumRange) {
    cumCount += model->counts[x++];
  }
  x--;

  // Update the arithmetic encoder

  // Compute the cumulative counts of x and x-1
  assert(x <= model->alphabetCard);

  for (i = 0; i < x; ++i) {
    cumCountX_1 += model->counts[i];
  }
  cumCountX = cumCountX_1 + model->counts[x];

  assert(cumCountX_1 < cumCountX);

  // update the arithmetic encoder
  arithmetic_decoder_step(as, cumCountX_1, cumCountX, model->n);

  if (DEBUG) {
    printf("cumCountX: %d, cumCountX_1: %d\n", cumCountX, cumCountX_1);
    // printf("x: %d\n", x);
  }
  return x;
}

} // namespace id_comp
} // namespace spring
