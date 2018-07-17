//
//  reads_stream.c
//  XC_s2fastqIO
//
//  Created by Mikel Hernaez on 11/5/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#include "algorithms/SPRING/ID_compression/include/sam_block.h"

namespace spring {
namespace id_comp {

enum BASEPAIR char2basepair(char c) {
  switch (c) {
    case 'A':
      return A;
    case 'C':
      return C;
    case 'G':
      return G;
    case 'T':
      return T;
    default:
      return N;
  }
}

int basepair2char(enum BASEPAIR c) {
  switch (c) {
    case 0:
      return 'A';
    case 1:
      return 'C';
    case 2:
      return 'G';
    case 3:
      return 'T';
    default:
      return 'N';
  }
}

char bp_complement(char c) {
  switch (c) {
    case 'A':
      return 'T';
    case 'C':
      return 'G';
    case 'G':
      return 'C';
    case 'T':
      return 'A';
    default:
      return c;
  }
}

enum token_type uint8t2token(uint8_t tok) {
  switch (tok) {
    case 0:
      return ID_ALPHA;
    case 1:
      return ID_DIGIT;
    case 2:
      return ID_CHAR;
    case 3:
      return ID_MATCH;
    case 4:
      return ID_ZEROS;
    case 5:
      return ID_DELTA;
    case 6:
      return ID_END;
  }
  printf("uint8t2token error: Passed in invalid number %d\n", tok);
  assert(0);
  return ID_END;
}

//////////////////////////////////////////////////////////////////////////////////////////
//                                                                                      //
//                                                                                      //
//                                  INITIALIZATION //
//                                                                                      //
//                                                                                      //
//////////////////////////////////////////////////////////////////////////////////////////

stream_model* initialize_stream_model_id(uint32_t rescale,
                                         uint32_t context_size,
                                         uint32_t alphabet_card) {
  uint32_t i = 0, j = 0;

  stream_model* s = (stream_model*)calloc(context_size, sizeof(stream_model));

  for (i = 0; i < context_size; i++) {
    s[i] = (stream_model)calloc(1, sizeof(struct stream_model_t));

    // Allocate memory
    s[i]->counts = (uint32_t*)calloc(alphabet_card + 1, sizeof(uint32_t));

    // An extra for the cumcounts
    s[i]->counts += 1;

    s[i]->alphabetCard = alphabet_card;

    s[i]->n = 0;
    for (j = 0; j < alphabet_card; j++) {
      s[i]->counts[j] = 1;
      s[i]->n += s[0]->counts[j];
    }

    // STEP
    s[i]->step = 10;

    // rescale bound
    s[i]->rescale = rescale;
  }

  return s;
}

stream_model* initialize_stream_model_flag(uint32_t rescale) {
  uint32_t i = 0, j = 0;
  uint32_t context_size = 1;
  uint32_t alphabetSize = 1 << 16;

  stream_model* s = (stream_model*)calloc(context_size, sizeof(stream_model));

  for (i = 0; i < context_size; i++) {
    s[i] = (stream_model)calloc(1, sizeof(struct stream_model_t));

    // Allocate memory
    s[i]->counts = (uint32_t*)calloc(alphabetSize + 1, sizeof(uint32_t));

    // An extra for the cumcounts
    s[i]->counts += 1;

    s[i]->alphabetCard = alphabetSize;

    for (j = 0; j < alphabetSize; j++) {
      s[i]->counts[j] = 1;
      s[i]->n++;
    }

    // STEP
    s[i]->step = 8;

    // rescale factor
    s[i]->rescale = rescale;
  }

  return s;
}

stream_model* initialize_stream_model_pos(uint32_t rescale) {
  stream_model* s;

  s = (stream_model*)calloc(1, sizeof(stream_model));

  s[0] = (stream_model)calloc(1, sizeof(struct stream_model_t));

  // Allocate memory
  s[0]->alphabet = (int32_t*)calloc(MAX_CARDINALITY, sizeof(int32_t));
  s[0]->counts = (uint32_t*)calloc(MAX_CARDINALITY, sizeof(uint32_t));
  s[0]->alphaExist = (uint8_t*)calloc(MAX_ALPHA, sizeof(uint8_t));
  s[0]->alphaMap = (int32_t*)calloc(MAX_ALPHA, sizeof(int32_t));

  // Initialize the alphabet assigning -1 to the bin 0 of the model
  s[0]->alphabetCard = 1;
  s[0]->alphabet[0] = -1;
  s[0]->alphaMap[0] = -1;

  s[0]->alphaExist[0] = 1;
  s[0]->counts[0] = 1;
  s[0]->n = 1;

  // STEP
  s[0]->step = 10;

  s[0]->rescale = rescale;

  return s;
}

stream_model* initialize_stream_model_pos_alpha(uint32_t rescale) {
  uint32_t i = 0, j = 0;
  uint32_t context_size = 4;

  stream_model* s = (stream_model*)calloc(context_size, sizeof(stream_model));

  for (i = 0; i < context_size; i++) {
    s[i] = (stream_model)calloc(1, sizeof(struct stream_model_t));

    // Allocate memory
    s[i]->counts = (uint32_t*)calloc(257, sizeof(uint32_t));

    // An extra for the cumcounts
    s[i]->counts += 1;

    s[i]->alphabetCard = 256;

    s[i]->n = 0;
    for (j = 0; j < 256; j++) {
      s[i]->counts[j] = 1;
      s[i]->n += s[0]->counts[j];
    }

    // STEP
    s[i]->step = 10;

    // rescale bound
    s[i]->rescale = rescale;
  }

  return s;
}

stream_model* initialize_stream_model_match(uint32_t rescale) {
  uint32_t i = 0;
  uint32_t context_size = 256;

  stream_model* s = (stream_model*)calloc(context_size, sizeof(stream_model));

  for (i = 0; i < context_size; i++) {
    s[i] = (stream_model)calloc(1, sizeof(struct stream_model_t));

    // Allocate memory
    s[i]->counts = (uint32_t*)calloc(8, sizeof(uint32_t));

    // An extra for the cumcounts
    s[i]->counts += 1;

    s[i]->alphabetCard = 2;

    s[i]->counts[0] = 1;
    s[i]->counts[1] = 1;

    s[i]->n = 2;

    // STEP
    s[i]->step = 1;

    // rescale bound
    s[i]->rescale = rescale;
  }

  return s;
}

stream_model* initialize_stream_model_snps(uint32_t readLength,
                                           uint32_t rescale) {
  stream_model* s;

  uint32_t i = 0;

  s = (stream_model*)calloc(1, sizeof(stream_model));

  s[0] = (stream_model)calloc(1, sizeof(struct stream_model_t));

  // Allocate memory
  s[0]->counts = (uint32_t*)calloc(readLength + 2, sizeof(uint32_t));

  // An extra for the cumcounts
  s[0]->counts += 1;

  s[0]->alphabetCard = readLength;

  s[0]->n = 0;
  for (i = 0; i < readLength; i++) {
    s[0]->counts[i] = 1;
    s[0]->n += s[0]->counts[i];
  }

  // STEP
  s[0]->step = 10;

  // rescale bound
  s[0]->rescale = rescale;

  return s;
}

stream_model* initialize_stream_model_indels(uint32_t readLength,
                                             uint32_t rescale) {
  stream_model* s;

  s = (stream_model*)calloc(1, sizeof(stream_model));

  uint32_t i = 0;

  s[0] = (stream_model)calloc(1, sizeof(struct stream_model_t));

  // Allocate memory
  s[0]->counts = (uint32_t*)calloc(readLength + 2, sizeof(uint32_t));

  // An extra for the cumcounts
  s[0]->counts += 1;

  s[0]->alphabetCard = readLength;

  s[0]->n = 0;
  for (i = 0; i < readLength; i++) {
    s[0]->counts[i] = 1;
    s[0]->n += s[0]->counts[i];
  }

  // STEP
  s[0]->step = 16;

  // rescale bound
  s[0]->rescale = rescale;

  return s;
}

stream_model* initialize_stream_model_var(uint32_t readLength,
                                          uint32_t rescale) {
  stream_model* s;

  uint32_t i = 0, j = 0;

  uint32_t num_models = 0xffff;

  s = (stream_model*)calloc(num_models, sizeof(stream_model));

  for (j = 0; j < num_models; j++) {
    s[j] = (stream_model)calloc(1, sizeof(struct stream_model_t));

    // Allocate memory
    s[j]->counts = (uint32_t*)calloc(readLength + 2, sizeof(uint32_t));

    // An extra for the cumcounts
    s[j]->counts += 1;

    s[j]->alphabetCard = readLength;

    s[j]->n = 0;
    for (i = 0; i < readLength; i++) {
      s[j]->counts[i] = 1;
      s[j]->n += s[0]->counts[i];
    }

    // STEP
    s[j]->step = 10;

    // rescale bound
    s[j]->rescale = rescale;
  }

  return s;
}

stream_model* initialize_stream_model_chars(uint32_t rescale) {
  stream_model* s;

  s = (stream_model*)calloc(6, sizeof(stream_model));

  uint32_t i = 0, j;

  for (j = 0; j < 6; j++) {
    s[j] = (stream_model)calloc(1, sizeof(struct stream_model_t));

    // Allocate memory
    s[j]->counts = (uint32_t*)calloc(16, sizeof(uint32_t));

    // An extra for the cumcounts
    s[j]->counts += 1;

    s[j]->alphabetCard = 5;

    s[j]->n = 0;
    for (i = 0; i < 4; i++) {
      s[j]->counts[i] = (i == j) ? 0 : 8;
      s[j]->n += s[j]->counts[i];
    }
    //{A,C,G,T} to N
    s[j]->counts[4] = 1;
    s[j]->n++;

    // STEP
    s[j]->step = 8;

    // rescale bound
    s[j]->rescale = rescale;
  }

  s[0]->counts[1] += 8;
  s[0]->counts[2] += 8;
  s[0]->n += 16;

  s[1]->counts[0] += 8;
  s[1]->counts[3] += 8;
  s[1]->n += 16;

  s[2]->counts[0] += 8;
  s[2]->counts[3] += 8;
  s[2]->n += 16;

  s[3]->counts[1] += 8;
  s[3]->counts[2] += 8;
  s[3]->n += 16;

  // Allow N to N
  // s[4]->counts[4] = 1, s[4]->n++;

  // We need to increase the total count of the insertion case
  // s[5]->n++;

  return s;
}

/**
 *
 */
id_models alloc_id_models_t() {
  uint32_t rescale = 1 << 20;

  id_models rtn = (id_models)calloc(1, sizeof(struct id_models_t));

  rtn->alpha_len =
      initialize_stream_model_id(rescale, MAX_NUMBER_TOKENS_ID, 256);
  rtn->alpha_value =
      initialize_stream_model_id(rescale, MAX_NUMBER_TOKENS_ID, 256);
  rtn->chars = initialize_stream_model_id(rescale, MAX_NUMBER_TOKENS_ID, 256);
  rtn->integer =
      initialize_stream_model_id(rescale, MAX_NUMBER_TOKENS_ID * 4, 256);
  rtn->delta = initialize_stream_model_id(rescale, MAX_NUMBER_TOKENS_ID, 256);
  rtn->zero_run =
      initialize_stream_model_id(rescale, MAX_NUMBER_TOKENS_ID, 256);
  rtn->token_type =
      initialize_stream_model_id(rescale, MAX_NUMBER_TOKENS_ID, 10);

  return rtn;
}

/**
 *
 */
stream_model* initialize_stream_model_codebook(uint32_t rescale) {
  // It is a byte-based model with the previous byte as context.
  uint32_t i = 0, j = 0;
  uint32_t context_size = 256 * 4;

  stream_model* s = (stream_model*)calloc(context_size, sizeof(stream_model));

  for (i = 0; i < context_size; i++) {
    s[i] = (stream_model)calloc(1, sizeof(struct stream_model_t));

    // Allocate memory
    s[i]->counts = (uint32_t*)calloc(257, sizeof(uint32_t));

    // An extra for the cumcounts
    s[i]->counts += 1;

    s[i]->alphabetCard = 256;

    s[i]->n = 0;
    for (j = 0; j < 256; j++) {
      s[i]->counts[j] = 1;
      s[i]->n += s[i]->counts[j];
    }

    // STEP
    s[i]->step = 1;

    // rescale bound
    s[i]->rescale = rescale;
  }

  return s;
}

} // namespace id_comp
} // namespace spring
