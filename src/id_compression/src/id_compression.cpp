//
//  id_compression.c
//  XC_s2fastqIO
//
//  Created by Mikel Hernaez on 12/10/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

// Compression of the IDs -- work based on the compression of ID in Samcomp by
// Mahonney and Bonfiled (2012)

#include <stdio.h>
#include "id_compression/include/sam_block.h"

namespace spring {
namespace id_comp {

uint32_t compute_num_digits(uint32_t x) {
  // Get the number of digits (We assume readLength < 1000)

  if (x < 10)
    return 1;
  else if (x < 100)
    return 2;
  else if (x < 1000)
    return 3;
  else if (x < 10000)
    return 4;
  else if (x < 100000)
    return 5;
  else if (x < 1000000)
    return 6;
  else if (x < 10000000)
    return 7;
  else if (x < 100000000)
    return 8;
  else if (x < 1000000000)
    return 9;
  else
    return 10;
}

uint8_t decompress_uint8t(Arithmetic_stream as, stream_model model) {
  // Send the value to the Arithmetic Stream
  uint8_t c = read_value_from_as(as, model);

  // Update model
  update_model(model, c);

  return c;
}

int compress_uint8t(Arithmetic_stream as, stream_model model, uint8_t c) {
  // Send the value to the Arithmetic Stream
  send_value_to_as(as, model, c);

  // Update model
  update_model(model, c);

  return 1;
}

int compress_id(Arithmetic_stream as, id_models models, char *id, char *prev_ID,
                uint32_t *prev_tokens_ptr) {
  uint8_t token_len = 0, match_len = 0;
  uint32_t i = 0, k = 0, tmp = 0, token_ctr = 0, digit_value = 0,
           digit_model = 0, prev_digit = 0;
  int delta = 0;

  char *id_ptr = id, *id_ptr_tok = NULL;

  while (*id_ptr != 0) {
    match_len += (*id_ptr == prev_ID[prev_tokens_ptr[token_ctr] + token_len]),
        token_len++;
    id_ptr_tok = id_ptr + 1;

    // Check if the token is a alphabetic word
    if (isalpha(*id_ptr)) {
      while (isalpha(*id_ptr_tok)) {
        // compare with the same token from previous ID
        match_len +=
            (*id_ptr_tok == prev_ID[prev_tokens_ptr[token_ctr] + token_len]),
            token_len++, id_ptr_tok++;
      }
      if (match_len == token_len &&
          !isalpha(prev_ID[prev_tokens_ptr[token_ctr] + token_len])) {
        // The token is the same as last ID
        // Encode a token_type ID_MATCH
        compress_uint8t(as, models->token_type[token_ctr], ID_MATCH);

      } else {
        // Encode a token type ID_ALPHA, the length of the string and the string
        compress_uint8t(as, models->token_type[token_ctr], ID_ALPHA);
        compress_uint8t(as, models->alpha_len[token_ctr], token_len);
        for (k = 0; k < token_len; k++) {
          compress_uint8t(as, models->alpha_value[token_ctr], *(id_ptr + k));
        }
      }

    }
    // check if the token is a run of zeros
    else if (*id_ptr == '0') {
      while (*id_ptr_tok == '0') {
        // compare with the same token from previous ID
        match_len += ('0' == prev_ID[prev_tokens_ptr[token_ctr] + token_len]),
            token_len++, id_ptr_tok++;
      }
      if (match_len == token_len &&
          prev_ID[prev_tokens_ptr[token_ctr] + token_len] != '0') {
        // The token is the same as last ID
        // Encode a token_type ID_MATCH
        compress_uint8t(as, models->token_type[token_ctr], ID_MATCH);

      } else {
        // Encode a token type ID_ZEROS and the length of the zeros
        compress_uint8t(as, models->token_type[token_ctr], ID_ZEROS);
        compress_uint8t(as, models->zero_run[token_ctr], token_len);
      }

    }
    // Check if the token is a number smaller than (1<<32-1)
    else if (isdigit(*id_ptr)) {
      digit_value = (*id_ptr - '0');
      bool prev_token_digit_flag =
          true;  // true if corresponding token in previous read is a digit
      if (*prev_ID != 0) {
        if (isdigit(prev_ID[prev_tokens_ptr[token_ctr] + token_len - 1]) &&
            prev_ID[prev_tokens_ptr[token_ctr] + token_len - 1] != '0')
          prev_digit =
              prev_ID[prev_tokens_ptr[token_ctr] + token_len - 1] - '0';
        else
          prev_token_digit_flag = false;
      }

      if (prev_token_digit_flag && *prev_ID != 0) {
        tmp = 1;
        while (isdigit(prev_ID[prev_tokens_ptr[token_ctr] + tmp]) &&
               prev_digit < (1 << 28)) {
          prev_digit = prev_digit * 10 +
                       (prev_ID[prev_tokens_ptr[token_ctr] + tmp] - '0');
          tmp++;
        }
      }

      while (isdigit(*id_ptr_tok) && digit_value < (1 << 28)) {
        digit_value = digit_value * 10 + (*id_ptr_tok - '0');
        // if (*prev_ID != 0){
        //    prev_digit = prev_digit * 10 + (prev_ID[prev_tokens_ptr[token_ctr]
        //    + token_len] - '0');
        //}
        // compare with the same token from previous ID
        match_len +=
            (*id_ptr_tok == prev_ID[prev_tokens_ptr[token_ctr] + token_len]),
            token_len++, id_ptr_tok++;
      }
      if (prev_token_digit_flag && match_len == token_len &&
          !isdigit(prev_ID[prev_tokens_ptr[token_ctr] + token_len])) {
        // The token is the same as last ID
        // Encode a token_type ID_MATCH
        compress_uint8t(as, models->token_type[token_ctr], ID_MATCH);

      } else if (prev_token_digit_flag &&
                 (delta = (digit_value - prev_digit)) < 256 && delta > 0) {
        compress_uint8t(as, models->token_type[token_ctr], ID_DELTA);
        compress_uint8t(as, models->delta[token_ctr], delta);

      } else {
        // Encode a token type ID_DIGIT and the value (byte-based)
        compress_uint8t(as, models->token_type[token_ctr], ID_DIGIT);
        digit_model = (token_ctr << 2);
        compress_uint8t(as, models->integer[digit_model | 0],
                        (digit_value >> 0) & 0xff);
        compress_uint8t(as, models->integer[digit_model | 1],
                        (digit_value >> 8) & 0xff);
        compress_uint8t(as, models->integer[digit_model | 2],
                        (digit_value >> 16) & 0xff);
        compress_uint8t(as, models->integer[digit_model | 3],
                        (digit_value >> 24) & 0xff);
      }
    } else {
      // compare with the same token from previous ID
      // match_len += (*id_ptr == prev_ID[prev_tokens_ptr[token_ctr]]);

      if (match_len == token_len) {
        // The token is the same as last ID
        // Encode a token_type ID_MATCH
        compress_uint8t(as, models->token_type[token_ctr], ID_MATCH);

      } else {
        // Encode a token type ID_CHAR and the char
        compress_uint8t(as, models->token_type[token_ctr], ID_CHAR);
        compress_uint8t(as, models->chars[token_ctr], *id_ptr);
      }
    }

    prev_tokens_ptr[token_ctr] = i;
    i += token_len;
    id_ptr = id_ptr_tok;
    match_len = 0;
    token_len = 0;
    token_ctr++;
  }
  strcpy(prev_ID, id);
  compress_uint8t(as, models->token_type[token_ctr], ID_END);

  return 1;
}

int decompress_id(Arithmetic_stream as, id_models model, char *id,
                  char *prev_ID, uint32_t *prev_tokens_ptr,
                  uint32_t *prev_tokens_len) {
  uint8_t token_len = 0;
  uint32_t i = 0, k = 0, token_ctr = 0, digit_value = 0;
  uint32_t delta = 0;

  enum token_type tok;

  id[0] = '\0';
  while ((tok = uint8t2token(
              decompress_uint8t(as, model->token_type[token_ctr]))) != ID_END) {
    switch (tok) {
      case ID_MATCH:
        memcpy(id + i, &(prev_ID[prev_tokens_ptr[token_ctr]]),
               prev_tokens_len[token_ctr]);
        token_len = prev_tokens_len[token_ctr];
        break;
      case ID_ALPHA:
        token_len = decompress_uint8t(as, model->alpha_len[token_ctr]);
        for (k = 0; k < token_len; k++) {
          id[i + k] = decompress_uint8t(as, model->alpha_value[token_ctr]);
        }
        break;
      case ID_DIGIT:
        digit_value = 0;
        digit_value |=
            ((decompress_uint8t(as, model->integer[(token_ctr << 2) | 0]) &
              0xff)
             << 0);
        digit_value |=
            ((decompress_uint8t(as, model->integer[(token_ctr << 2) | 1]) &
              0xff)
             << 8);
        digit_value |=
            ((decompress_uint8t(as, model->integer[(token_ctr << 2) | 2]) &
              0xff)
             << 16);
        digit_value |=
            ((decompress_uint8t(as, model->integer[(token_ctr << 2) | 3]) &
              0xff)
             << 24);
        sprintf(id + i, "%u", digit_value);
        token_len = compute_num_digits(digit_value);
        break;
      case ID_DELTA:
        digit_value = 0;
        delta = decompress_uint8t(as, model->delta[token_ctr]);
        memcpy(id + i, &(prev_ID[prev_tokens_ptr[token_ctr]]),
               prev_tokens_len[token_ctr]);
        id[i + prev_tokens_len[token_ctr]] = '\0';
        digit_value = atoi(id + i) + delta;
        sprintf(id + i, "%u", digit_value);
        token_len = compute_num_digits(digit_value);
        break;
      case ID_ZEROS:
        token_len = decompress_uint8t(as, model->zero_run[token_ctr]);
        memset(id + i, '0', token_len);
        break;
      case ID_CHAR:
        id[i] = decompress_uint8t(as, model->chars[token_ctr]);
        token_len = 1;
        break;
      default:
        break;
    }

    prev_tokens_ptr[token_ctr] = i;
    prev_tokens_len[token_ctr] = token_len;
    i += token_len;
    id[i] = '\0';
    token_len = 0;
    token_ctr++;
  }
  id[i] = '\0';
  strcpy(prev_ID, id);
  // for (kk=i;kk<=1024;kk++){
  //    prev_ID[kk]='\0';
  //}
  // id[i++] = '\n';
  // putc('@', fs);
  // fwrite(id, i, sizeof(char), fs);

  return 1;
}

}  // namespace id_comp
}  // namespace spring
