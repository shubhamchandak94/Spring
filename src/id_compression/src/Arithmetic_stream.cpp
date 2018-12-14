//
//  Arithmetic_stream.c
//  XC_s2fastqIO
//
//  Created by Mikel Hernaez on 11/4/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#include "id_compression/include/Arithmetic_stream.h"

namespace spring {
namespace id_comp {

int file_available = 0;

/************************************
 *
 *  io_stream functions
 *
 ************************************/

/**
 * Allocates a file stream wrapper for the arithmetic encoder, with a given
 * already opened file handle
 */
struct io_stream_t *alloc_io_stream(uint8_t mode, FILE *fp) {
  struct io_stream_t *rtn =
      (struct io_stream_t *)calloc(1, sizeof(struct io_stream_t));

  rtn->mode = mode;

  rtn->buf = (uint8_t *)calloc(IO_STREAM_BUF_LEN + 1, sizeof(uint8_t));

  // char lockname[256];

  //  int rc;

  switch (mode) {
    case COMPRESSION:
      rtn->fp = fp;
      break;
    case DECOMPRESSION:
      rtn->fp = fp;
      fread(rtn->buf, sizeof(uint8_t), IO_STREAM_BUF_LEN, rtn->fp);
      break;

    default:
      break;
  }

  /*if (mode == DECOMPRESSION) {

//        clean_compressed_dir(rtn);
//        file_available = 100;
      rtn->fileCtr = 0;
      open_new_iofile(rtn);
  }
  else{

      clean_compressed_dir(rtn);

      rtn->fileCtr = 0;
  }*/

  rtn->bufPos = 0;
  rtn->bitPos = 0;
  rtn->written = 0;

  return rtn;
}

/**
 * Deallocate the output stream.
 */
void free_os_stream(struct io_stream_t *os) {
  free(os->buf);
  free(os);
}

/**
 * Reads a single bit from the stream
 */
uint8_t stream_read_bit(struct io_stream_t *is) {
  uint8_t rtn = is->buf[is->bufPos] >> 7;

  is->buf[is->bufPos] = is->buf[is->bufPos] << 1;
  is->bitPos += 1;

  if (is->bitPos == 8) {
    is->bitPos = 0;
    is->bufPos += 1;
    if (is->bufPos == IO_STREAM_BUF_LEN) {
      stream_fill_buffer(is);
    }
  }

  return rtn;
}

/**
 * Reads a grouping of bits to be interpreted as a single integer, regardless of
 * length
 * Bits are implicitly written in bit endian order ONE AT A TIME elsewhere in
 * the code,
 * so this must read that way too
 */
uint32_t stream_read_bits(struct io_stream_t *is, uint8_t len) {
  uint32_t rtn = 0;
  int8_t bit;

  for (bit = len - 1; bit >= 0; --bit) {
    rtn |= stream_read_bit(is) << bit;
  }

  return rtn;
}

/**
 * Writes a single bit to the stream
 */
void stream_write_bit(struct io_stream_t *os, uint8_t bit) {
  bit = (bit & 1);
  os->buf[os->bufPos] |= bit;

  os->bitPos += 1;

  if (os->bitPos == 8) {
    os->bitPos = 0;
    os->bufPos += 1;
    if (os->bufPos == IO_STREAM_BUF_LEN) {
      stream_write_buffer(os);
    }
  } else {
    os->buf[os->bufPos] <<= 1;
  }
}

/**
 * Writes a grouping of bits to be interpreted as a single integer and read back
 * the
 * same way. Bits need to be written msb first
 */
void stream_write_bits(struct io_stream_t *os, uint32_t dw, uint8_t len) {
  int8_t bit;

  for (bit = len - 1; bit >= 0; --bit) {
    stream_write_bit(os, (uint8_t)(dw >> bit));
  }
}

/**
 * Finishes the current byte in progress and writes the buffer out
 */
void stream_finish_byte(struct io_stream_t *os) {
  os->buf[os->bufPos] <<= (7 - os->bitPos);
  os->bitPos = 0;
  os->bufPos += 1;
  stream_write_buffer(os);
}

/**
 * Writes a number of bytes to the stream
 */
void stream_write_bytes(struct io_stream_t *is, char *ch, uint32_t len) {
  uint32_t i = 0;

  for (i = 0; i < len; i++) {
    stream_write_bits(is, ch[i], 8);
  }
}

/**
 * Reads a number of bytes from the stream
 */
void stream_read_bytes(struct io_stream_t *is, char *ch, uint32_t len) {
  uint32_t i = 0;

  for (i = 0; i < len; i++) {
    ch[i] = stream_read_bits(is, 8);
  }
}

/**
 *
 */
void stream_read_line(struct io_stream_t *is, char *line, uint32_t len) {
  char ch = 0;
  uint32_t i = 0;
  for (i = 0; i < len; i++) {
    stream_read_bytes(is, &ch, 1);
    if (ch == '\n') {
      return;
    }
    line[i] = ch;
  }
}

/******************************
 *
 *  arithmetic_stream functions
 *
 ********************************/

Arithmetic_stream alloc_arithmetic_stream(uint8_t direction, FILE *fp) {
  Arithmetic_stream a;

  uint32_t m = ARITHMETIC_WORD_LENGTH;

  a = (Arithmetic_stream)calloc(1, sizeof(struct Arithmetic_stream_t));

  a->m = m;
  a->r = 1 << (m - 3);
  a->l = 0;
  a->u = (1 << m) - 1;

  a->ios = alloc_io_stream(direction, fp);

  if (direction == DECOMPRESSION || direction == DOWNLOAD ||
      direction == REMOTE_DECOMPRESSION) {
    // Read the tag
    a->t = stream_read_bits(a->ios, ARITHMETIC_WORD_LENGTH);
  }

  return a;
}

void free_arithmetic_stream(Arithmetic_stream a) {
  free_os_stream(a->ios);
  free(a);
}

/**
 * E1/E2 check for the MSB of the lower and upper regions being the same,
 * indicating that a bit has
 * been determined and must be sent to the output stream
 * E3 checks for upper being 10xxxx... and lower being 01xxxx... indicating that
 * after rescaling the
 * range we are still in the indetermined central region
 */
void arithmetic_encoder_step(Arithmetic_stream a, uint32_t cumCountX_1,
                             uint32_t cumCountX, uint32_t n) {
  uint64_t range = 0;
  uint8_t msbU = 0, msbL = 0, E1_E2 = 0, E3 = 0, smsbL = 0, smsbU = 0;

  // These are actually constants, need to lift a->m out of the struct because
  // it is compile-time constant
  uint32_t msb_shift = a->m - 1;
  uint32_t smsb_shift = a->m - 2;
  uint32_t msb_clear_mask = (1 << msb_shift) - 1;

  range = a->u - a->l + 1;

  // assert(x < stats->alphabetCard);

  // cumCountX_1 = 0;
  // for (i = 0; i < x; ++i) {
  //    cumCountX_1 += stats->counts[i];
  //}
  // cumCountX = cumCountX_1 + stats->counts[x];

  assert(cumCountX_1 < cumCountX);

  a->u = a->l + (uint32_t)((range * cumCountX) / n) - 1;
  a->l = a->l + (uint32_t)((range * cumCountX_1) / n);

  assert(a->l <= a->u);

  // Check the rescaling conditions
  msbL = a->l >> msb_shift;
  msbU = a->u >> msb_shift;
  E1_E2 = (msbL == msbU);
  E3 = 0;

  if (!E1_E2) {
    smsbL = a->l >> smsb_shift;
    smsbU = a->u >> smsb_shift;
    E3 = (smsbL == 0x01 && smsbU == 0x02);
  }

  // While the bounds need rescaling
  while (E1_E2 || E3) {
    if (E1_E2) {
      // We are in one half of the integer range so the next bit is fixed as the
      // MSB
      stream_write_bit(a->ios, msbL);

      // Clear the msb from both bounds and rescale them
      a->l = (a->l & msb_clear_mask) << 1;
      a->u = ((a->u & msb_clear_mask) << 1) + 1;

      // Write any extra bits based on the number of rescalings without an
      // output before now
      while (a->scale3 > 0) {
        stream_write_bit(a->ios, !msbL);
        a->scale3 -= 1;
      }
    } else {  // E3 is true
      a->scale3 += 1;
      a->u = (((a->u << 1) & msb_clear_mask) | (1 << msb_shift)) + 1;
      a->l = (a->l << 1) & msb_clear_mask;
    }

    msbL = a->l >> msb_shift;
    msbU = a->u >> msb_shift;
    E1_E2 = (msbL == msbU);
    E3 = 0;

    if (!E1_E2) {
      smsbL = a->l >> smsb_shift;
      smsbU = a->u >> smsb_shift;
      E3 = (smsbL == 0x01 && smsbU == 0x02);
    }
  }
}

uint64_t encoder_last_step(Arithmetic_stream a) {
  uint8_t msbL = a->l >> (a->m - 1);

  // Write the msb of the tag (l)
  stream_write_bit(a->ios, msbL);

  // write as many !msbL as scale3 left
  while (a->scale3 > 0) {
    stream_write_bit(a->ios, !msbL);
    a->scale3 -= 1;
  }

  // write the rest of the tag (l)
  stream_write_bits(a->ios, a->l, a->m - 1);
  stream_finish_byte(a->ios);

  // Create a dummy 0B file
  //    open_new_iofile(a->ios);
  //    fclose(a->ios->fp);
  //    file_available++;

  return a->ios->written;
}

uint32_t arithmetic_get_symbol_range(Arithmetic_stream a, uint32_t n) {
  uint64_t range, tagGap;

  range = a->u - a->l + 1;
  tagGap = a->t - a->l + 1;

  return (uint32_t)((tagGap * n - 1) / range);

  // while (subRange >= cumCount)
  //  cumCount += stats->counts[k++];

  // x = --k;
}

void arithmetic_decoder_step(Arithmetic_stream a, uint32_t cumCountX_1,
                             uint32_t cumCountX, uint32_t n) {
  uint64_t range = 0;

  uint8_t msbU = 0, msbL = 0, E1_E2 = 0, E3 = 0, smsbL = 0, smsbU = 0;

  // Again, these are actually constants
  uint32_t msb_shift = a->m - 1;
  uint32_t smsb_shift = a->m - 2;
  uint32_t msb_clear_mask = (1 << msb_shift) - 1;

  range = a->u - a->l + 1;

  // tagGap = a->t - a->l + 1;
  // subRange = (uint32_t)((tagGap * stats->n - 1) / range);
  // while (subRange >= cumCount)
  //  cumCount += stats->counts[k++];
  // x = --k;

  // cumCountX_1 = 0;
  // for (i = 0; i < x; ++i) {
  //    cumCountX_1 += stats->counts[i];
  //}
  // cumCountX = cumCountX_1 + stats->counts[x];

  a->u = a->l + (uint32_t)((range * cumCountX) / n) - 1;
  a->l = a->l + (uint32_t)((range * cumCountX_1) / n);

  // Check the rescaling conditions.
  msbL = a->l >> msb_shift;
  msbU = a->u >> msb_shift;

  E1_E2 = (msbL == msbU);
  E3 = 0;

  // If E1 or E2 doen't hold, check E3
  if (!E1_E2) {
    smsbL = a->l >> smsb_shift;
    smsbU = a->u >> smsb_shift;
    E3 = (smsbL == 0x01 && smsbU == 0x02);
  }

  // While any of E conditions hold
  while (E1_E2 || E3) {
    if (E1_E2) {
      a->l = (a->l & msb_clear_mask) << 1;
      a->u = ((a->u & msb_clear_mask) << 1) + 1;
      a->t = ((a->t & msb_clear_mask) << 1) + stream_read_bit(a->ios);
    } else {  // E3 is true
      a->l = (a->l << 1) & msb_clear_mask;
      a->u = (((a->u << 1) & msb_clear_mask) | (1 << msb_shift)) + 1;
      a->t = (((a->t & msb_clear_mask) << 1) ^ (1 << msb_shift)) +
             stream_read_bit(a->ios);
    }

    msbL = a->l >> msb_shift;
    msbU = a->u >> msb_shift;
    E1_E2 = (msbL == msbU);
    E3 = 0;

    if (!E1_E2) {
      smsbL = a->l >> smsb_shift;
      smsbU = a->u >> smsb_shift;
      E3 = (smsbL == 0x01 && smsbU == 0x02);
    }
  }
}

}  // namespace id_comp
}  // namespace spring
