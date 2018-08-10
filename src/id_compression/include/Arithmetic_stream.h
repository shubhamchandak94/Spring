//
//  Arithmetic_stream.h
//  XC_s2fastqIO
//
//  Created by Mikel Hernaez on 11/4/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#ifndef SPRING_XC_s2fastqIO_Arithmetic_stream_h
#define SPRING_XC_s2fastqIO_Arithmetic_stream_h

#define IO_STREAM_BUF_LEN 2048 * 2048

#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <inttypes.h>

#include <fcntl.h>
#include <sys/stat.h>

#define COMPRESSION 0
#define DECOMPRESSION 1
#define UPLOAD 2
#define DOWNLOAD 3
#define STREAMING 4
#define REMOTE_DECOMPRESSION 5

#define LOSSLESS 8
#define LOSSY 9

#define ARITHMETIC_WORD_LENGTH 26

#ifndef IDOFILE_PATH_ROOT
#define IDOFILE_PATH_ROOT "/tmp/idoFiles/idoFile."
#endif

namespace spring {
namespace id_comp {

extern int file_available;

struct remote_file_info {
  char host_name[1024];
  char username[1024];
  char filename[1024];
};

typedef struct io_stream_t {
  char filePath[1024];

  FILE *fp;

  uint8_t *buf;
  uint32_t bufPos;
  uint8_t bitPos;
  uint64_t written;

  uint32_t fileCtr;

  uint8_t mode;

} * io_stream;

typedef struct Arithmetic_stream_t {
  int32_t scale3;

  uint32_t l;  // lower limit
  uint32_t u;  // upper limit
  uint32_t t;  // the tag

  uint32_t m;  // size of the arithmetic word
  uint32_t r;  // Rescaling condition

  io_stream ios;  // input/output stream of encoded bits
} * Arithmetic_stream;

// Function Prototypes
struct io_stream_t *alloc_io_stream(uint8_t mode, FILE *fp);
void free_os_stream(struct io_stream_t *os);
uint8_t stream_read_bit(struct io_stream_t *is);
uint32_t stream_read_bits(struct io_stream_t *is, uint8_t len);
void stream_write_bit(struct io_stream_t *os, uint8_t bit);
void stream_write_bits(struct io_stream_t *os, uint32_t dw, uint8_t len);
void stream_finish_byte(struct io_stream_t *os);
void stream_write_buffer(struct io_stream_t *os);

void stream_write_bytes(struct io_stream_t *is, char *ch, uint32_t len);
void stream_read_bytes(struct io_stream_t *is, char *ch, uint32_t len);
void stream_read_line(struct io_stream_t *is, char *line, uint32_t len);

Arithmetic_stream alloc_arithmetic_stream(uint8_t direction, FILE *fp);
void free_arithmetic_stream(Arithmetic_stream a);
void arithmetic_encoder_step(Arithmetic_stream a, uint32_t cumCountX_1,
                             uint32_t cumCountX, uint32_t n);
uint64_t encoder_last_step(Arithmetic_stream a);
void arithmetic_decoder_step(Arithmetic_stream a, uint32_t cumCountX,
                             uint32_t cumCountX_1, uint32_t n);
uint32_t arithmetic_get_symbol_range(Arithmetic_stream a, uint32_t n);

void *download(void *remote_info);
void *upload(void *remote_info);
void *remote_decompression(void *remote_info);

int clean_compressed_dir(struct io_stream_t *ios);
void open_new_iofile(struct io_stream_t *ios);
void stream_fill_buffer(struct io_stream_t *os);

} // namespace id_comp
} // namespace spring

#endif
