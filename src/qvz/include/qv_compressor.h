#ifndef SPRING_QVZ_QV_COMPRESSOR_H_
#define SPRING_QVZ_QV_COMPRESSOR_H_

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <string.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>
 
#include "qvz/include/codebook.h"

#define m_arith 22

#define OS_STREAM_BUF_LEN (4096 * 4096)

#define COMPRESSION 0
#define DECOMPRESSION 1

namespace spring {
namespace qvz {

typedef struct Arithmetic_code_t {
  int32_t scale3;

  uint32_t l;
  uint32_t u;
  uint32_t t;

  uint32_t m;
  uint32_t r;  // Rescaling condition
} * Arithmetic_code;

typedef struct os_stream_t {
  FILE *fp;
  uint8_t *buf;
  uint32_t bufPos;
  uint8_t bitPos;
  uint64_t written;
} * osStream;

typedef struct stream_stats_t {
  uint32_t *counts;
  uint32_t alphabetCard;
  uint32_t step;
  uint32_t n;
} * stream_stats_ptr_t;

typedef struct arithStream_t {
  stream_stats_ptr_t cluster_stats;
  stream_stats_ptr_t ***stats;
  Arithmetic_code a;
  osStream os;
} * arithStream;

typedef struct qv_compressor_t { arithStream Quals; } * qv_compressor;

// Stream interface
struct os_stream_t *alloc_os_stream(FILE *fp, uint8_t in);
void free_os_stream(struct os_stream_t *);
uint8_t stream_read_bit(struct os_stream_t *);
uint32_t stream_read_bits(struct os_stream_t *os, uint8_t len);
void stream_write_bit(struct os_stream_t *, uint8_t);
void stream_write_bits(struct os_stream_t *os, uint32_t dw, uint8_t len);
void stream_finish_byte(struct os_stream_t *);
void stream_write_buffer(struct os_stream_t *);

// Arithmetic ncoder interface
Arithmetic_code initialize_arithmetic_encoder(uint32_t m);
void arithmetic_encoder_step(Arithmetic_code a, stream_stats_ptr_t stats,
                             int32_t x, osStream os);
int encoder_last_step(Arithmetic_code a, osStream os);
uint32_t arithmetic_decoder_step(Arithmetic_code a, stream_stats_ptr_t stats,
                                 osStream is);
uint32_t decoder_last_step(Arithmetic_code a, stream_stats_ptr_t stats);

// Encoding stats management
stream_stats_ptr_t **initialize_stream_stats(
    struct cond_quantizer_list_t *q_list);
void update_stats(stream_stats_ptr_t stats, uint32_t x, uint32_t r);

// Quality value compression interface
void compress_qv(arithStream as, uint32_t x, uint8_t cluster, uint32_t column,
                 uint32_t idx);
void qv_write_cluster(arithStream as, uint8_t cluster);
uint32_t decompress_qv(arithStream as, uint8_t cluster, uint32_t column,
                       uint32_t idx);
uint8_t qv_read_cluster(arithStream as);

qv_compressor initialize_qv_compressor(FILE *fout, uint8_t streamDirection,
                                       struct quality_file_t *info);

uint32_t start_qv_compression(struct quality_file_t *info, FILE *fout,
                              double *dis, FILE *funcompressed);

uint32_t start_qv_compression_lossless(struct quality_file_t *info, FILE *fout);

void start_qv_decompression(FILE *fout, FILE *fin, struct quality_file_t *info,
                            uint16_t *read_lengths);

void start_qv_decompression_lossless(FILE *fin, struct quality_file_t *info, uint16_t *read_lengths);

} // namespace qvz
} // namespace spring

#endif
