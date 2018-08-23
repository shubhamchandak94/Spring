//
//  sam_line.h
//  XC_s2fastqIO
//
//  Created by Mikel Hernaez on 11/5/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#ifndef SPRING_XC_s2fastqIO_sam_line_h
#define SPRING_XC_s2fastqIO_sam_line_h

#include <stdio.h>

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <fcntl.h>
#include <sys/stat.h>

#include <pthread.h>

#include "id_compression/include/stream_model.h"
//#include "pmf.h"
//#include "qv_codebook.h"
//#include "util.h"
//#include "well.h"
//#include "quantizer.h"
//#include "aux_data.h"

#include <fstream>
#include <string>

#define MAX_READ_LENGTH 1024
#define MAX_NUMBER_TOKENS_ID 128

// This limits us to chunks that aren't too big to fit into a modest amount of
// memory at a time
//#define MAX_LINES_PER_BLOCK			1
#define MAX_READS_PER_LINE 1022
#define READ_LINEBUF_LENGTH (MAX_READS_PER_LINE + 2)

// Error codes for reading a line block
#define LF_ERROR_NONE 0
#define LF_ERROR_NOT_FOUND 1
#define LF_ERROR_NO_MEMORY 2
#define LF_ERROR_TOO_LONG 4

#define get_qv_model_index(a, b) ((a & 0xff) << 8 | (b & 0xff))

#define MAX_ALPHA 5000000
#define MAX_CARDINALITY 50000000

namespace spring {
namespace id_comp {

struct sam_line_t {
  char ID[1024];
};

struct compressor_info_t {
  FILE *f_id;
  FILE *fcomp;
  std::string *id_array;
  std::ifstream *f_order;
  uint32_t numreads;
  uint8_t mode;
};

typedef struct id_models_t {
  stream_model *token_type;
  stream_model *integer;
  stream_model *delta;
  stream_model *alpha_len;
  stream_model *alpha_value;
  stream_model *chars;
  stream_model *zero_run;
} * id_models;

enum token_type {
  ID_ALPHA,
  ID_DIGIT,
  ID_CHAR,
  ID_MATCH,
  ID_ZEROS,
  ID_DELTA,
  ID_END
};

// To store the model of the chars both in ref and target
enum BASEPAIR { A, C, G, T, N, O };

/**
 *
 */
typedef struct id_block_t {
  char **IDs;
  id_models models;
  uint32_t block_length;
} * id_block;

/**
 *
 */
typedef struct sam_block_t {
  id_block IDs;
  std::string *id_array;
  std::ifstream *f_order;
  uint32_t numreads;
  uint32_t current_read_number;
  stream_model *codebook_model;
} * sam_block;

// Function Prototypes

enum BASEPAIR char2basepair(char c);
int basepair2char(enum BASEPAIR c);
char bp_complement(char c);
enum token_type uint8t2token(uint8_t tok);

stream_model *initialize_stream_model_flag(uint32_t rescale);
stream_model *initialize_stream_model_pos(uint32_t rescale);
stream_model *initialize_stream_model_pos_alpha(uint32_t rescale);
stream_model *initialize_stream_model_match(uint32_t rescale);
stream_model *initialize_stream_model_snps(uint32_t readLength,
                                           uint32_t rescale);
stream_model *initialize_stream_model_indels(uint32_t readLength,
                                             uint32_t rescale);
stream_model *initialize_stream_model_var(uint32_t readLength,
                                          uint32_t rescale);
stream_model *initialize_stream_model_chars(uint32_t rescale);
void initialize_stream_model_qv(stream_model *s,
                                struct cond_quantizer_list_t *q_list);
stream_model *initialize_stream_model_codebook(uint32_t rescale);

// void alloc_stream_model_qv(qv_block qvBlock);

stream_model *alloc_stream_model_qv(uint32_t read_length,
                                    uint32_t input_alphabet_size,
                                    uint32_t rescale);

sam_block alloc_sam_models(  // Arithmetic_stream as,
    std::string *id_array, std::ifstream *f_order, uint32_t numreads);
void free_sam_models(sam_block sb);
uint32_t load_sam_block(sam_block sb);

id_models alloc_id_models_t();
void free_id_models_t(id_models rtn);

int compress_id(Arithmetic_stream as, id_models models, char *id, char *prev_ID,
                uint32_t *prev_tokens_ptr);
int decompress_id(Arithmetic_stream as, id_models model, char *id,
                  char *prev_ID, uint32_t *prev_tokens_ptr,
                  uint32_t *prev_tokens_len);

int compress_line(Arithmetic_stream as, sam_block samBlock);
int decompress_line(Arithmetic_stream as, sam_block samBlock, FILE *f_id);
void *compress(void *thread_info);
void *decompress(void *thread_info);

uint32_t compute_num_digits(uint32_t x);

uint32_t load_sam_line(sam_block sb);

uint8_t create_most_common_list(sam_block sb);
uint8_t get_most_common_token(char **list, uint32_t list_size, char *aux_field);

}  // namespace id_comp
}  // namespace spring

#endif
