//
//  compression.c
//  XC_s2fastqIO
//
//  Created by Mikel Hernaez on 12/4/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>
#include "id_compression/include/sam_block.h"

namespace spring {
namespace id_comp {

int print_line(struct sam_line_t *sline, FILE *fs) {
  fprintf(fs, "%s\n", sline->ID);
  return 0;
}

void *compress(void *thread_info) {

  unsigned long long lineCtr = 0;


  struct compressor_info_t info = *((struct compressor_info_t *)thread_info);

  // Allocs the Arithmetic and the I/O stream
  Arithmetic_stream as = alloc_arithmetic_stream(info.mode, info.fcomp);

  // Allocs the different blocks and all the models for the Arithmetic
  sam_block samBlock = alloc_sam_models(info.id_array, info.f_order,
                                        info.numreads);
  char prev_ID[1024] = {0};  // these were static before. That didn't play well
                             // with parallelization
  uint32_t prev_tokens_ptr[1024] = {0};
  while (!load_sam_line(samBlock)) {
    compress_id(as, samBlock->IDs->models, *samBlock->IDs->IDs, prev_ID,
                prev_tokens_ptr);
    ++lineCtr;
  }

  // end the compression
  encoder_last_step(as);

  free_arithmetic_stream(as);
  free_sam_models(samBlock);
  return NULL;
}

void *decompress(void *thread_info) {

  struct compressor_info_t *info = (struct compressor_info_t *)thread_info;

  Arithmetic_stream as = alloc_arithmetic_stream(info->mode, info->fcomp);

  sam_block samBlock = alloc_sam_models(NULL, NULL, 0);

  char prev_ID[1024] = {0};
  uint32_t prev_tokens_ptr[1024] = {0};
  uint32_t prev_tokens_len[1024] = {0};
  struct sam_line_t sline;
  // Decompress the blocks
  for (uint32_t n = 0; n < info->numreads; n++) {
    decompress_id(as, samBlock->IDs->models, sline.ID, prev_ID, prev_tokens_ptr,
                  prev_tokens_len);
    info->id_array[n] = sline.ID;
//    print_line(&sline, info->f_id);
  }
  free_arithmetic_stream(as);
  free_sam_models(samBlock);
  return NULL;
}

} // namespace id_comp
} // namespace spring
