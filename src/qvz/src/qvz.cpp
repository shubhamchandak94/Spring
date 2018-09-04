#include "qvz/include/util.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <string>

#include "qvz/include/qvz.h"
#include "qvz/include/cluster.h"
#include "qvz/include/codebook.h"
#include "qvz/include/qv_compressor.h"

namespace spring {
namespace qvz {

/**
 *
 */
void encode(struct qv_options_t *opts, uint32_t max_readlen, uint32_t numreads, std::string *quality_string_array, uint32_t *str_len_array) {
  struct quality_file_t qv_info;
  struct distortion_t *dist;
  struct alphabet_t *alphabet = alloc_alphabet(ALPHABET_SIZE);


  if (opts->distortion == DISTORTION_CUSTOM) {
    dist = gen_custom_distortion(ALPHABET_SIZE, opts->dist_file);
  } else {
    dist = generate_distortion_matrix(ALPHABET_SIZE, opts->distortion);
  }

  qv_info.alphabet = alphabet;
  qv_info.dist = dist;
  qv_info.cluster_count = opts->clusters;
  qv_info.columns = max_readlen;
  qv_info.lines = numreads;
  // from alloc_lines & alloc_blocks - we'll allocate single block
  qv_info.block_count = 1;
  qv_info.blocks = (struct line_block_t *)calloc(qv_info.block_count,
                                                 sizeof(struct line_block_t));
  qv_info.blocks[0].count = qv_info.lines;
  qv_info.blocks[0].quality_array = quality_string_array;
  qv_info.blocks[0].read_lengths = str_len_array;

  // Set up clustering data structures
  qv_info.clusters = alloc_cluster_list(&qv_info);
  qv_info.opts = opts;

  // Then find stats and generate codebooks for each cluster
  calculate_statistics(&qv_info);
  generate_codebooks(&qv_info);

  start_qv_quantization(&qv_info);

  free_blocks(&qv_info);
  free_cluster_list(qv_info.clusters);

}

} // namespace qvz
} // namespace spring
