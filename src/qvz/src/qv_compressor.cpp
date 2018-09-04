#include "qvz/include/qv_compressor.h"
#include <fstream>
#include <string>

namespace spring {
namespace qvz {

/**
 * Compress a sequence of quality scores including dealing with organization by
 * cluster
 */
void start_qv_quantization(struct quality_file_t *info) {

  uint32_t s = 0, idx = 0;
  uint8_t qv = 0, prev_qv = 0;
  uint32_t cur_readlen;
  struct quantizer_t *q;
  struct cond_quantizer_list_t *qlist;

  uint32_t block_idx, line_idx;
  uint8_t cluster_id;

  char *line;
  symbol_t data;

  // Start compressing the file
  block_idx = 0;
  line_idx = 0;
  while(line_idx < info->blocks[block_idx].count) {
    line = &(info->blocks[block_idx].quality_array[line_idx][0]);
    cur_readlen = info->blocks[block_idx].read_lengths[line_idx];

    cluster_id = 0;
    qlist = info->clusters->clusters[cluster_id].qlist;
    // Select first column's codebook with no left context
    q = choose_quantizer(qlist, &info->well, 0, 0, &idx);

    // Quantize, compress and calculate error simultaneously
    data = line[0] - 33;
    qv = q->q[(uint8_t)data];

    line[0] = qv+33; 

    prev_qv = qv;

    for (s = 1; s < cur_readlen; ++s) {
      q = choose_quantizer(qlist, &info->well, s, prev_qv, &idx);
      data = line[s] - 33;
      qv = q->q[(uint8_t)data];
      line[s] = qv+33;
      prev_qv = qv;
    }

    line_idx += 1;
   }
  return;
}

} // namespace qvz
} // namespace spring
