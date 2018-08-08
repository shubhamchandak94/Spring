#include "qvz/include/qv_compressor.h"
#include "qvz/include/qvz.h"

namespace spring {
namespace qvz {

/**
 * Update stats structure used for adaptive arithmetic coding
 * @param stats Pointer to stats structure
 * @param x Symbol to update
 * @param r Rescaling condition (if n > r, rescale all stats)
 */
void update_stats(stream_stats_ptr_t stats, uint32_t x, uint32_t r) {
  uint32_t i = 0;

  stats->counts[x] += stats->step;
  stats->n += stats->step;

  if (stats->n > r) {
    stats->n = 0;
    for (i = 0; i < stats->alphabetCard; ++i) {
      if (stats->counts[i]) {
        stats->counts[i] >>= 1;
        stats->counts[i] += 1;
        stats->n += stats->counts[i];
      }
    }
  }
}

/**
 * Initialize stats structures used for adaptive arithmetic coding based on
 * the number of contexts required to handle the set of conditional quantizers
 * that we have (one context per quantizer)
 */

/* // commented as this depends on qlist etc. (see version below)
stream_stats_ptr_t **initialize_stream_stats(
    struct cond_quantizer_list_t *q_list) {
  stream_stats_ptr_t **s;
  uint32_t i = 0, j = 0, k = 0;

  s = (stream_stats_ptr_t **)calloc(q_list->columns,
                                    sizeof(stream_stats_ptr_t *));

  // Allocate jagged array, one set of stats per column
  for (i = 0; i < q_list->columns; ++i) {
    // And for each column, one set of stats per low/high quantizer per previous
    // context
    s[i] = (stream_stats_ptr_t *)calloc(2 * q_list->input_alphabets[i]->size,
                                        sizeof(stream_stats_ptr_t));

    // Finally each individual stat structure needs to be filled in uniformly
    for (j = 0; j < 2 * q_list->input_alphabets[i]->size; ++j) {
      s[i][j] = (stream_stats_ptr_t)calloc(1, sizeof(struct stream_stats_t));
      s[i][j]->counts = (uint32_t *)calloc(
          q_list->q[i][j]->output_alphabet->size, sizeof(uint32_t));

      // Initialize the quantizer's stats uniformly
      for (k = 0; k < q_list->q[i][j]->output_alphabet->size; k++) {
        s[i][j]->counts[k] = 1;
      }
      s[i][j]->n = q_list->q[i][j]->output_alphabet->size;
      s[i][j]->alphabetCard = q_list->q[i][j]->output_alphabet->size;

      // Step size is 8 counts per symbol seen to speed convergence
      s[i][j]->step = 8;
    }
  }

  return s;
}
*/

stream_stats_ptr_t **initialize_stream_stats(
    uint32_t columns) {
  stream_stats_ptr_t **s;
  uint32_t i = 0, j = 0, k = 0;

  s = (stream_stats_ptr_t **)calloc(columns,
                                    sizeof(stream_stats_ptr_t *));

  // Allocate jagged array, one set of stats per column
  for (i = 0; i < columns; ++i) {
    // And for each column, one set of stats per low/high quantizer per previous
    // context
    s[i] = (stream_stats_ptr_t *)calloc(2 * ALPHABET_SIZE,
                                        sizeof(stream_stats_ptr_t));

    // Finally each individual stat structure needs to be filled in uniformly
    for (j = 0; j < 2 * ALPHABET_SIZE; ++j) {
      s[i][j] = (stream_stats_ptr_t)calloc(1, sizeof(struct stream_stats_t));
      s[i][j]->counts = (uint32_t *)calloc(
          ALPHABET_SIZE, sizeof(uint32_t));

      // Initialize the quantizer's stats uniformly
      for (k = 0; k < ALPHABET_SIZE; k++) {
        s[i][j]->counts[k] = 1;
      }
      s[i][j]->n = ALPHABET_SIZE;
      s[i][j]->alphabetCard = ALPHABET_SIZE;

      // Step size is 8 counts per symbol seen to speed convergence
      s[i][j]->step = 8;
    }
  }

  return s;
}
/**
 * @todo add cluster stats
 */
arithStream initialize_arithStream(FILE *fout, uint8_t decompressor_flag,
                                   struct quality_file_t *info) {
  arithStream as;
  uint32_t i;

  memset(&info->well, 0, sizeof(struct well_state_t));

  if (decompressor_flag) {
    fread(info->well.state, sizeof(uint32_t), 32, fout);
  } else {
    // Initialize WELL state vector with libc rand
    srand((uint32_t)time(0));
    for (i = 0; i < 32; ++i) {
#ifndef DEBUG
      info->well.state[i] = rand();
#else
      info->well.state[i] = 0x55555555;
#endif
    }

    // Write the initial WELL state vector to the file first (fixed size of 32
    // bytes)
    // @todo strictly this needs to be stored in network order because we're
    // interpreting it as a 32 bit int
    // but I am a bit too lazy for that right now
    fwrite(info->well.state, sizeof(uint32_t), 32, fout);
  }

  // Must start at zero
  info->well.n = 0;

  as = (arithStream)calloc(1, sizeof(struct arithStream_t));

  as->cluster_stats =
      (stream_stats_ptr_t)calloc(1, sizeof(struct stream_stats_t));
  as->cluster_stats->step = 8;
  as->cluster_stats->counts =
      (uint32_t *)calloc(info->cluster_count, sizeof(uint32_t));
  as->cluster_stats->alphabetCard = info->cluster_count;
  as->cluster_stats->n = info->cluster_count;

  as->stats = (stream_stats_ptr_t ***)calloc(info->cluster_count,
                                             sizeof(stream_stats_ptr_t **));
  for (i = 0; i < info->cluster_count; ++i) {
//    as->stats[i] = initialize_stream_stats(info->clusters->clusters[i].qlist);
    as->stats[i] = initialize_stream_stats(info->columns);
    as->cluster_stats->counts[i] = 1;
  }

  as->a = initialize_arithmetic_encoder(m_arith);
  as->os = alloc_os_stream(fout, decompressor_flag);

  if (decompressor_flag)
    as->a->t = stream_read_bits(as->os, as->a->m);
  else
    as->a->t = 0;

  return as;
}

qv_compressor initialize_qv_compressor(FILE *fout, uint8_t streamDirection,
                                       struct quality_file_t *info) {
  qv_compressor s;
  s = (struct qv_compressor_t *)calloc(1, sizeof(struct qv_compressor_t));
  s->Quals = initialize_arithStream(fout, streamDirection, info);
  return s;
}

} // namespace qvz
} // namespace spring
