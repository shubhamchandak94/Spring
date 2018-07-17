#ifndef SPRING_REORDER_COMPRESS_QUALITY_ID_H_
#define SPRING_REORDER_COMPRESS_QUALITY_ID_H_

#include <string>
#include "input/fastq/FastqFileReader.h"
#include "algorithms/SPRING/reorder_compress_quality_id.h"

namespace spring {

struct reorder_compress_quality_id_global {
  uint8_t paired_id_code;
  bool preserve_order, paired_end, preserve_quality, preserve_id;
  std::string quality_mode;
  double quality_ratio;
  uint32_t numreads, numreads_by_2;
  int max_readlen, num_thr;
  char illumina_binning_table[128];

  std::string infile_id_1;
  std::string infile_id_2;
  std::string infile_order;
  std::string outfile_order;
  std::string infilenumreads;
  std::string basedir;
};

void generate_order(reorder_compress_quality_id_global &rg);
// generate reordering information for the two separate files (pairs) from
// read_order.bin

void reorder_quality(dsg::input::fastq::FastqFileReader *fastqFileReader1,
                     dsg::input::fastq::FastqFileReader *fastqFileReader2,
                     reorder_compress_quality_id_global &rg);
void reorder_id(reorder_compress_quality_id_global &rg);

namespace qvz {
void encode(FILE *fout, struct qv_options_t *opts, uint32_t max_readlen,
            uint32_t numreads, char *quality_array, uint8_t *read_lengths,
            std::string &infile_order, uint64_t startpos);
} // namespace qvz

void illumina_binning(char *quality, uint8_t readlen, reorder_compress_quality_id_global &rg);
void generate_illumina_binning_table(reorder_compress_quality_id_global &rg);

int reorder_compress_quality_id(std::string &working_dir, int max_readlen,
                                int num_thr, bool paired_end,
                                bool preserve_order, bool preserve_quality,
                                bool preserve_id,
                                dsg::input::fastq::FastqFileReader *fastqFileReader1,
                                dsg::input::fastq::FastqFileReader *fastqFileReader2,
                                std::string quality_mode,
                                double quality_ratio);

} // namespace spring

#endif // SPRING_REORDER_COMPRESS_QUALITY_ID_H_
