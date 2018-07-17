#ifndef SPRING_SPRING_H_
#define SPRING_SPRING_H_

#include <string>
#include "generation.h"
#include "input/fastq/FastqFileReader.h"

namespace spring {

void call_reorder(const std::string &working_dir, int max_readlen, int num_thr);
void call_encoder(const std::string &working_dir, int max_readlen, int num_thr);

void generate_streams_SPRING(
    dsg::input::fastq::FastqFileReader *fastqFileReader1,
    dsg::input::fastq::FastqFileReader *fastqFileReader2, int num_thr,
    bool paired_end);

}  // namespace spring

#endif  // SPRING_SPRING_H_
