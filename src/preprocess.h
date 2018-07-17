#ifndef SPRING_PREPROCESS_H_
#define SPRING_PREPROCESS_H_

#include <string>
#include "input/fastq/FastqFileReader.h"

namespace spring {

uint8_t find_id_pattern(const std::string &id_1, const std::string &id_2);

bool check_id_pattern(const std::string &id_1, const std::string &id_2,
                      uint8_t paired_id_code);

int preprocess(dsg::input::fastq::FastqFileReader *fastqFileReader1,
               dsg::input::fastq::FastqFileReader *fastqFileReader2,
               const std::string &working_dir, bool paired_end, bool preserve_id,
               bool preserve_quality);

}  // namespace spring

#endif  // SPRING_PREPROCESS_H_
