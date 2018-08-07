#include <fstream>
#include <iostream>
#include <string>
#include "omp.h"
#include "params.h"
#include "preprocess.h"
#include "util.h"

namespace spring {

void preprocess(std::string &infile_1, std::string &infile_2,
               std::string &temp_dir, bool &paired_end, bool &preserve_id,
               bool &preserve_quality, bool &preserve_order, bool &ill_bin_flag, std::string &quality_compressor, bool &long_flag, compression_params &p)
{}
	


//
//  std::string outfileclean;
//  std::string outfileN;
//  std::string outfileorderN;
//  std::string outfileid[2];
//  std::string outfilenumreads;
//  std::string outfile_meta;
//  int max_readlen = -1;
//  std::string basedir = temp_dir;
//  outfileclean = basedir + "/input_clean.dna";
//  outfileN = basedir + "/input_N.dna";
//  outfileorderN = basedir + "/read_order_N.bin";
//  outfileid[0] = basedir + "/input_1.id";
//  outfileid[1] = basedir + "/input_2.id";
//  outfilenumreads = basedir + "/numreads.bin";
//  outfile_meta = basedir + "/read_meta.txt";
//
//  std::string id_1;
//  std::ofstream f_clean(outfileclean);
//  std::ofstream f_N(outfileN);
//  std::ofstream f_order_N(outfileorderN, std::ios::binary);
//
//  uint64_t total_reads[2] = {0, 0};
//  uint64_t readnum = 0, num_clean = 0;
//  uint8_t paired_id_code = 0;
//  bool paired_id_match = false;
//  int current_readlen;
//  std::vector<dsg::input::fastq::FastqRecord> fastqRecord;
//  // code 0: no pattern found
//  // code 1: */1 and */2 where * are same in both
//  // code 2: * and * where * are same in both
//  // code 3: * 1:# and * 2:# where * and # are common to both and * contains no
//  // space (used in new versions)
//  for (int j = 0; j < 2; j++) {
//    dsg::input::fastq::FastqFileReader *fastqFileReader = fastqFileReader1;
//    if (j == 1) fastqFileReader = fastqFileReader2;
//    if (j == 1 && paired_end == false) continue;
//
//    std::ofstream f_id;
//    std::ifstream fin_id_1;
//
//    if (preserve_id == true) {
//      f_id.open(outfileid[j]);
//      if (j == 1) {
//        fin_id_1.open(outfileid[0]);
//        // check first ids to detect patterns
//        size_t numRecordsRead = fastqFileReader->readRecords(1, &fastqRecord);
//        if (numRecordsRead == 1) {
//          std::string id_2 = fastqRecord[0].title;
//          std::getline(fin_id_1, id_1);
//          paired_id_code = find_id_pattern(id_1, id_2);
//          if (paired_id_code != 0) paired_id_match = true;
//        }
//        fastqFileReader->seekFromSet(0);
//        fin_id_1.close();
//        fin_id_1.open(outfileid[0]);
//      }
//    }
//    while (true) {
//      size_t numRecordsRead = fastqFileReader->readRecords(1, &fastqRecord);
//      if (numRecordsRead != 1) break;
//      if (preserve_id == true) {
//        f_id << fastqRecord[0].title << "\n";
//        if (j == 1 && paired_id_match) {
//          std::getline(fin_id_1, id_1);
//          if (fin_id_1.eof())
//            paired_id_match = false;
//          else
//            paired_id_match =
//                check_id_pattern(id_1, fastqRecord[0].title, paired_id_code);
//        }
//      }
//      current_readlen = (int)fastqRecord[0].sequence.length();
//      if (current_readlen >= 256) {
//        std::cout << "Read length cannot exceed 255. Read with length "
//                  << current_readlen << " found\n";
//        return -1;
//      }
//      if (current_readlen > max_readlen) max_readlen = current_readlen;
//      if (fastqRecord[0].sequence.find('N') != std::string::npos) {
//        f_N << fastqRecord[0].sequence << "\n";
//        f_order_N.write((char *)&readnum, sizeof(uint32_t));
//      } else {
//        num_clean++;
//        f_clean << fastqRecord[0].sequence << "\n";
//      }
//      if (preserve_quality == true) {
//        if ((int)fastqRecord[0].qualityScores.length() != current_readlen) {
//          std::cout << "Quality length does not match read length: "
//                    << current_readlen << " and "
//                    << fastqRecord[0].qualityScores.length() << " found.\n";
//          return -1;
//        }
//      }
//      readnum++;
//    }
//    total_reads[j] = readnum;
//  }
//  total_reads[1] = total_reads[1] - total_reads[0];
//  if (readnum > MAX_NUM_READS) {
//    std::cout << "Too many reads. HARC supports at most " << MAX_NUM_READS
//              << " reads\n";
//    return -1;
//  } else if (total_reads[1] != total_reads[0] && paired_end == true) {
//    std::cout << "Number of reads in the two paired files are not equal\n";
//    return -1;
//  } else {
//    std::ofstream f_numreads(outfilenumreads, std::ios::binary);
//    uint32_t num_clean_32 = (uint32_t)num_clean;
//    uint32_t readnum_32 = (uint32_t)readnum;
//    f_numreads.write((char *)&num_clean_32, sizeof(uint32_t));
//    f_numreads.write((char *)&readnum_32, sizeof(uint32_t));
//    if (paired_id_match == true)
//      f_numreads.write((char *)&paired_id_code, sizeof(uint8_t));
//    else {
//      paired_id_code = 0;
//      f_numreads.write((char *)&paired_id_code, sizeof(uint8_t));
//    }
//    std::cout << "Max Read length: " << max_readlen << "\n";
//    std::cout << "Total number of reads: " << readnum << "\n";
//    std::cout << "Total number of reads without N: " << num_clean << "\n";
//    if (preserve_id == true && paired_end == true)
//      std::cout << "Paired id match code: " << (int)paired_id_code << "\n";
//    f_numreads.close();
//    std::ofstream f_meta(outfile_meta);
//    f_meta << max_readlen << "\n";
//    f_meta.close();
//  }
//  std::cout << "Preprocessing Done!\n";
//
//  return 0;
//}

}  // namespace spring
