#include <fstream>
#include <iostream>
#include <string>
#include <cstdio>
#include <algorithm>
#include <omp.h>
#include "params.h"
#include "preprocess.h"
#include "util.h"
#include "bcm/bcm.h"

namespace spring {

void preprocess(std::string &infile_1, std::string &infile_2,
               std::string &temp_dir, bool &paired_end, bool &preserve_id,
               bool &preserve_quality, bool &preserve_order, bool &ill_bin_flag, std::string &quality_compressor, bool &long_flag, compression_params &cp)
{
  std::string outfileclean[2];
  std::string outfileN;
  std::string outfileorderN;
  std::string outfileid[2];
  std::string outfilequality[2];
  std::string outfileread[2];
  std::string outfilereadlength[2];
  std::string basedir = temp_dir;
  outfileclean[0] = basedir + "/input_clean_1.bin";
  outfileclean[1] = basedir + "/input_clean_2.bin";
  outfileN = basedir + "/input_N.dna";
  outfileorderN = basedir + "/read_order_N.bin";
  outfileid[0] = basedir + "/id_1";
  outfileid[1] = basedir + "/id_2";
  outfilequality[0] = basedir + "/quality_1";
  outfilequality[1] = basedir + "/quality_2";
  outfileread[0] = basedir + "/read_1";
  outfileread[1] = basedir + "/read_2";
  outfilereadlength[0] = basedir + "/readlength_1";
  outfilereadlength[1] = basedir + "/readlength_2";
  
  std::ifstream fin_1, fin_2;
  std::ofstream fout_clean_1, fout_clean_2;
  std::ofstream fout_N_1, fout_N_2;
  std::ofstream fout_order_N_1, fout_order_N_2;
  std::ofstream fout_id_1, fout_id_2;
  std::ofstream fout_quality_1, fout_quality_2;
  std::ofstream fout_readlength_1, fout_readlength_2;
  if(long_flag) {
    fin_1.open(infile_1);
    if(paired_end)
      fin_2.open(infile_2); 
  }
  else {
    fin_1.open(infile_1);
    fout_clean_1.open(outfileclean[0],std::ios::binary);
    fout_N_1.open(outfileN+".1");
    fout_order_N_1.open(outfileorderN+".1",std::ios::binary);
    if(!preserve_order) {
      if(preserve_id)
        fout_id_1.open(outfileid[0]);
      if(preserve_quality)
        fout_quality_1.open(outfilequality[0]);
	}
	if(paired_end) {
	  fin_2.open(infile_2);
	  fout_clean_2.open(outfileclean[1],std::ios::binary);
	  fout_N_2.open(outfileN+".2");
	  fout_order_N_2.open(outfileorderN+".2",std::ios::binary);
	  if(!preserve_order) {
		if(preserve_id)
		  fout_id_2.open(outfileid[1]);
		if(preserve_quality)
		  fout_quality_2.open(outfilequality[1]);
	  }
	}
  }
  
  uint32_t max_readlen = 0;
  uint64_t num_reads_1 = 0;
  uint64_t num_reads_2 = 0;
  uint64_t num_reads_clean_1 = 0;
  uint64_t num_reads_clean_2 = 0;
  uint32_t num_reads_per_chunk;
  if(long_flag)
	num_reads_per_chunk = NUM_READS_PER_CHUNK_LONG;
  else
	num_reads_per_chunk = NUM_READS_PER_CHUNK;
  uint8_t paired_id_code = 0;
  bool paired_id_match = false;
  
  char *illumina_binning_table = new char[128];
  if(ill_bin_flag && preserve_order)
	generate_illumina_binning_table(illumina_binning_table);
  
  // Check that we were able to open the input files and also look for 
  // paired end matching ids if relevant
  if(!fin_1.is_open())
	throw std::runtime_error("Error opening input file");
  if(paired_end) {
    if(!fin_2.is_open())
      throw std::runtime_error("Error opening input file");
    if(preserve_id) {
	  std::string id_1, id_2;
      std::getline(fin_1,id_1);
	  std::getline(fin_2,id_2);
	  paired_id_code = find_id_pattern(id_1,id_2);
	  if(paired_id_code != 0)
		paired_id_match = true;  
	  fin_1.seekg(0);
	  fin_2.seekg(0);
	}
  }	
  uint64_t num_reads_per_step = (uint64_t)cp.num_thr*num_reads_per_chunk;
  std::string *read_array = new std::string[num_reads_per_step];
  std::string *id_array_1 = new std::string[num_reads_per_step];
  std::string *id_array_2 = new std::string[num_reads_per_step];
  std::string *quality_array = new std::string[num_reads_per_step];
  bool *read_contains_N_array = new bool[num_reads_per_step];
  uint32_t *read_lengths_array = new uint32_t[num_reads_per_step];
  
  omp_set_num_threads(cp.num_thr);
  
  uint32_t num_chunks_done = 0;
  
  while(true) {
	bool done_1 = false, done_2 = false;
    uint32_t num_reads_read = read_fastq_block(fin_1, id_array_1, read_array, quality_array, num_reads_per_step);
	if(num_reads_read < num_reads_per_step)
	  done_1 = true;
	#pragma omp parallel
	{
		bool done = false;
		uint64_t tid = omp_get_thread_num();
		if(tid*num_reads_per_chunk >= num_reads_read) 
		  done = true;
		uint32_t num_reads_thr = std::min(num_reads_read, (tid+1)*num_reads_per_chunk) - tid*num_reads_per_chunk;
		if(!done) {
		  if(long_flag)	
	     	    fout_readlength_1.open(outfilereadlength[0]+"."+std::to_string(num_chunks_done+tid), std::ios::binary);
		  // check if reads and qualities have equal lengths 
		  for(uint32_t i = tid*num_reads_per_chunk; i < tid*num_reads_per_chunk + num_reads_thr; i++) {
			size_t len = read_array[i].size();
			if(read_array.size() == 0) 
			  throw std::runtime_error("Read of length 0 detected.");
			if(long_flag && len > MAX_READ_LEN_LONG) {
			  std::cout << "Max read length for long mode is " << MAX_READ_LEN_LONG << ", but found read of length " << len << "\n";
			  throw std::runtime_error("Too long read length");
			}
			if(!long_flag && len > MAX_READ_LEN) {
			  std::cout << "Max read length without long mode is " << MAX_READ_LEN << ", but found read of length " << len << "\n";
			  throw std::runtime_error("Too long read length (please try --long/-l flag).");
			}
			if(preserve_quality && (quality_array[i].size() != len))
			  throw std::runtime_error("Read length does not match quality length.");
			if(preserve_id && (id_array_1.size() == 0))
			  throw std::runtime_error("Identifier of length 0 detected.");
			read_lengths_array[i] = (uint32_t)len;
			
			// mark reads with N
			if(!long_flag) 
			  read_contains_N_array[i] = (read_array[i].find('N') != std::string::npos);  
			
			// Write read length to a file (for long mode)
			if(long_flag)
			  fout_readlength_1.write((char*)&read_lengths_array[i], sizeof(uint32_t));
		  }
		  if(long_flag)
		    fout_readlength_1.close();
		  // apply Illumina binning (if asked to do so)
		  if(preserve_quality && ill_bin_flag)
			quantize_quality(quality_array + tid*num_reads_per_chunk, num_reads_thr, illumina_binning_table);
		
		  if(!long_flag) {
			if(preserve_order) {
			  // Compress ids
			  if(preserve_id) {
				std::string outfile_name = outfileid[0] + "." + std::to_string(num_chunks_done+tid);
				compress_id_block(outfile_name.c_str(), id_array + tid*num_reads_per_chunk, num_reads_thr);  
			  }
			  // Compress qualities
			  if(preserve_quality) {
				std::string outfile_name = outfilequality[0] + "." + std::to_string(num_chunks_done+tid);
				if(quality_compressor == "qvz")
				  compress_quality_block_qvz(outfile_name.c_str(), quality_array + tid*num_reads_per_chunk, num_reads_thr, read_lengths_array + tid*num_reads_per_chunk);
				else
				  bcm::bcm_str_array_compress(outfile_name.c_str(), quality_array + tid*num_reads_per_chunk, num_reads_thr, read_lengths_array + tid*num_reads_per_chunk);
			  }	
			}
		  }
		  else {
			// Compress read lengths file
			std::string infile_name = outfilereadlength[0]+"."+std::to_string(num_chunks_done+tid);
			std::string outfile_name = outfilereadlength[0]+"."+std::to_string(num_chunks_done+tid) + ".bcm";
			bcm::bcm_compress(infile_name.c_str(), outfile_name.c_str());
			remove(infile_name.c_str());
			// Compress ids
			if(preserve_id) {
			  std::string outfile_name = outfileid[0] + "." + std::to_string(num_chunks_done+tid);
			  compress_id_block(outfile_name.c_str(), id_array + tid*num_reads_per_chunk, num_reads_thr);
			}
			// Compress qualities
			if(preserve_quality) {
			  std::string outfile_name = outfilequality[0] + "." + std::to_string(num_chunks_done+tid);
			  bcm::bcm_str_array_compress(outfile_name.c_str(), quality_array + tid*num_reads_per_chunk, num_reads_thr, read_lengths_array + tid*num_reads_per_chunk);
			}
			// Compress reads
			outfile_name = outfileread[0] + "." + std::to_string(num_chunks_done+tid);
			bcm::bcm_str_array_compress(outfile_name.c_str(), read_array + tid*num_reads_per_chunk, num_reads_thr, read_lengths_array + tid*num_reads_per_chunk);
		  }
		}
	} // omp parallel
	if(!long_flag) {
	  // write reads and read_order_N to respective files

	  if(!preserve_order) {
		// write qualities to file (after Illumina binning if needed)
		
		// write ids to file
	  }
	}
	
  }
  
  delete[] read_array;
  delete[] id_array_1;
  delete[] id_array_2;
  delete[] quality_array;
  delete[] read_contains_N_array;
  delete[] read_lengths_array;
  
  if(!long_flag && paired_end) {
	// merge input_N and input_order_N for the two files
    
  }
  
  if(paired_id_match) {
	// delete id files for second file since we found a pattern
	
  }	
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
