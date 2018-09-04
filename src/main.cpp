#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <csignal>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "spring.h"

std::string temp_dir_global;  // for interrupt handling
bool temp_dir_flag_global = false;

void signalHandler(int signum) {
  std::cout << "Interrupt signal (" << signum << ") received.\n";
  std::cout << "Program terminated unexpectedly\n";
  if (temp_dir_flag_global) {
    std::cout << "Deleting temporary directory:" << temp_dir_global << "\n";
    boost::filesystem::remove_all(temp_dir_global);
  }
  exit(signum);
}

int main(int argc, char** argv) {
  // register signal SIGINT and signal handler
  signal(SIGINT, signalHandler);
  namespace po = boost::program_options;
  bool help_flag = false, compress_flag = false, decompress_flag = false,
       pairing_only_flag = false, no_quality_flag = false, no_ids_flag = false,
       long_flag = false;
  std::vector<std::string> infile_vec, outfile_vec, quality_opts;
  std::vector<uint64_t> decompress_range_vec;
  std::string working_dir;
  int num_thr;
  po::options_description desc("Allowed options");
  desc.add_options()("help,h", po::bool_switch(&help_flag),
                     "produce help message")(
      "compress,c", po::bool_switch(&compress_flag), "compress")(
      "decompress,d", po::bool_switch(&decompress_flag), "decompress")(
      "decompress_range", po::value<std::vector<uint64_t> >(&decompress_range_vec)->multitoken(),
      "--decompress_range start end\n(optional) decompress only reads (or read pairs for PE datasets) from start to end (both inclusive) (1 <= start <= end <= num_reads (or num_read_pairs for PE)). If -r was specified during compression, the range of reads does not correspond to the original order of reads in the FASTQ file.")(
      "input-file,i", po::value<std::vector<std::string> >(&infile_vec)->multitoken(),
      "input file name (two files for paired end)")(
      "output-file,o", po::value<std::vector<std::string> >(&outfile_vec)->multitoken(),
      "output file name (for paired end decompression, if only one file is "
      "specified, two output files will be created by suffixing .1 and .2.)")(
      "working-dir,w", po::value<std::string>(&working_dir)->default_value("."),
      "directory to create temporary files (default current directory)")(
      "num-threads,t", po::value<int>(&num_thr)->default_value(8),
      "number of threads (default 8)")(
      "allow_read_reordering,r", po::bool_switch(&pairing_only_flag),
      "do not retain read order during compression (paired reads still remain "
      "paired)")("no-quality", po::bool_switch(&no_quality_flag),
                  "do not retain quality values during compression")(
      "no-ids", po::bool_switch(&no_ids_flag),
      "do not retain read identifiers during compression")(
      "quality_opts,q", po::value<std::vector<std::string> >(&quality_opts)->multitoken(),
      "quality mode: possible modes are\n1. -q lossless (default)\n2. -q qvz qv_ratio (QVZ lossy compression, parameter qv_ratio roughly corresponds to bits used per quality value)\n3. -q ill_bin (Illumina 8-level binning)\n4. -q binary thr high low (binary (2-level) thresholding, quality binned to high if >= thr and to low if < thr)")(
      "long,l", po::bool_switch(&long_flag),
      "Use for compression of arbitrarily long read lengths. Can also provide "
      "better compression for reads with significant number of indels. "
      "-r disabled in this mode. For Illumina short "
      "reads, compression is better without -l flag.");
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  if (help_flag) {
    std::cout << desc << "\n";
    return 0;
  }

  if ((!compress_flag && !decompress_flag) ||
      (compress_flag && decompress_flag)) {
    std::cout
        << "Exactly one of compress or decompress needs to be specified \n";
    std::cout << desc << "\n";
    return 1;
  }
  // generate randomly named temporary directory in the working directory
  std::string temp_dir;
  while (true) {
    std::string random_str = "tmp." + spring::random_string(10);
    temp_dir = working_dir + "/" + random_str + '/';
    if (!boost::filesystem::exists(temp_dir)) break;
  }
  if (!boost::filesystem::create_directory(temp_dir)) {
    throw std::runtime_error("Cannot create temporary directory.");
  }
  std::cout << "Temporary directory: " << temp_dir << "\n";

  temp_dir_global = temp_dir;
  temp_dir_flag_global = true;

  if (compress_flag && long_flag) {
    std::cout << "Long flag detected.\n";
    if (pairing_only_flag) {
      std::cout << "For long mode: allow_read_reordering flag is disabled.\n";
      pairing_only_flag = false;
    }
  }
  try {
    if (compress_flag)
      spring::compress(temp_dir, infile_vec, outfile_vec, num_thr,
                       pairing_only_flag, no_quality_flag, no_ids_flag,
                       quality_opts, long_flag);
    else
      spring::decompress(temp_dir, infile_vec, outfile_vec, num_thr, decompress_range_vec);

  }
  // Error handling
  catch (std::runtime_error& e) {
    std::cout << "Program terminated unexpectedly with error: " << e.what()
              << "\n";
    std::cout << "Deleting temporary directory...\n";
    boost::filesystem::remove_all(temp_dir);
    temp_dir_flag_global = false;
    std::cout << desc << "\n";
    return 1;
  } catch (...) {
    std::cout << "Program terminated unexpectedly\n";
    std::cout << "Deleting temporary directory...\n";
    boost::filesystem::remove_all(temp_dir);
    temp_dir_flag_global = false;
    std::cout << desc << "\n";
    return 1;
  }
  boost::filesystem::remove_all(temp_dir);
  temp_dir_flag_global = false;
  return 0;
}
