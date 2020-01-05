/*
* Copyright 2018 University of Illinois Board of Trustees and Stanford
University. All Rights Reserved.
* Licensed under the “Non-exclusive Research Use License for SPRING Software”
license (the "License");
* You may not use this file except in compliance with the License.
* The License is included in the distribution as license.pdf file.

* Software distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
limitations under the License.
*/

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
  bool help_flag = false, gzipped_input_flag = false, gzipped_output_flag = false,
        fasta_flag = false;
  std::vector<std::string> infile_vec, outfile_vec;
  std::string working_dir;
  int num_thr;
  po::options_description desc("Allowed options");
  desc.add_options()("help,h", po::bool_switch(&help_flag),
                     "produce help message")(
      "input-file,i",
      po::value<std::vector<std::string> >(&infile_vec)->multitoken(),
      "input FASTQ file name (two files for paired end)")(
      "output-file,o",
      po::value<std::vector<std::string> >(&outfile_vec)->multitoken(),
      "output FASTQ file name (for paired end, if only one file is "
      "specified, two output files will be created by suffixing .1 and .2.)")(
      "working-dir,w", po::value<std::string>(&working_dir)->default_value("."),
      "directory to create temporary files (default current directory)")(
      "num-threads,t", po::value<int>(&num_thr)->default_value(8),
      "number of threads (default 8)")(
      "gzipped-input", po::bool_switch(&gzipped_input_flag),
      "enable if compression input is gzipped fastq")(
      "gzipped-output", po::bool_switch(&gzipped_output_flag),
      "enable to output gzipped fastq ")(
      "fasta-input", po::bool_switch(&fasta_flag),
      "enable if input is fasta file (i.e., no qualities)");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  if (help_flag) {
    std::cout << desc << "\n";
    return 0;
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

  try {
    spring::spring_reorder(temp_dir, infile_vec, outfile_vec, num_thr,
                         gzipped_input_flag, gzipped_output_flag, fasta_flag);
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
