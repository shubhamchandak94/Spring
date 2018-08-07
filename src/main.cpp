#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "spring.h"

int main(int argc, char** argv)
{
	namespace po = boost::program_options;
	bool help_flag, compress_flag, decompress_flag, pairing_only_flag, no_quality_flag, no_ids_flag, ill_bin_flag, long_flag;
	std::vector<std::string> infile_vec, outfile_vec;
	std::string working_dir, quality_compressor;
	int num_thr;
	po::options_description desc("Allowed options");
	desc.add_options()
    	("help,h",po::bool_switch(&help_flag), "produce help message")
    	("compress,c", po::bool_switch(&compress_flag), "compress")
    	("decompress,d", po::bool_switch(&decompress_flag), "decompress")
    	("input-file,i", po::value<std::vector<std::string> >(&infile_vec), "input file name (specify two files for paired end)")
    	("output-file,o", po::value<std::vector<std::string> >(&outfile_vec), "output file name (for paired end decompression, if only one file is specified, two output files will be created by suffixing .1 and .2.)")
    	("num-threads,t", po::value<int>(&num_thr)->default_value(8),
 "number of threads (default 8)")
    	("allow_read_reordering,r", po::bool_switch(&pairing_only_flag), "do not retain read order during compression (paired reads still remain paired). For single end files, this leads to arbitrary reordering of the reads.")
    	("no-quality", po::bool_switch(&no_quality_flag), "do not retain quality values during compression")
    	("no-ids", po::bool_switch(&no_ids_flag), "do not retain read identifiers during compression")
    	("working-dir,w", po::value<std::string>(&working_dir)->default_value(""), "directory to create temporary files (default pwd)")
	("ill-bin",po::bool_switch(&ill_bin_flag), "apply Illumina binning to quality scores before compression")
    	("quality-compressor,q", po::value<std::string>(&quality_compressor)->default_value("qvz"), "compressor to use for quality values: bcm or qvz (default qvz). For long reads, only bcm supported.")
	("long,l",po::bool_switch(&long_flag), "Use for compression of arbitrarily long read lengths. Can also provide better compression for reads with significant number of indels. Some other options might be disabled in this mode.")
	;
	po::parse_command_line(argc, argv, desc);
	if(help_flag) {
		std::cout <<desc << "\n";
		return 0;
	}
	if(quality_compressor == "Bcm" || quality_compressor == "BCM")
		quality_compressor = "bcm";
	if(quality_compressor == "Qvz" || quality_compressor == "QVZ")
		quality_compressor = "qvz";

	if((!compress_flag && !decompress_flag) || (compress_flag && decompress_flag)) {
		std::cout << "Exactly one of compress or decompress needs to be specified \n";
		std::cout <<desc << "\n";
		return 1;
	}
	// generate randomly named temporary directory in the working directory
	std::string temp_dir;
	while(true) {
		std::string random_str = "tmp."+spring::random_string(10);
		temp_dir = working_dir+"/"+random_str + '/';
		if(!boost::filesystem::exists(temp_dir))
			break;
	}
	if(!boost::filesystem::create_directory(temp_dir)) {
		throw std::runtime_error("Cannot create temporary directory.");
	}
	std::cout << "Temporary directory: " << temp_dir << "\n";
	
	if(compress_flag && long_flag) {
		std::cout << "Long flag detected.\n";	
		std::cout << "For long mode: allow_read_reordering flag is disabled and quality compressor is fixed to bcm.\n";
		quality_compressor = "bcm";
		pairing_only_flag = false;
	}
	try {
	if(compress_flag) 
		spring::compress(temp_dir, infile_vec, outfile_vec, num_thr, pairing_only_flag, no_quality_flag, no_ids_flag, ill_bin_flag, quality_compressor, long_flag);
	else 
		spring::decompress();

	}
	// Error handling
	catch(std::runtime_error& e) {
		std::cout << "Program terminated unexpectedly with error: "<< e.what() << "\n";
		std::cout << "Deleting temporary directory...\n";
		boost::filesystem::remove_all(temp_dir);
		std::cout << desc << "\n";
		return 1;
	}
	catch(...) {	
		std::cout << "Program terminated unexpectedly\n";
		std::cout << "Deleting temporary directory...\n";
		boost::filesystem::remove_all(temp_dir);
		std::cout << desc << "\n";
		return 1;
	}
	return 0;	
}
