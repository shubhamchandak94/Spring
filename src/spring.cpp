#include <boost/filesystem.hpp>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
//#include "encoder.h"
//#include "pe_encode.h"
#include "preprocess.h"
//#include "reorder.h"
//#include "reorder_compress_quality_id.h"
#include "spring.h"

namespace spring {

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
    	("input-file,i", po::value< vector<string> >(&infile_vec), "input file name (specify two files for paired end)")
    	("output-file,o", po::value< vector<string> >(&outfile_vec), "output file name (for paired end decompression, if only one file is specified, two output files will be created by suffixing .1 and .2.)")
    	("num-threads,t", po::value<int>(&num_thr), "number of threads (default 8)")->default_value(8)
    	("allow_read_reordering,r", po::bool_switch(&pairing_only_flag), "do not retain read order during compression (paired reads still remain paired). For single end files, this leads to arbitrary reordering of the reads.")
    	("no-quality", po::bool_switch(&no_quality_flag), "do not retain quality values during compression")
    	("no-ids", po::bool_switch(&no_ids_flag), "do not retain read identifiers during compression")
    	("working-dir,w", po::value<std::string>(&working_dir), "directory to create temporary files (default pwd)")->default_value("")
	("ill-bin",po::bool_switch(&ill_bin_flag), "apply Illumina binning to quality scores before compression")
    	("quality-compressor,q", po::value<std::string>(&quality_compressor), "compressor to use for quality values: bcm or qvz (default qvz). For long reads, only bcm supported.")->default("qvz")
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

	if((!compress_flag && !decompress_flag) || compress_flag && decompress_flag) {
		std::cout << "Exactly one of compress or decompress needs to be specified \n";
		std::cout <<desc << "\n";
		return 1;
	}
	// generate randomly named temporary directory in the working directory
	while(true) {
		std::string random_str = "tmp."+random_string(10);
		std::string temp_dir = working_dir+"/"+random_str + '/';
		if(!boost::filesystem::exists(temp_dir))
			break;
	}
	if(!boost::filesystem::create_directory(temp_dir)) {
		throw std::runtime_error("Cannot create temporary directory.");
	}
	std::cout << "Temporary directory: " << temp_dir << "\n";
	if(compress_flag) 
		compress(temp_dir, infile_vec, outfile_vec, num_thr, pairing_only_flag, no_quality_flag, no_ids_flag, ill_bin_flag, quality_compressor, long_flag);
	else 
		decompress();

	if(compress_flag && long_flag) {
		std::cout << "Long flag detected.\n";	
		std::cout << "For long mode: allow_read_reordering flag is disabled and quality compressor is fixed to bcm.\n";
		quality_compressor = "bcm";
		pairing_only_flag = false;
	}
	// Error handling
	catch(std::runtime_error& e) {
		std::cout << "Program terminated unexpectedly with error: "<< e.what() << "\n";
		std::cout << "Deleting temporary directory...\n";
		boost::filesystem::remove_all(temp_dir);
		std::cout << desc << "\n";
		return 1;
	}	
	return 0;	
}

void compress(std::string &temp_dir, std::vector<std::string>& infile_vec, std::vector<std::string>& outfile_vec, int &num_thr, bool &pairing_only_flag, bool &no_quality_flag, bool &no_ids_flag, bool &ill_bin_flag, std::string &quality_compressor, bool &long_flag) {

	std::string infile_1, infile_2, outfile;
	bool paired_end, preserve_quality, preserve_id, preserve_order;
	// Check options
	preserve_order = !paired_end_only_flag;
	preserve_id = !no_ids_flag;
	preserve_quality = !no_quality_flag;
	switch(infile_vec.size()) {
		case 0: throw std::runtime_error("No input file specified";
			break;
		case 1: paired_end = false;
			infile_1 = infile_vec[0];
			break;
		case 2: paired_end = true;
			infile_1 = infile_vec[0];
			infile_2 = infile_vec[1];
			break;
		default: throw std::runtime_error("Too many input files specified");
	}
	if(outfile_vec.size() == 1)
		outfile = outfile_vec[0];
	else
		throw std::runtime_error("Number of output files not equal to 1");
	if(quality_compressor != "bcm" && quality_compressor != "qvz")
		throw std::runtime_error("Invalid quality compressor");

	preprocess(infile_1, infile_2, temp_dir, paired_end, preserve_id, preserve_quality, preserve_order, ill_bin_flag, quality_compressor, long_flag);
	/*		
	if (status != 0) throw std::runtime_error("Bad input file");
	std::ifstream f_meta(temp_dir + "/read_meta.txt");
	std::string max_readlen_str;
	std::getline(f_meta, max_readlen_str);
	int max_readlen = std::stoi(max_readlen_str);

	call_reorder(temp_dir, max_readlen, num_thr);
	call_encoder(temp_dir, max_readlen, num_thr);
	if (paired_end == true) pe_encode_main(temp_dir, false);
	fastqFileReader1->seekFromSet(0);
	if (paired_end == true)
	     fastqFileReader2->seekFromSet(0);
	reorder_compress_quality_id(temp_dir, max_readlen, num_thr,
	paired_end, false, true, true, fastqFileReader1, fastqFileReader2, "bsc",
	8.0);
	*/
	return;
}

void decompress() {
	throw std::runtime_error("Not implemented");
}

void call_reorder(const std::string &temp_dir, int max_readlen, int num_thr) {
  size_t bitset_size_reorder = (2 * max_readlen - 1) / 64 * 64 + 64;
  switch (bitset_size_reorder) {
    case 64:
      reorder_main<64>(temp_dir, max_readlen, num_thr);
      break;
    case 128:
      reorder_main<128>(temp_dir, max_readlen, num_thr);
      break;
    case 192:
      reorder_main<192>(temp_dir, max_readlen, num_thr);
      break;
    case 256:
      reorder_main<256>(temp_dir, max_readlen, num_thr);
      break;
    case 320:
      reorder_main<320>(temp_dir, max_readlen, num_thr);
      break;
    case 384:
      reorder_main<384>(temp_dir, max_readlen, num_thr);
      break;
    case 448:
      reorder_main<448>(temp_dir, max_readlen, num_thr);
      break;
    case 512:
      reorder_main<512>(temp_dir, max_readlen, num_thr);
      break;
    case 576:
      reorder_main<576>(temp_dir, max_readlen, num_thr);
      break;
    case 640:
      reorder_main<640>(temp_dir, max_readlen, num_thr);
      break;
    case 704:
      reorder_main<704>(temp_dir, max_readlen, num_thr);
      break;
    case 768:
      reorder_main<768>(temp_dir, max_readlen, num_thr);
      break;
    case 832:
      reorder_main<832>(temp_dir, max_readlen, num_thr);
      break;
    case 896:
      reorder_main<896>(temp_dir, max_readlen, num_thr);
      break;
    case 960:
      reorder_main<960>(temp_dir, max_readlen, num_thr);
      break;
    case 1024:
      reorder_main<1024>(temp_dir, max_readlen, num_thr);
      break;
    default:
      throw std::runtime_error("Wrong bitset size.");
  }
}

void call_encoder(const std::string &temp_dir, int max_readlen, int num_thr) {
  size_t bitset_size_encoder = (3 * max_readlen - 1) / 64 * 64 + 64;
  switch (bitset_size_encoder) {
    case 64:
      encoder_main<64>(temp_dir, max_readlen, num_thr);
      break;
    case 128:
      encoder_main<128>(temp_dir, max_readlen, num_thr);
      break;
    case 192:
      encoder_main<192>(temp_dir, max_readlen, num_thr);
      break;
    case 256:
      encoder_main<256>(temp_dir, max_readlen, num_thr);
      break;
    case 320:
      encoder_main<320>(temp_dir, max_readlen, num_thr);
      break;
    case 384:
      encoder_main<384>(temp_dir, max_readlen, num_thr);
      break;
    case 448:
      encoder_main<448>(temp_dir, max_readlen, num_thr);
      break;
    case 512:
      encoder_main<512>(temp_dir, max_readlen, num_thr);
      break;
    case 576:
      encoder_main<576>(temp_dir, max_readlen, num_thr);
      break;
    case 640:
      encoder_main<640>(temp_dir, max_readlen, num_thr);
      break;
    case 704:
      encoder_main<704>(temp_dir, max_readlen, num_thr);
      break;
    case 768:
      encoder_main<768>(temp_dir, max_readlen, num_thr);
      break;
    case 832:
      encoder_main<832>(temp_dir, max_readlen, num_thr);
      break;
    case 896:
      encoder_main<896>(temp_dir, max_readlen, num_thr);
      break;
    case 960:
      encoder_main<960>(temp_dir, max_readlen, num_thr);
      break;
    case 1024:
      encoder_main<1024>(temp_dir, max_readlen, num_thr);
      break;
    case 1088:
      encoder_main<1088>(temp_dir, max_readlen, num_thr);
      break;
    case 1152:
      encoder_main<1152>(temp_dir, max_readlen, num_thr);
      break;
    case 1216:
      encoder_main<1216>(temp_dir, max_readlen, num_thr);
      break;
    case 1280:
      encoder_main<1280>(temp_dir, max_readlen, num_thr);
      break;
    case 1344:
      encoder_main<1344>(temp_dir, max_readlen, num_thr);
      break;
    case 1408:
      encoder_main<1408>(temp_dir, max_readlen, num_thr);
      break;
    case 1472:
      encoder_main<1472>(temp_dir, max_readlen, num_thr);
      break;
    case 1536:
      encoder_main<1536>(temp_dir, max_readlen, num_thr);
      break;
    default:
      throw std::runtime_error("Wrong bitset size.");
  }
}

std::string random_string( size_t length )
{
    auto randchar = []() -> char
    {
        const char charset[] =
        "0123456789"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz";
        const size_t max_index = (sizeof(charset) - 1);
        return charset[ rand() % max_index ];
    };
    std::string str(length,0);
    std::generate_n( str.begin(), length, randchar );
    return str;
}

}  // namespace spring
