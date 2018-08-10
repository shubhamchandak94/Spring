#ifndef SPRING_UTIL_H_
#define SPRING_UTIL_H_

#include <fstream>
#include <string>

namespace spring {

struct compression_params {
//	std::string quality_compressor;
	bool paired_end;
	bool preserve_order;
	bool preserve_quality;
	bool preserve_id;
	bool long_flag;
	bool ill_bin_flag;
	uint32_t num_reads;
	uint32_t num_reads_clean[2];
	uint32_t max_readlen;
	uint8_t paired_id_code;
	bool paired_id_match;
	int num_reads_per_chunk;
	int num_reads_per_chunk_long;
	int bcm_block_size;
	int num_thr;
};

int read_fastq_block(std::ifstream &fin, std::string *id_array, std::string *read_array, std::string *quality_array, int num_reads);

void write_fastq_block(std::ofstream &fout, std::string *id_array, std::string *read_array, std::string *quality_array, int num_reads);

void compress_id_block(const char* outfile_name, std::string *id_array, int num_ids);

void decompress_id_block(const char* infile_name, std::string *id_array, int num_ids);

void compress_quality_block_qvz(const char* outfile_name, std::string *quality_array, int num_lines, uint32_t *read_lengths);

void decompress_quality_block_qvz(const char* infile_name, std::string *quality_array, int num_lines, uint16_t *read_lengths);

void quantize_quality(std::string *quality_array, int num_lines, char *quantization_table);

void generate_illumina_binning_table(char *illumina_binning_table);

uint8_t find_id_pattern(const std::string &id_1, const std::string &id_2);

bool check_id_pattern(const std::string &id_1, const std::string &id_2,
                      uint8_t paired_id_code);

void write_dna_in_bits(std::string &read, std::ofstream &fout);

void read_dna_from_bits(std::string &read, std::ifstream &fin);
void reverse_complement(char *s, char *s1, int readlen,
                        char chartorevchar[128]);

std::string reverse_complement(std::string s, int readlen,
                               char chartorevchar[128]);

}  // namespace spring

#endif  // SPRING_UTIL_H_
