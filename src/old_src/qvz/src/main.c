
#include "util.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <fstream>
#include <iostream>
#include <string>

#include "codebook.h"
#include "qv_compressor.h"
#include "cluster.h"

#define ALPHABET_SIZE 72

/**
 *
 */
void encode(FILE *fout, struct qv_options_t *opts, uint32_t max_readlen, uint32_t numreads, char *quality_array, uint8_t* read_lengths, std::string &infile_order ,uint64_t startpos) {
	struct quality_file_t qv_info;
	struct distortion_t *dist;
	struct alphabet_t *alphabet = alloc_alphabet(ALPHABET_SIZE);
	uint32_t status;
	struct hrtimer_t stats, encoding, total;
//	FILE *fout, *funcompressed = NULL;
	FILE *funcompressed = NULL;
	uint64_t bytes_used;
	double distortion;

	start_timer(&total);

	if (opts->distortion == DISTORTION_CUSTOM) {
		dist = gen_custom_distortion(ALPHABET_SIZE, opts->dist_file);
	}
	else {
		dist = generate_distortion_matrix(ALPHABET_SIZE, opts->distortion);
	}
    
	qv_info.alphabet = alphabet;
	qv_info.dist = dist;
	qv_info.cluster_count = opts->clusters;
	qv_info.columns = max_readlen;
	qv_info.lines = numreads;
	//from alloc_lines & alloc_blocks - we'll allocate single block
	qv_info.block_count = 1;
	qv_info.blocks = (struct line_block_t *) calloc(qv_info.block_count, sizeof(struct line_block_t));
	qv_info.blocks[0].count = qv_info.lines;
	qv_info.blocks[0].quality_array = quality_array;
	qv_info.blocks[0].read_lengths = read_lengths;
	qv_info.blocks[0].infile_order = new std::string;
	*(qv_info.blocks[0].infile_order) = infile_order;
	qv_info.blocks[0].startpos = startpos;
//	qv_info.blocks[0].lines = (struct line_t *) calloc(qv_info.blocks[0].count, sizeof(struct line_t));
//	uint32_t order;
//	for(uint32_t i = 0; i < qv_info.lines; i++)
//	{
//		f_order->read((char*)&order, sizeof(uint32_t));
//		qv_info.blocks[0].lines[i].m_data = quality_array + (uint64_t)(order)*(readlen+1);
//	}
		
	// Set up clustering data structures
	qv_info.clusters = alloc_cluster_list(&qv_info);
	qv_info.opts = opts;

    
	// Then find stats and generate codebooks for each cluster
	start_timer(&stats);
	calculate_statistics(&qv_info);
	generate_codebooks(&qv_info);
	stop_timer(&stats);
    
	if (opts->verbose) {
		printf("Stats and codebook generation took %.4f seconds\n", get_timer_interval(&stats));
		// @todo expected distortion is inaccurate due to lack of pmf
		//printf("Expected distortion: %f\n", opts->e_dist);
	}
    
/*	std::cout << output_name << "\n";
	// Note that we want \r\n translation in the input
	// but we do not want it in the output
	fout = fopen(output_name, "wb");
	if (!fout) {
		perror("Unable to open output file");
		exit(1);
	}
*/    
    if (opts->uncompressed) {
        funcompressed = fopen(opts->uncompressed_name, "w");
        if (!funcompressed) {
            perror("Unable to open uncompressed file");
            exit(1);
        }
    }
	
	// @todo qv_compression should use quality_file structure with data in memory, now
	start_timer(&encoding);
	write_codebooks(fout, &qv_info);
    bytes_used = start_qv_compression(&qv_info, fout, &distortion, funcompressed);
	stop_timer(&encoding);
	stop_timer(&total);

	fclose(fout);
	delete qv_info.blocks[0].infile_order;
   	free_blocks(&qv_info); 
	free_cluster_list(qv_info.clusters);
	// Verbose stats
	if (opts->verbose) {
		// @todo add cluster info here
        switch (opts->distortion) {
            case DISTORTION_MANHATTAN:
                printf("L1 distortion: %f\n", distortion);
                break;
            case DISTORTION_MSE:
                printf("MSE distortion: %f\n", distortion);
                break;
            case DISTORTION_LORENTZ:
                printf("log(1+L1) distortion: %f\n", distortion);
                break;
			case DISTORTION_CUSTOM:
				printf("Custom distortion: %f\n", distortion);
				break;
            default:
                break;
        }
		printf("Lines: %llu\n", qv_info.lines);
		printf("Columns: %u\n", qv_info.columns);
		printf("Total bytes used: %llu\n", bytes_used);
		printf("Encoding took %.4f seconds.\n", get_timer_interval(&total));
		printf("Total time elapsed: %.4f seconds.\n", get_timer_interval(&total));
	}

	// Parse-able stats
	if (opts->stats) {
		printf("rate, %.4f, distortion, %.4f, time, %.4f, size, %llu \n", (bytes_used*8.)/((double)(qv_info.lines)*qv_info.columns), distortion, get_timer_interval(&total), bytes_used);
	}
}

/**
 *
 */
void decode(char *input_file, char *output_file, struct qv_options_t *opts, uint8_t *read_lengths) {
	FILE *fin, *fout;
	struct hrtimer_t timer;
	struct quality_file_t qv_info;
	struct alphabet_t *A = alloc_alphabet(ALPHABET_SIZE);
    
	qv_info.alphabet = A;
	qv_info.opts = opts;

	start_timer(&timer);

	fin = fopen(input_file, "rb");
	fout = fopen(output_file, "wt");
	if (!fin || !fout) {
		perror("Unable to open input or output files");
		exit(1);
	}

	read_codebooks(fin, &qv_info);
    start_qv_decompression(fout, fin, &qv_info, read_lengths);

	fclose(fout);
	fclose(fin);
	stop_timer(&timer);

	if (opts->verbose) {
		printf("Decoded %llu lines in %f seconds.\n", qv_info.lines, get_timer_interval(&timer));
	}
}

/**
 * Displays a usage name
 * @param name Program name string
 */
void usage(char *name) {
	printf("Usage: %s (options) [input file] [output file]\n", name);
	printf("Options are:\n");
	printf("   -q           : Store quality values in compressed file (default)\n");
	printf("   -x           : Extract quality values from compressed file\n");
	printf("   -f [ratio]   : Compress using [ratio] bits per bit of input entropy per symbol\n");
	printf("   -r [rate]    : Compress using fixed [rate] bits per symbol\n");
    printf("   -d [M|L|A]   : Optimize for MSE, Log(1+L1), L1 distortions, respectively (default: MSE)\n");
	printf("   -D [FILE]    : Optimize using the custom distortion matrix specified in FILE\n");
	printf("   -c [#]       : Compress using [#] clusters (default: 1)\n");
	printf("   -T [#]       : Use [#] as a threshold for cluster center movement (L2 norm) to declare a stable solution (default: 4).\n");
    printf("   -u [FILE]    : Write the uncompressed lossy values to FILE (default: off)\n");
	printf("   -h           : Print this help\n");
	printf("   -s           : Print summary stats\n");
	printf("   -v           : Enable verbose output\n");
	printf("\nFor custom distortion matrices, a 72x72 matrix of values must be provided as the cost of reconstructing\n");
	printf("the x-th row as the y-th column, where x and y range from 0 to 71 (inclusive) corresponding to the possible\n");
	printf("Phred scores.\n");
}

/**
 *
 */
/*
int main(int argc, char **argv) {
    char *input_name = 0;
	char *output_name = 0;
	struct qv_options_t opts;
	uint32_t i;

	uint8_t extract = 0;
	uint8_t file_idx = 0;

	opts.verbose = 0;
	opts.stats = 0;
	opts.ratio = 0.5;
	opts.clusters = 1;
    opts.uncompressed = 0;
    opts.distortion = DISTORTION_MSE;
	opts.cluster_threshold = 4;

	// No dependency, cross-platform command line parsing means no getopt
	// So we need to settle for less than optimal flexibility (no combining short opts, maybe that will be added later)
	i = 1;
	while (i < argc) {
		// Handle file names and reject any other untagged arguments
		if (argv[i][0] != '-') {
			switch (file_idx) {
				case 0:
					input_name = argv[i];
					file_idx = 1;
					break;
				case 1:
					output_name = argv[i];
					file_idx = 2;
					break;
				default:
					printf("Garbage argument \"%s\" detected.\n", argv[i]);
					usage(argv[0]);
					exit(1);
			}
			i += 1;
			continue;
		}

		// Flags for options
		switch(argv[i][1]) {
			case 'x':
				extract = 1;
				i += 1;
				break;
			case 'q':
				extract = 0;
				i += 1;
				break;
			case 'f':
				extract = 0;
				opts.ratio = atof(argv[i+1]);
				opts.mode = MODE_RATIO;
				i += 2;
				break;
			case 'r':
				extract = 0;
				opts.ratio = atof(argv[i+1]);
				opts.mode = MODE_FIXED;
				i += 2;
				break;
			case 'c':
				opts.clusters = atoi(argv[i+1]);
				i += 2;
				break;
			case 'v':
				opts.verbose = 1;
				i += 1;
				break;
			case 'h':
				usage(argv[0]);
				exit(0);
			case 's':
				opts.stats = 1;
				i += 1;
				break;
            case 'u':
                opts.uncompressed = 1;
                opts.uncompressed_name = argv[i+1];
                i += 2;
                break;
			case 'T':
				opts.cluster_threshold = atoi(argv[i+1]);
				i += 2;
				break;
            case 'd':
                switch (argv[i+1][0]) {
                    case 'M':
                        opts.distortion = DISTORTION_MSE;
                        break;
                    case 'L':
                        opts.distortion = DISTORTION_LORENTZ;
                        break;
                    case 'A':
                        opts.distortion = DISTORTION_MANHATTAN;
                        break;
                    default:
                        printf("Distortion measure not supported, using MSE.\n");
                        break;
                }
                i += 2;
                break;
			case 'D':
				opts.distortion = DISTORTION_CUSTOM;
				opts.dist_file = argv[i+1];
				i += 2;
				break;
			default:
				printf("Unrecognized option -%c.\n", argv[i][1]);
				usage(argv[0]);
				exit(1);
		}
	}

	if (file_idx != 2) {
		printf("Missing required filenames.\n");
		usage(argv[0]);
		exit(1);
	}

	if (opts.verbose) {
		if (extract) {
			printf("%s will be decoded to %s.\n", input_name, output_name);
		}
		else {
			printf("%s will be encoded as %s.\n", input_name, output_name);
			if (opts.mode == MODE_RATIO)
				printf("Ratio mode selected, targeting %f compression ratio.\n", opts.ratio);
			else if (opts.mode == MODE_FIXED)
				printf("Fixed-rate mode selected, targeting %f bits per symbol.\n", opts.ratio);
			else if (opts.mode == MODE_FIXED_MSE)
				printf("Fixed-MSE mode selected, targeting %f average distortion per context.\n", opts.ratio);

			switch (opts.distortion) {
				case DISTORTION_MSE:
					printf("MSE will be used as a distortion metric.\n");
					break;
				case DISTORTION_LORENTZ:
					printf("log(1+L1) will be used as a distortion metric.\n");
					break;
				case DISTORTION_MANHATTAN:
					printf("L1 will be used as a distortion metric.\n");
					break;
				case DISTORTION_CUSTOM:
					printf("A custom distortion metric stored in %s will be used.\n", opts.dist_file);
					break;
			}

			printf("Compression will use %d clusters, with a movement threshold of %.0f.\n", opts.clusters, opts.cluster_threshold);
		}
	}

	if (extract) {
		decode(input_name, output_name, &opts);
	}
	else {
		encode(input_name, output_name, &opts);
	}

#ifdef _WIN32
	system("pause");
#endif

	return 0;
}
*/
