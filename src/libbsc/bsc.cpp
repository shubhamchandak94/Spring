/*-----------------------------------------------------------*/
/* Block Sorting, Lossless Data Compression Library.         */
/* Block Sorting Compressor                                  */
/*-----------------------------------------------------------*/

/*--

This file is a part of bsc and/or libbsc, a program and a library for
lossless, block-sorting data compression.

   Copyright (c) 2009-2012 Ilya Grebnov <ilya.grebnov@gmail.com>

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

Please see the file LICENSE for full copyright information and file AUTHORS
for full list of contributors.

See also the bsc and libbsc web site:
  http://libbsc.com/ for more information.

--*/

/*--

Sort Transform is patented by Michael Schindler under US patent 6,199,064.
However for research purposes this algorithm is included in this software.
So if you are of the type who should worry about this (making money) worry away.
The author shall have no liability with respect to the infringement of
copyrights, trade secrets or any patents by this software. In no event will
the author be liable for any lost revenue or profits or other special,
indirect and consequential damages.

Sort Transform is disabled by default and can be enabled by defining the
preprocessor macro LIBBSC_SORT_TRANSFORM_SUPPORT at compile time.

--*/

#define _CRT_SECURE_NO_WARNINGS

#include <string>
#include <vector>
#include <stdexcept>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "libbsc/libbsc/libbsc.h"
#include "libbsc/libbsc/filters.h"
#include "libbsc/libbsc/platform/platform.h"
#include "libbsc/bsc.h"
#include "params.h"

namespace spring {
namespace bsc {

class bsc_class {

#pragma pack(push, 1)

#define LIBBSC_CONTEXTS_AUTODETECT   3

typedef struct BSC_BLOCK_HEADER
{
    long long       blockOffset;
    signed char     recordSize;
    signed char     sortingContexts;
} BSC_BLOCK_HEADER;

#pragma pack(pop)

int paramBlockSize                = 25 * 1024 * 1024;
int paramBlockSorter              = LIBBSC_BLOCKSORTER_BWT;
int paramCoder                    = LIBBSC_CODER_QLFC_STATIC;
int paramSortingContexts          = LIBBSC_CONTEXTS_FOLLOWING;

int paramEnableParallelProcessing = 0;
int paramEnableMultiThreading     = 1;
int paramEnableFastMode           = 1;
int paramEnableLargePages         = 0;
int paramEnableCUDA               = 0;
int paramEnableSegmentation       = 0;
int paramEnableReordering         = 0;
int paramEnableLZP                = 1;
int paramLZPHashSize              = 16;
int paramLZPMinLen                = 128;
int paramnumthr			  = 1;

int paramFeatures()
{
    int features =
        (paramEnableFastMode       ? LIBBSC_FEATURE_FASTMODE       : LIBBSC_FEATURE_NONE) |
        (paramEnableMultiThreading ? LIBBSC_FEATURE_MULTITHREADING : LIBBSC_FEATURE_NONE) |
        (paramEnableLargePages     ? LIBBSC_FEATURE_LARGEPAGES     : LIBBSC_FEATURE_NONE) |
        (paramEnableCUDA           ? LIBBSC_FEATURE_CUDA           : LIBBSC_FEATURE_NONE)
    ;

    return features;
}

#define BSC_FSEEK fseek
#define BSC_FTELL ftell
#define BSC_FILEOFFSET long


int segmentedBlock[256];

void Compression(char * argv[])
{
    if (!paramEnableLZP)
    {
        paramLZPHashSize = 0;
        paramLZPMinLen = 0;
    }

    FILE * fInput = fopen(argv[2], "rb");
    if (fInput == NULL)
    {
        fprintf(stderr, "Can't open input file: %s!\n", argv[2]);
	throw std::runtime_error("BSC error.");
    }

    FILE * fOutput = fopen(argv[3], "wb");
    if (fOutput == NULL)
    {
        fprintf(stderr, "Can't create output file: %s!\n", argv[3]);
	throw std::runtime_error("BSC error.");
    }

    if (BSC_FSEEK(fInput, 0, SEEK_END))
    {
        fprintf(stderr, "IO error on file: %s!\n", argv[2]);
	throw std::runtime_error("BSC error.");
    }

    BSC_FILEOFFSET fileSize = BSC_FTELL(fInput);
    if (fileSize < 0)
    {
        fprintf(stderr, "IO error on file: %s!\n", argv[2]);
	throw std::runtime_error("BSC error.");
    }

    if (BSC_FSEEK(fInput, 0, SEEK_SET))
    {
        fprintf(stderr, "IO error on file: %s!\n", argv[2]);
	throw std::runtime_error("BSC error.");
    }

    if (paramBlockSize > fileSize)
    {
        paramBlockSize = (int)fileSize;
    }

    int nBlocks = paramBlockSize > 0 ? (int)((fileSize + paramBlockSize - 1) / paramBlockSize) : 0;
    if (fwrite(&nBlocks, sizeof(nBlocks), 1, fOutput) != 1)
    {
        fprintf(stderr, "IO error on file: %s!\n", argv[3]);
	throw std::runtime_error("BSC error.");
    }

    int segmentationStart = 0, segmentationEnd = 0;

    {
        unsigned char * buffer = (unsigned char *)bsc_malloc(paramBlockSize + LIBBSC_HEADER_SIZE);
        if (buffer == NULL)
        {
            {

                fprintf(stderr, "Not enough memory! Please check README file for more information.\n");
		throw std::runtime_error("BSC error.");
            }
        }

        while (true)
        {
            BSC_FILEOFFSET  blockOffset     = 0;
            int             dataSize        = 0;

            {
                if ((feof(fInput) == 0) && (BSC_FTELL(fInput) != fileSize))
                {

                    blockOffset = BSC_FTELL(fInput);

                    int currentBlockSize = paramBlockSize;
                    if (paramEnableSegmentation)
                    {
                        if (segmentationEnd - segmentationStart > 1) currentBlockSize = segmentedBlock[segmentationStart];
                    }

                    dataSize = (int)fread(buffer, 1, currentBlockSize, fInput);
                    if (dataSize <= 0)
                    {
                        fprintf(stderr, "\nIO error on file: %s!\n", argv[2]);
			throw std::runtime_error("BSC error.");
                    }

                    if (paramEnableSegmentation)
                    {
                        bool bSegmentation = false;

                        if (segmentationStart == segmentationEnd) bSegmentation = true;
                        if ((segmentationEnd - segmentationStart == 1) && (dataSize != segmentedBlock[segmentationStart])) bSegmentation = true;

                        if (bSegmentation)
                        {
                            segmentationStart = 0; segmentationEnd = bsc_detect_segments(buffer, dataSize, segmentedBlock, 256, paramFeatures());
                            if (segmentationEnd <= LIBBSC_NO_ERROR)
                            {
                                switch (segmentationEnd)
                                {
                                    case LIBBSC_NOT_ENOUGH_MEMORY   : fprintf(stderr, "\nNot enough memory! Please check README file for more information.\n"); break;
                                    default                         : fprintf(stderr, "\nInternal program error, please contact the author!\n");
                                }
				throw std::runtime_error("BSC error.");
                            }
                        }

                        int newDataSize = segmentedBlock[segmentationStart++];
                        if (dataSize != newDataSize)
                        {
                            BSC_FILEOFFSET pos = BSC_FTELL(fInput) - dataSize + newDataSize;
                            BSC_FSEEK(fInput, pos, SEEK_SET);
                            dataSize = newDataSize;
                        }
                    }
                }
            }

            if (dataSize == 0) break;

            signed char recordSize = 1;
            if (paramEnableReordering)
            {
                recordSize = bsc_detect_recordsize(buffer, dataSize, paramFeatures());
                if (recordSize < LIBBSC_NO_ERROR)
                {
                    {
                        switch (recordSize)
                        {
                            case LIBBSC_NOT_ENOUGH_MEMORY   : fprintf(stderr, "\nNot enough memory! Please check README file for more information.\n"); break;
                            default                         : fprintf(stderr, "\nInternal program error, please contact the author!\n");
                        }
			throw std::runtime_error("BSC error.");
                    }
                }
                if (recordSize > 1)
                {
                    int result = bsc_reorder_forward(buffer, dataSize, recordSize, paramFeatures());
                    if (result != LIBBSC_NO_ERROR)
                    {
                        {
                            switch (result)
                            {
                                case LIBBSC_NOT_ENOUGH_MEMORY   : fprintf(stderr, "\nNot enough memory! Please check README file for more information.\n"); break;
                                default                         : fprintf(stderr, "\nInternal program error, please contact the author!\n");
                            }
			    throw std::runtime_error("BSC error.");
                        }
                    }
                }
            }

            signed char sortingContexts = paramSortingContexts;
            if (paramSortingContexts == LIBBSC_CONTEXTS_AUTODETECT)
            {
                sortingContexts = bsc_detect_contextsorder(buffer, dataSize, paramFeatures());
                if (sortingContexts < LIBBSC_NO_ERROR)
                {
                    {
                        switch (sortingContexts)
                        {
                            case LIBBSC_NOT_ENOUGH_MEMORY   : fprintf(stderr, "\nNot enough memory!\n"); break;
                            default                         : fprintf(stderr, "\nInternal program error, please contact the author!\n");
                        }
			throw std::runtime_error("BSC error.");
                    }
                }
            }
            if (sortingContexts == LIBBSC_CONTEXTS_PRECEDING)
            {
                int result = bsc_reverse_block(buffer, dataSize, paramFeatures());
                if (result != LIBBSC_NO_ERROR)
                {
                    {
                        fprintf(stderr, "\nInternal program error, please contact the author!\n");
			throw std::runtime_error("BSC error.");
                    }
                }
            }

            int blockSize = bsc_compress(buffer, buffer, dataSize, paramLZPHashSize, paramLZPMinLen, paramBlockSorter, paramCoder, paramFeatures());
            if (blockSize == LIBBSC_NOT_COMPRESSIBLE)
            {
                {
                    sortingContexts = LIBBSC_CONTEXTS_FOLLOWING; recordSize = 1;

                    BSC_FILEOFFSET pos = BSC_FTELL(fInput);
                    {
                        BSC_FSEEK(fInput, blockOffset, SEEK_SET);
                        if (dataSize != (int)fread(buffer, 1, dataSize, fInput))
                        {
                            fprintf(stderr, "\nInternal program error, please contact the author!\n");
			    throw std::runtime_error("BSC error.");
                        }
                    }
                    BSC_FSEEK(fInput, pos, SEEK_SET);
                }

                blockSize = bsc_store(buffer, buffer, dataSize, paramFeatures());
            }
            if (blockSize < LIBBSC_NO_ERROR)
            {
                {
                    switch (blockSize)
                    {
                        case LIBBSC_NOT_ENOUGH_MEMORY       : fprintf(stderr, "\nNot enough memory! Please check README file for more information.\n"); break;
                        case LIBBSC_NOT_SUPPORTED           : fprintf(stderr, "\nSpecified compression method is not supported on this platform!\n"); break;
                        case LIBBSC_GPU_ERROR               : fprintf(stderr, "\nGeneral GPU failure! Please check README file for more information.\n"); break;
                        case LIBBSC_GPU_NOT_SUPPORTED       : fprintf(stderr, "\nYour GPU is not supported! Please check README file for more information.\n"); break;
                        case LIBBSC_GPU_NOT_ENOUGH_MEMORY   : fprintf(stderr, "\nNot enough GPU memory! Please check README file for more information.\n"); break;

                        default                             : fprintf(stderr, "\nInternal program error, please contact the author!\n");
                    }
		    throw std::runtime_error("BSC error.");
                }
            }

            {
                BSC_BLOCK_HEADER header = {blockOffset, recordSize, sortingContexts};

                if (fwrite(&header, sizeof(BSC_BLOCK_HEADER), 1, fOutput) != 1)
                {
                    fprintf(stderr, "\nIO error on file: %s!\n", argv[3]);
		    throw std::runtime_error("BSC error.");
                }

                if ((int)fwrite(buffer, 1, blockSize, fOutput) != blockSize)
                {
                    fprintf(stderr, "\nIO error on file: %s!\n", argv[3]);
		    throw std::runtime_error("BSC error.");
                }
            }

        }

        bsc_free(buffer);
    }

    fclose(fInput); fclose(fOutput);
}

void Decompression(char * argv[])
{
    FILE * fInput = fopen(argv[2], "rb");
    if (fInput == NULL)
    {
        fprintf(stderr, "Can't open input file: %s!\n", argv[2]);
	throw std::runtime_error("BSC error.");
    }

    FILE * fOutput = fopen(argv[3], "wb");
    if (fOutput == NULL)
    {
        fprintf(stderr, "Can't create output file: %s!\n", argv[3]);
	throw std::runtime_error("BSC error.");
    }

    if (BSC_FSEEK(fInput, 0, SEEK_END))
    {
        fprintf(stderr, "IO error on file: %s!\n", argv[2]);
	throw std::runtime_error("BSC error.");
    }

    BSC_FILEOFFSET fileSize = BSC_FTELL(fInput);
    if (fileSize < 0)
    {
        fprintf(stderr, "IO error on file: %s!\n", argv[2]);
	throw std::runtime_error("BSC error.");
    }

    if (BSC_FSEEK(fInput, 0, SEEK_SET))
    {
        fprintf(stderr, "IO error on file: %s!\n", argv[2]);
	throw std::runtime_error("BSC error.");
    }

    int nBlocks = 0;
    if (fread(&nBlocks, sizeof(nBlocks), 1, fInput) != 1)
    {
        fprintf(stderr, "This is not bsc archive!\n");
	throw std::runtime_error("BSC error.");
    }

    {
        int bufferSize = -1; unsigned char * buffer = NULL;

        while (true)
        {
            BSC_FILEOFFSET  blockOffset     = 0;

            signed char     sortingContexts = 0;
            signed char     recordSize      = 0;
            int             blockSize       = 0;
            int             dataSize        = 0;

            {
                if ((feof(fInput) == 0) && (BSC_FTELL(fInput) != fileSize))
                {

                    BSC_BLOCK_HEADER header = {0, 0, 0};
                    if (fread(&header, sizeof(BSC_BLOCK_HEADER), 1, fInput) != 1)
                    {
                        fprintf(stderr, "\nUnexpected end of file: %s!\n", argv[2]);
			throw std::runtime_error("BSC error.");
                    }

                    recordSize = header.recordSize;
                    if (recordSize < 1)
                    {
                        fprintf(stderr, "\nThis is not bsc archive or invalid compression method!\n");
			throw std::runtime_error("BSC error.");
                    }

                    sortingContexts = header.sortingContexts;
                    if ((sortingContexts != LIBBSC_CONTEXTS_FOLLOWING) && (sortingContexts != LIBBSC_CONTEXTS_PRECEDING))
                    {
                        fprintf(stderr, "\nThis is not bsc archive or invalid compression method!\n");
			throw std::runtime_error("BSC error.");
                    }

                    blockOffset = (BSC_FILEOFFSET)header.blockOffset;

                    unsigned char bscBlockHeader[LIBBSC_HEADER_SIZE];

                    if (fread(bscBlockHeader, LIBBSC_HEADER_SIZE, 1, fInput) != 1)
                    {
                        fprintf(stderr, "\nUnexpected end of file: %s!\n", argv[2]);
			throw std::runtime_error("BSC error.");
                    }

                    if (bsc_block_info(bscBlockHeader, LIBBSC_HEADER_SIZE, &blockSize, &dataSize, paramFeatures()) != LIBBSC_NO_ERROR)
                    {
                        fprintf(stderr, "\nThis is not bsc archive or invalid compression method!\n");
			throw std::runtime_error("BSC error.");
                    }

                    if ((blockSize > bufferSize) || (dataSize > bufferSize))
                    {
                        if (blockSize > bufferSize) bufferSize = blockSize;
                        if (dataSize  > bufferSize) bufferSize = dataSize;

                        if (buffer != NULL) bsc_free(buffer); buffer = (unsigned char *)bsc_malloc(bufferSize);
                    }

                    if (buffer == NULL)
                    {
                        fprintf(stderr, "\nNot enough memory! Please check README file for more information.\n");
			throw std::runtime_error("BSC error.");
                    }

                    memcpy(buffer, bscBlockHeader, LIBBSC_HEADER_SIZE);

                    if (fread(buffer + LIBBSC_HEADER_SIZE, blockSize - LIBBSC_HEADER_SIZE, 1, fInput) != 1)
                    {
                        fprintf(stderr, "\nUnexpected end of file: %s!\n", argv[2]);
			throw std::runtime_error("BSC error.");
                    }
                }
            }

            if (dataSize == 0) break;

            int result = bsc_decompress(buffer, blockSize, buffer, dataSize, paramFeatures());
            if (result < LIBBSC_NO_ERROR)
            {
                {
                    switch (result)
                    {
                        case LIBBSC_DATA_CORRUPT            : fprintf(stderr, "\nThe compressed data is corrupted!\n"); break;
                        case LIBBSC_NOT_ENOUGH_MEMORY       : fprintf(stderr, "\nNot enough memory! Please check README file for more information.\n"); break;
                        case LIBBSC_GPU_ERROR               : fprintf(stderr, "\nGeneral GPU failure! Please check README file for more information.\n"); break;
                        case LIBBSC_GPU_NOT_SUPPORTED       : fprintf(stderr, "\nYour GPU is not supported! Please check README file for more information.\n"); break;
                        case LIBBSC_GPU_NOT_ENOUGH_MEMORY   : fprintf(stderr, "\nNot enough GPU memory! Please check README file for more information.\n"); break;

                        default                             : fprintf(stderr, "\nInternal program error, please contact the author!\n");
                    }
		    throw std::runtime_error("BSC error.");
                }
            }

            if (sortingContexts == LIBBSC_CONTEXTS_PRECEDING)
            {
                result = bsc_reverse_block(buffer, dataSize, paramFeatures());
                if (result != LIBBSC_NO_ERROR)
                {
                    {
                        fprintf(stderr, "\nInternal program error, please contact the author!\n");
			throw std::runtime_error("BSC error.");
                    }
                }
            }

            if (recordSize > 1)
            {
                result = bsc_reorder_reverse(buffer, dataSize, recordSize, paramFeatures());
                if (result != LIBBSC_NO_ERROR)
                {
                    {
                        switch (result)
                        {
                            case LIBBSC_NOT_ENOUGH_MEMORY   : fprintf(stderr, "\nNot enough memory! Please check README file for more information.\n"); break;
                            default                         : fprintf(stderr, "\nInternal program error, please contact the author!\n");
                        }
			throw std::runtime_error("BSC error.");
                    }
                }
            }

            {
                if (BSC_FSEEK(fOutput, blockOffset, SEEK_SET))
                {
                    fprintf(stderr, "\nIO error on file: %s!\n", argv[3]);
		    throw std::runtime_error("BSC error.");
                }

                if ((int)fwrite(buffer, 1, dataSize, fOutput) != dataSize)
                {
                    fprintf(stderr, "\nIO error on file: %s!\n", argv[3]);
		    throw std::runtime_error("BSC error.");
                }
            }
        }

        if (buffer != NULL) bsc_free(buffer);
    }

    if (BSC_FSEEK(fOutput, 0, SEEK_END))
    {
        fprintf(stderr, "IO error on file: %s!\n", argv[3]);
	throw std::runtime_error("BSC error.");
    }

    fclose(fInput); fclose(fOutput);
}

void ShowUsage(void)
{
    fprintf(stdout, "Usage: bsc <e|d> inputfile outputfile <options>\n\n");

    fprintf(stdout, "Block sorting options:\n");
    fprintf(stdout, "  -b<size> Block size in megabytes, default: -b25\n");
    fprintf(stdout, "             minimum: -b1, maximum: -b1024\n");
    fprintf(stdout, "  -m<algo> Block sorting algorithm, default: -m0\n");
    fprintf(stdout, "             -m0 Burrows Wheeler Transform (default)\n");
#ifdef LIBBSC_SORT_TRANSFORM_SUPPORT
    fprintf(stdout, "             -m3..8 Sort Transform of order n\n");
#endif
    fprintf(stdout, "  -c<ctx>  Contexts for sorting, default: -cf\n");
    fprintf(stdout, "             -cf Following contexts (default)\n");
    fprintf(stdout, "             -cp Preceding contexts\n");
    fprintf(stdout, "             -ca Autodetect (experimental)\n");
    fprintf(stdout, "  -e<algo> Entropy encoding algorithm, default: -e1\n");
    fprintf(stdout, "             -e1 Static Quantized Local Frequency Coding (default)\n");
    fprintf(stdout, "             -e2 Adaptive Quantized Local Frequency Coding (best compression)\n");
   
    fprintf(stdout, "\nPreprocessing options:\n");
    fprintf(stdout, "  -p       Disable all preprocessing techniques\n");
    fprintf(stdout, "  -s       Enable segmentation (adaptive block size), default: disable\n");
    fprintf(stdout, "  -r       Enable structured data reordering, default: disable\n");
    fprintf(stdout, "  -l       Enable Lempel-Ziv preprocessing, default: enable\n");
    fprintf(stdout, "  -H<size> LZP dictionary size in bits, default: -H16\n");
    fprintf(stdout, "             minimum: -H10, maximum: -H28\n");
    fprintf(stdout, "  -M<size> LZP minimum match length, default: -M128\n");
    fprintf(stdout, "             minimum: -M4, maximum: -M255\n");

    fprintf(stdout, "\nPlatform specific options:\n");

    fprintf(stdout,"\nOptions may be combined into one, like -b128p -m5e1\n");
    throw std::runtime_error("BSC error.");
}

void ProcessSwitch(char * s)
{
    if (*s == 0)
    {
        ShowUsage();
    }

    for (; *s != 0; )
    {
        switch (*s++)
        {
            case 'b':
            {
                char * strNum = s; while ((*s >= '0') && (*s <= '9')) s++;
                paramBlockSize = atoi(strNum) * 1024 * 1024;
                if ((paramBlockSize < 1024 * 1024) || (paramBlockSize > 1024 * 1024 * 1024)) ShowUsage();
                break;
            }

            case 'm':
            {
                char * strNum = s; while ((*s >= '0') && (*s <= '9')) s++;
                switch (atoi(strNum))
                {
                    case 0   : paramBlockSorter = LIBBSC_BLOCKSORTER_BWT; break;
                    default  : ShowUsage();
                }
                break;
            }

            case 'c':
            {
                switch (*s++)
                {
                    case 'f' : paramSortingContexts = LIBBSC_CONTEXTS_FOLLOWING;  break;
                    case 'p' : paramSortingContexts = LIBBSC_CONTEXTS_PRECEDING;  break;
                    case 'a' : paramSortingContexts = LIBBSC_CONTEXTS_AUTODETECT; break;
                    default  : ShowUsage();
                }
                break;
            }

            case 'e':
            {
                switch (*s++)
                {
                    case '1' : paramCoder = LIBBSC_CODER_QLFC_STATIC;   break;
                    case '2' : paramCoder = LIBBSC_CODER_QLFC_ADAPTIVE; break;
                    default  : ShowUsage();
                }
                break;
            }

            case 'H':
            {
                char * strNum = s; while ((*s >= '0') && (*s <= '9')) s++;
                paramLZPHashSize = atoi(strNum);
                if ((paramLZPHashSize < 10) || (paramLZPHashSize > 28)) ShowUsage();
                break;
            }

            case 'M':
            {
                char * strNum = s; while ((*s >= '0') && (*s <= '9')) s++;
                paramLZPMinLen = atoi(strNum);
                if ((paramLZPMinLen < 4) || (paramLZPMinLen > 255)) ShowUsage();
                break;
            }

            case 'l': paramEnableLZP            = 1; break;
            case 's': paramEnableSegmentation   = 1; break;
            case 'r': paramEnableReordering     = 1; break;

            case 'p': paramEnableLZP = paramEnableSegmentation = paramEnableReordering = 0; break;

            default : ShowUsage();
        }
    }
}

void ProcessCommandline(int argc, char * argv[])
{
    if (argc < 4 || strlen(argv[1]) != 1)
    {
        ShowUsage();
    }

    for (int i = 4; i < argc; ++i)
    {
        if (argv[i][0] == '-')
        {
            ProcessSwitch(&argv[i][1]);
        }
        else
        {
            ShowUsage();
        }
    }
}
public:
int bsc_main(int argc, char * argv[])
{

    ProcessCommandline(argc, argv);

    if (bsc_init(paramFeatures()) != LIBBSC_NO_ERROR)
    {
        fprintf(stderr, "\nInternal program error, please contact the author!\n");
	throw std::runtime_error("BSC error.");
    }

    switch (*argv[1])
    {
        case 'e' : case 'E' : Compression(argv); break;
        case 'd' : case 'D' : Decompression(argv); break;
        default  : ShowUsage();
    }

    return 0;
}

}; //end class bsc_class

void BSC_compress(const char *infile, const char *outfile, const int bsize/* = BSC_BLOCK_SIZE*/) {
  bsc_class b;
  // from https://stackoverflow.com/questions/39883433/create-argc-argv-in-the-code
  std::vector<std::string> arguments = {"", "e", infile, outfile, "-b" + std::to_string(bsize), "-p", "-e1"};
  std::vector<char*> argv;
  for (const auto& arg : arguments)
    argv.push_back((char*)arg.data());
  argv.push_back(nullptr);
  b.bsc_main(argv.size() - 1, argv.data());
}

void BSC_decompress(const char *infile, const char *outfile) {
  bsc_class b;
  // from https://stackoverflow.com/questions/39883433/create-argc-argv-in-the-code
  std::vector<std::string> arguments = {"", "d", infile, outfile};
  std::vector<char*> argv;
  for (const auto& arg : arguments)
    argv.push_back((char*)arg.data());
  argv.push_back(nullptr);
  b.bsc_main(argv.size() - 1, argv.data());
}

} // namespace bsc
} // namespace spring

/*-----------------------------------------------------------*/
/* End                                               bsc.cpp */
/*-----------------------------------------------------------*/
