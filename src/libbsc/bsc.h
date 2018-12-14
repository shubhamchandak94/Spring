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

#ifndef SPRING_LIBBSC_BSC_H_
#define SPRING_LIBBSC_BSC_H_

#include "params.h"

namespace spring {
namespace bsc {

void BSC_compress(const char *infile, const char *outfile,
                  const int bsize = BSC_BLOCK_SIZE);

void BSC_decompress(const char *infile, const char *outfile);

void BSC_str_array_compress(const char *outfile, std::string *str_array_param,
                            const uint32_t size_str_array_param,
                            uint32_t *str_lengths_param,
                            const int bsize = BSC_BLOCK_SIZE);

void BSC_str_array_decompress(const char *infile, std::string *str_array_param,
                              const uint32_t size_str_array_param,
                              uint32_t *str_lengths_param);

}  // namespace bsc
}  // namespace spring

#endif  // SPRING_LIBBSC_BSC_H_

/*-----------------------------------------------------------*/
/* End                                               bsc.cpp */
/*-----------------------------------------------------------*/
