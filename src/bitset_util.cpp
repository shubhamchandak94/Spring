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

#include "bitset_util.h"
#include "params.h"

namespace spring {

void bbhashdict::findpos(int64_t *dictidx, const uint64_t &startposidx) {
  dictidx[0] = startpos[startposidx];
  auto endidx = startpos[startposidx + 1];
  if (read_id[endidx - 1] ==
      MAX_NUM_READS)  // means exactly one read has been removed
    dictidx[1] = endidx - 1;
  else if (read_id[endidx - 1] ==
           MAX_NUM_READS + 1)  // means two or more reads have
                               // been removed (in this case
                               // second last entry stores the
                               // number of reads left)
    dictidx[1] = dictidx[0] + read_id[endidx - 2];
  else
    dictidx[1] = endidx;  // no read deleted
  return;
}

void bbhashdict::remove(int64_t *dictidx, const uint64_t &startposidx,
                        const int64_t current) {
  auto size = dictidx[1] - dictidx[0];
  if (size == 1)  // just one read left in bin
  {
    empty_bin[startposidx] = 1;
    return;  // need to keep one read to check during matching
  }
  int64_t pos =
      std::lower_bound(read_id + dictidx[0], read_id + dictidx[1], current) -
      (read_id + dictidx[0]);

  for (int64_t i = dictidx[0] + pos; i < dictidx[1] - 1; i++)
    read_id[i] = read_id[i + 1];
  auto endidx = startpos[startposidx + 1];
  if (dictidx[1] == endidx)  // this is first read to be deleted
    read_id[endidx - 1] = MAX_NUM_READS;
  else if (read_id[endidx - 1] ==
           MAX_NUM_READS)  // exactly one read has been deleted till now
  {
    read_id[endidx - 1] = MAX_NUM_READS + 1;
    read_id[endidx - 2] = (uint32_t)(size - 1);  // number of reads left in bin
  } else  // more than two reads have been deleted
    read_id[endidx - 2]--;

  return;
}

}  // namespace spring
