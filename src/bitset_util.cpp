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
