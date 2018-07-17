#include "algorithms/SPRING/qvz/include/util.h"

namespace spring {
namespace qvz {

/**
 * Finds the ceiling of the log2 of a number iteratively
 */
int cb_log2(int x) {
  int res = 0;
  int x2 = x;

  while (x2 > 1) {
    x2 >>= 1;
    res += 1;
  }

  if ((1 << res) == x) return res;
  return res + 1;
}

} // namespace qvz
} // namespace spring
