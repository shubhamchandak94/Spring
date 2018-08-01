#ifndef SPRING_PARAMS_H_
#define SPRING_PARAMS_H_ 

#include <string>

namespace spring {

const int MAX_READ_LEN = 511;
const uint32_t MAX_NUM_READS = 4294967290;
const int NUM_DICT_REORDER = 2;
const int MAX_SEARCH_REORDER = 1000;
const int THRESH_REORDER = 4;
const int NUM_LOCKS_REORDER = 0x1000000;// limits on number of locks (power of 2 for fast mod)
const int NUM_DICT_ENCODER = 2;
const int MAX_SEARCH_ENCODER = 1000;
const int THRESH_ENCODER = 24;

} // namespace spring

#endif // SPRING_PARAMS_H_
