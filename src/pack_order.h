#ifndef SPRING_PACK_ORDER_H_
#define SPRING_PACK_ORDER_H_ 

#include <string>

namespace spring {

int pack_order(std::string &temp_dir, bool &paired_end);
//pack order into least number of bits possible

} // namespace spring

#endif // SPRING_PACK_ORDER_H_
