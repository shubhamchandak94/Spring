#ifndef SPRING_UTIL_H_
#define SPRING_UTIL_H_

#include <fstream>
#include <string>

namespace spring {

void reverse_complement(char *s, char *s1, int readlen,
                        char chartorevchar[128]);

std::string reverse_complement(std::string s, int readlen,
                               char chartorevchar[128]);

}  // namespace spring

#endif  // SPRING_UTIL_H_
