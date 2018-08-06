#include <boost/lambda/lambda.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <string>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <omp.h>
#include "bcm/bcm.h"

int main()
{
   omp_set_num_threads(8);
#pragma omp parallel
{
int tid = omp_get_thread_num();
std::string file_name = "abc."+std::to_string(tid);
std::ofstream fout(file_name);
fout << tid;
fout.close();
spring::bcm::bcm_compress(file_name.c_str(),(file_name+".bcm").c_str());
spring::bcm::bcm_decompress((file_name+".bcm").c_str(),(file_name+".decomp").c_str());
}
std::string *str_array = new std::string[2000];
uint32_t *str_lengths = new uint32_t[2000];
for(int i = 0; i < 2000; i++) {
  str_array[i] = std::string(20,'A');
  str_lengths[i] = 20;
}
spring::bcm::bcm_str_array_compress("def",str_array,2000,str_lengths,127);
std::string *str_array_1 = new std::string[2000];
spring::bcm::bcm_str_array_decompress("def",str_array_1,2000,str_lengths);
std::cout << str_array_1[0]<<"\n";
std::cout << str_array_1[1999]<<"\n";
  
   std::cout << "Hello World" << std::endl;
    using namespace boost::lambda;
    typedef std::istream_iterator<int> in;

    std::for_each(
        in(std::cin), in(), std::cout << (_1 * 3) << " " );
    boost::filesystem::path p1 = boost::filesystem::unique_path();
    boost::filesystem::create_directory(p1);
    std::cout << p1.native();
    boost::filesystem::remove(p1);
    return 0;
}
