#include <boost/lambda/lambda.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <string>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <omp.h>
#include "libbsc/bsc.h"

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
spring::bsc::BSC_compress(file_name.c_str(),(file_name+".bsc").c_str());
spring::bsc::BSC_decompress((file_name+".bsc").c_str(),(file_name+".decomp").c_str());
}

uint32_t num_str = 2000000;
uint32_t str_size = 200;
std::string *str_array = new std::string[num_str];
uint32_t *str_lengths = new uint32_t[num_str];
for(int i = 0; i < num_str; i++) {
  str_array[i] = std::string(str_size,'A');
  str_lengths[i] = str_size;
}
spring::bsc::BSC_str_array_compress("def",str_array,num_str,str_lengths,64);

std::string *str_array_1 = new std::string[num_str];
spring::bsc::BSC_str_array_decompress("def",str_array_1,num_str,str_lengths);
std::cout << str_array_1[0]<<"\n";
std::cout << str_array_1[num_str-1]<<"\n";
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
