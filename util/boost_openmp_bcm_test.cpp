#include <boost/lambda/lambda.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <string>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <omp.h>
#include "bcm/bcm.cpp"

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
spring::bcm::bcm_class b;
b.bcm_main(file_name.c_str(),(file_name+".bcm").c_str());
spring::bcm::bcm_class b1;
b1.bcm_main((file_name+".bcm").c_str(),(file_name+".decomp").c_str(),true);
}	
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
