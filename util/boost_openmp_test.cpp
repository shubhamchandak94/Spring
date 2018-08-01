#include <boost/lambda/lambda.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <omp.h>

int main()
{
    omp_set_num_threads(1);
   std::cout << "Hello World" << std::endl;
    using namespace boost::lambda;
    typedef std::istream_iterator<int> in;

    std::for_each(
        in(std::cin), in(), std::cout << (_1 * 3) << " " );
    boost::filesystem::path p1 = boost::filesystem::unique_path();
    boost::filesystem::create_directory(p1);
    std::cout << p1.native();
    boost::filesystem::remove(p1);
}
