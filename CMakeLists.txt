project (SPRING)
cmake_minimum_required (VERSION 2.8)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (NOT UNIX)
	message(FATAL_ERROR "Non Unix platforms not supported currently.")
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	message(STATUS "Detected compiler: GNU Compiler Collection")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
	message(STATUS "Detected compiler: LLVM Clang or Apple Clang")
else ()
	message(FATAL_ERROR "Compiler currently not supported by CMake build: ${CMAKE_CXX_COMPILER_ID}")
endif()

message(STATUS "Setting compiler flags")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")


FIND_PACKAGE( OpenMP REQUIRED)
message(STATUS "OPENMP FOUND")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")


message(STATUS "** CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
message(STATUS "** CMAKE_C_FLAGS: ${CMAKE_C_FLAGS}")

add_subdirectory(boost-cmake)

set(source_dir ${CMAKE_SOURCE_DIR}/src)
set(include_dir ${CMAKE_SOURCE_DIR}/src)
set(source_files ${source_files} ${source_dir}/../util/boost_openmp_bcm_test.cpp)
set(source_files ${source_files} ${source_dir}/bcm/bcm.cpp)
set(source_files ${source_files} ${source_dir}/bcm/divsufsort.cpp)

add_executable (spring ${source_files})
target_include_directories(spring PRIVATE ${include_dir})
target_link_libraries(spring PUBLIC Boost::filesystem)
