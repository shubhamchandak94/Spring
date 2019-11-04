project (SPRING_REORDER)
cmake_minimum_required (VERSION 3.9)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include(CheckIPOSupported)

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

check_ipo_supported(RESULT result)
if(result)
  message(STATUS "Support for Interprocedural Optimization detected.")
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

message(STATUS "Setting compiler flags")
#set(FLAGS "${FLAGS} -g")
#set(FLAGS "${FLAGS} -Wextra")
#set(FLAGS "${FLAGS} -Wpedantic")
#set(FLAGS "${FLAGS} -Wall")
set(FLAGS "${FLAGS} -O3")
set(FLAGS "${FLAGS} -march=native")
set(FLAGS "${FLAGS} -fomit-frame-pointer")
set(FLAGS "${FLAGS} -fstrict-aliasing")
set(FLAGS "${FLAGS} -ffast-math")

set(CMAKE_CXX_FLAGS "${FLAGS}")
set(CMAKE_C_FLAGS "${FLAGS}")

FIND_PACKAGE( OpenMP REQUIRED)
message(STATUS "OPENMP FOUND")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")

message(STATUS "** CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
message(STATUS "** CMAKE_C_FLAGS: ${CMAKE_C_FLAGS}")

add_subdirectory(boost-cmake)

set(source_dir ${CMAKE_SOURCE_DIR}/src)
set(include_dir ${CMAKE_SOURCE_DIR}/src)

# spring-reorder
set(source_files ${source_files} ${source_dir}/main.cpp)
set(source_files ${source_files} ${source_dir}/spring.cpp)
set(source_files ${source_files} ${source_dir}/util.cpp)
set(source_files ${source_files} ${source_dir}/bitset_util.cpp)
set(source_files ${source_files} ${source_dir}/preprocess.cpp)
set(source_files ${source_files} ${source_dir}/encoder.cpp)
set(source_files ${source_files} ${source_dir}/pe_encode.cpp)
set(source_files ${source_files} ${source_dir}/generate_reordered_fastq.cpp)
set(source_files ${source_files} ${source_dir}/call_template_functions.cpp)

add_executable (spring-reorder ${source_files})

target_include_directories(spring-reorder PRIVATE ${include_dir})
target_link_libraries(spring-reorder PUBLIC Boost::filesystem)
target_link_libraries(spring-reorder PUBLIC Boost::program_options)
target_link_libraries(spring-reorder PUBLIC Boost::iostreams)
