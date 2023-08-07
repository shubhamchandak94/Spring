project (SPRING)
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

option(spring_optimize_for_native "Build spoa with -march=native" ON)
option(spring_optimize_for_portability "Build spoa with -msse4.1" OFF)

message(STATUS "Setting compiler flags")
set(FLAGS "${FLAGS} -g")
##set(FLAGS "${FLAGS} -Wextra")
#set(FLAGS "${FLAGS} -Wpedantic")
#set(FLAGS "${FLAGS} -Wall")
set(FLAGS "${FLAGS} -O3")
if (spring_optimize_for_portability)
    set(FLAGS "${FLAGS} -msse4.1")
elseif (spring_optimize_for_native)
    set(FLAGS "${FLAGS} -march=native")
endif ()
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

#######################################################
# build cloudflare zlib and set environment variable ZLIB_ROOT
# the build runs during configuration time so boost can use the
# correct zlib

message(STATUS "Building zlib")
message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
execute_process(
    COMMAND ./configure
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/cloudflare_zlib/
)
execute_process(
    COMMAND make clean
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/cloudflare_zlib/
)
execute_process(
    COMMAND make static
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/cloudflare_zlib/
)

set( ENV{ZLIB_ROOT} ${CMAKE_SOURCE_DIR}/cloudflare_zlib/ )

#######################################################

add_subdirectory(boost-cmake)

set(source_dir ${CMAKE_SOURCE_DIR}/src)
set(include_dir ${CMAKE_SOURCE_DIR}/src)

# spring
set(source_files ${source_files} ${source_dir}/main.cpp)
set(source_files ${source_files} ${source_dir}/spring.cpp)
set(source_files ${source_files} ${source_dir}/util.cpp)
set(source_files ${source_files} ${source_dir}/bitset_util.cpp)
set(source_files ${source_files} ${source_dir}/preprocess.cpp)
set(source_files ${source_files} ${source_dir}/encoder.cpp)
set(source_files ${source_files} ${source_dir}/reorder_compress_streams.cpp)
set(source_files ${source_files} ${source_dir}/pe_encode.cpp)
set(source_files ${source_files} ${source_dir}/reorder_compress_quality_id.cpp)
set(source_files ${source_files} ${source_dir}/decompress.cpp)
set(source_files ${source_files} ${source_dir}/call_template_functions.cpp)

# id compression
set(source_files ${source_files} ${source_dir}/id_compression/src/Arithmetic_stream.cpp)
set(source_files ${source_files} ${source_dir}/id_compression/src/compression.cpp)
set(source_files ${source_files} ${source_dir}/id_compression/src/id_compression.cpp)
set(source_files ${source_files} ${source_dir}/id_compression/src/io_functions.cpp)
set(source_files ${source_files} ${source_dir}/id_compression/src/sam_file_allocation.cpp)
set(source_files ${source_files} ${source_dir}/id_compression/src/sam_models.cpp)
set(source_files ${source_files} ${source_dir}/id_compression/src/stream_model.cpp)

# bsc
set(source_files ${source_files} ${source_dir}/libbsc/bsc.cpp)
set(source_files ${source_files} ${source_dir}/libbsc/bsc_str_array.cpp)
set(source_files ${source_files} ${source_dir}/libbsc/libbsc/adler32/adler32.cpp)
set(source_files ${source_files} ${source_dir}/libbsc/libbsc/bwt/divsufsort/divsufsort.c)
set(source_files ${source_files} ${source_dir}/libbsc/libbsc/bwt/bwt.cpp)
set(source_files ${source_files} ${source_dir}/libbsc/libbsc/coder/coder.cpp)
set(source_files ${source_files} ${source_dir}/libbsc/libbsc/coder/qlfc/qlfc.cpp)
set(source_files ${source_files} ${source_dir}/libbsc/libbsc/coder/qlfc/qlfc_model.cpp)
set(source_files ${source_files} ${source_dir}/libbsc/libbsc/filters/detectors.cpp)
set(source_files ${source_files} ${source_dir}/libbsc/libbsc/filters/preprocessing.cpp)
set(source_files ${source_files} ${source_dir}/libbsc/libbsc/libbsc/libbsc.cpp)
set(source_files ${source_files} ${source_dir}/libbsc/libbsc/lzp/lzp.cpp)
set(source_files ${source_files} ${source_dir}/libbsc/libbsc/platform/platform.cpp)

# qvz
set(source_files ${source_files} ${source_dir}/qvz/src/cluster.cpp)
set(source_files ${source_files} ${source_dir}/qvz/src/codebook.cpp)
set(source_files ${source_files} ${source_dir}/qvz/src/distortion.cpp)
set(source_files ${source_files} ${source_dir}/qvz/src/lines.cpp)
set(source_files ${source_files} ${source_dir}/qvz/src/pmf.cpp)
set(source_files ${source_files} ${source_dir}/qvz/src/quantizer.cpp)
set(source_files ${source_files} ${source_dir}/qvz/src/qv_compressor.cpp)
set(source_files ${source_files} ${source_dir}/qvz/src/qvz.cpp)
set(source_files ${source_files} ${source_dir}/qvz/src/util.cpp)
set(source_files ${source_files} ${source_dir}/qvz/src/well.cpp)

add_executable (spring ${source_files})
install(TARGETS spring DESTINATION bin)

target_include_directories(spring PRIVATE ${include_dir})
target_link_libraries(spring PUBLIC Boost::filesystem)
target_link_libraries(spring PUBLIC Boost::program_options)
target_link_libraries(spring PUBLIC Boost::iostreams)
