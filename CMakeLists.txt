cmake_minimum_required (VERSION 3.5.1)
project (SPRING)

if (UNIX)
  message(STATUS "Setting GCC flags")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
else()
  message(FATAL_ERROR "Non Unix platforms not supported currently.")
endif()

message(STATUS "** CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")

set(boost_root_dir ${CMAKE_CURRENT_BINARY_DIR}/boost)
set(Boost_INCLUDE_DIR ${boost_root_dir}/include)
set(Boost_LIBRARY_DIRS ${boost_root_dir}/lib)
if (UNIX)
	include(ExternalProject)

	set(BOOST_ROOT ${boost_root_dir})
	set(Boost_NO_SYSTEM_PATHS "ON")
	find_package(Boost REQUIRED COMPONENTS filesystem system)

	if (NOT Boost_FOUND)
	    set(boost_git_repository https://github.com/boostorg/boost.git)
	    set(boost_git_tag "boost-1.67.0")
	    ExternalProject_Add(Boost
		PREFIX ${boost_root_dir}
		GIT_REPOSITORY ${boost_git_repository}
		GIT_TAG ${boost_git_tag}
		CONFIGURE_COMMAND ./bootstrap.sh --with-libraries=filesystem,system --prefix=${boost_root_dir}
		BUILD_COMMAND ./b2 link=static install
		BUILD_IN_SOURCE 1
		INSTALL_COMMAND ""
	    )

	    add_library(boost STATIC IMPORTED DEPENDS Boost)
	else ()
	    include_directories (SYSTEM ${Boost_INCLUDE_DIR})
	    message(STATUS "Boost already installed - rebuild not required")
	endif()
else()
	message(FATAL_ERROR "Boost build currently not supported on non UNIX-like systems")
endif()

set(spring "spring")
set(source_dir ${CMAKE_SOURCE_DIR}/src)
set(include_dir ${CMAKE_SOURCE_DIR}/src)
set(source_files ${source_files} ${source_dir}/test.cpp)
add_executable (${spring} ${source_files})
target_include_directories(${spring} PRIVATE ${include_dir})
target_include_directories(${spring} PRIVATE ${Boost_INCLUDE_DIR})
target_link_libraries(${spring} "-fopenmp")
target_link_libraries(${spring} ${Boost_LIBRARY_DIRS}/libboost_filesystem.a)
target_link_libraries(${spring} ${Boost_LIBRARY_DIRS}/libboost_system.a)
