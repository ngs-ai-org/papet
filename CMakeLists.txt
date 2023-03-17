cmake_minimum_required(VERSION 3.0.0)


# user defined paths
## list of directories containing libraries headers
set(INCLUDE_DIRECTORIES "/usr/local/include/")
## list of directories in which the required libraries are installed
set(LINK_DIRECTORY      "/usr/local/lib")
## path in which the final executable will be installed
set(INSTALL_DIRECTORY   "/usr/local/bin")



# project name and it is C++ only
project(papet CXX)

# contains source code
add_subdirectory(src)

# check required libraries
## pthread
find_library(THREADS_LIB pthread)
if(NOT THREADS_LIB)
  message(FATAL_ERROR "pthread library not found")
endif()
## boost_program_options
find_library(BOOST_PROG_OPT_LIB boost_program_options)
if(NOT BOOST_PROG_OPT_LIB)
  message(FATAL_ERROR "boost_program_options library not found")
endif()
## boost_serialization
find_library(BOOST_SERIALIZATION_LIB boost_serialization)
if(NOT BOOST_SERIALIZATION_LIB)
  message(FATAL_ERROR "boost_serialization library not found")
endif()
## pbbam
find_library(PBBAM_LIB pbbam)
if(NOT PBBAM_LIB)
  message(FATAL_ERROR "pbbam library not found")
endif()
## pbcopper
find_library(PBCOPPER_LIB pbcopper)
if(NOT PBCOPPER_LIB)
  message(FATAL_ERROR "pbcopper library not found")
endif()
## ngsaipp
find_library(NGSAIPP_LIB ngsaipp)
if(NOT NGSAIPP_LIB)
  message(FATAL_ERROR "ngsaipp library not found")
endif()
