cmake_minimum_required(VERSION 3.11)
set(targetName "BayesOmics")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -O3 -lboost_system -lboost_thread")
# Set C++14 standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
project(${targetName})

# Set build type to RELEASE by default:
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type provided, default to RELEASE.")
  set(CMAKE_BUILD_TYPE "RELEASE")
endif()

# Options can be turned on and off by providing cmake with a -DOPTIONNAME=ON (or OFF) argument.
# For example, to turn off benchmarks, build a shared library, and use sanitize flags in a DEBUG build:
#    cmake .. -DCMAKE_BUILD_TYPE=DEBUG -DGOOGLE_BENCH=OFF -DBENCH_STATIC=OFF -DIS_SANITIZE=ON
option(GOOGLE_TEST      "Build tests" ON)
option(GOOGLE_BENCH     "Build benchmarks" OFF)
option(LIB_STATIC    "Builds a static library instead of a shared one" OFF)
option(IS_SANITIZE  "Adds sanitive flags" OFF)



if(LIB_STATIC)
  set(CUR_LIB_NAME ${targetName}) ## current library name
  set(CUR_LIB_TYPE STATIC)
else()
  set(CUR_LIB_NAME ${targetName})
  set(CUR_LIB_TYPE SHARED)
endif()

# # Detect operating system
# if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
#     # For Linux, compile static software
#     set(BUILD_SHARED_LIBS OFF)
#     # Set any additional flags you might need for static linking
#     set(CMAKE_EXE_LINKER_FLAGS "-static")
# elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
#     # For macOS (Darwin), compile dynamic software
#     set(BUILD_SHARED_LIBS ON)
#     # Set any macOS-specific flags if needed
# endif()



include(FetchContent)


if(IS_SANITIZE)
  set(CMAKE_CXX_FLAGS "-fsanitize=address -fsanitize=undefined ${CMAKE_CXX_FLAGS}")
endif()

add_subdirectory(src)

if(GOOGLE_BENCH)
  add_subdirectory(bench)
endif()

message(STATUS "")
message(STATUS "BUILD SUMMARY")
message(STATUS "  CMAKE_GENERATOR      : ${CMAKE_GENERATOR}")
message(STATUS "  Compiler ID          : ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "  Build type           : ${CMAKE_BUILD_TYPE}")
message(STATUS "  Library name         : ${CUR_LIB_NAME}")
message(STATUS "  Library type         : ${CUR_LIB_TYPE}")
message(STATUS "  Build tests          : ${GOOGLE_TEST}")
message(STATUS "  Build benchmarks     : ${GOOGLE_BENCH}")
message(STATUS "  Sanitize flags       : ${IS_SANITIZE}")
message(STATUS "")
