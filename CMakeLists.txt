cmake_minimum_required(VERSION 3.0)
project(tdigestcpp)

set(CMAKE_CXX_STANDARD 11)
option(BUILD_TEST "Build tests" OFF)
add_library(tdigestcpp SHARED src/TDigestHistogram.cpp src/TDigestHistogram.h)
target_include_directories(tdigestcpp PRIVATE src)

if (CMAKE_BUILD_TYPE MATCHES Debug)
    add_executable(tdigest-test test/compare_test.cpp)
    target_include_directories(tdigest-test PRIVATE src)
    find_library(LIBTDIGEST libtdigest.so)
    target_link_libraries(tdigest-test tdigestcpp ${LIBTDIGEST})
ENDIF()