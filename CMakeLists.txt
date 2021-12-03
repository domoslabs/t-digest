cmake_minimum_required(VERSION 3.0)
project(tdigestcpp)

set(CMAKE_CXX_STANDARD 11)
add_library(tdigestcpp SHARED src/libtdigestcpp.cpp include/libtdigestcpp.h)
target_include_directories(tdigestcpp PRIVATE src)

set_target_properties(tdigestcpp PROPERTIES
        PUBLIC_HEADER include/libtdigestcpp.h)

if (CMAKE_BUILD_TYPE MATCHES Debug)
    add_executable(tdigest-test test/compare_test.cpp)
    target_include_directories(tdigest-test PRIVATE src)
    find_library(LIBTDIGEST libtdigest.so)
    target_link_libraries(tdigest-test tdigestcpp ${LIBTDIGEST})
ENDIF()

install(TARGETS tdigestcpp
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})