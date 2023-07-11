cmake_minimum_required(VERSION 3.7)
project(BioParsers)

set(CMAKE_CXX_STANDARD 17)

add_library(biopar
        lib/dna_string.cpp
        lib/fastx_handler.cpp include/kseq.h)

target_compile_options(biopar
        PUBLIC
        -Wall -Wextra -Wno-ignored-qualifiers -pedantic -Wno-deprecated-copy -Wno-unused-function
        -O3 -funroll-loops -fomit-frame-pointer -ffast-math
        )

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message(STATUS "Disabling vla-extension due to ${CMAKE_CXX_COMPILER_ID}")
    target_compile_options(biopar PUBLIC -Wno-vla-extension -Wno-undefined-var-template)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    #message(STATUS "Adding -march=native compiler flag due to ${CMAKE_CXX_COMPILER_ID}")
    #target_compile_options(biopar PUBLIC -march=native)
endif()

target_link_libraries(biopar LINK_PUBLIC z cdt)

# there seems to be problem with the msse4.2 compiler flag and the new Apple chips
if(NOT CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "arm64")
    target_compile_options(biopar PUBLIC -msse4.2)
else()
    message(STATUS "Disabling SSE4.2 instructions due to conflict with host system processor ${CMAKE_HOST_SYSTEM_PROCESSOR}")
endif()

target_include_directories(biopar PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)