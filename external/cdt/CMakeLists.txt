cmake_minimum_required(VERSION 3.7)
project(CDT)

set(CMAKE_CXX_STANDARD 17)

add_library(cdt
        lib/cdt_common.cpp
        lib/utils.cpp
        include/xxHash-dev/xxhash.c)

target_compile_options(cdt
        PUBLIC
        -Wall -Wextra -Wno-ignored-qualifiers -pedantic -Wno-deprecated-copy -Wno-unused-function
        -O3 -funroll-loops -fomit-frame-pointer -ffast-math
        )

# this extra flags only applies to clang compilers
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message(STATUS "Disabling vla-extension due to ${CMAKE_CXX_COMPILER_ID}")
    target_compile_options(cdt PUBLIC -Wno-vla-extension -Wno-undefined-var-template)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(STATUS "Adding -march=native compiler flag due to ${CMAKE_CXX_COMPILER_ID}")
    target_compile_options(cdt PUBLIC -march=native)
endif()

# there seems to be problem with the msse4.2 compiler flag and the new Apple chips
if(NOT CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "arm64")
    target_compile_options(cdt PUBLIC -msse4.2)
else()
    message(STATUS "Disabling SSE4.2 instructions due to conflict with host system processor ${CMAKE_HOST_SYSTEM_PROCESSOR}")
endif()

target_link_libraries(cdt LINK_PUBLIC ${LIBSDSL_LIBRARIES})

target_include_directories(cdt PUBLIC
                           ${CMAKE_CURRENT_SOURCE_DIR}/include
                           ${CMAKE_CURRENT_SOURCE_DIR}/external)

target_include_directories(cdt SYSTEM PUBLIC ${LIBSDSL_INCLUDE_DIRS})
