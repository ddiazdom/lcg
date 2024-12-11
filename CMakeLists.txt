cmake_minimum_required(VERSION 3.10)
project(lcg)
set(CMAKE_CXX_STANDARD 17)

include(CheckCXXSourceCompiles)

set(SOURCE_FILES main.cpp
        external/xxHash-dev/xxhash.c
        lib/utils.cpp
        lib/cdt_common.cpp
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -msse4.2 -Wall")

set(SSE42_CXX_CODE "
        #include <x86intrin.h>
        int main(){
        #if __SSE4_2__
        return 0;
        #else
        #error \"SSE4.2 is not supported\"
        #endif
        }
")

set(AVX2_CXX_CODE "
        #include <immintrin.h>
        int main(){
        #if __AVX2__
        return 0;
        #else
        #error \"AVX2 is not supported\"
        #endif
        }
")

set(NEON_CXX_CODE "
        #include <arm_neon.h>
        int main(){
        #if __ARM_NEON__
        return 0;
        #else
        #error \"NEON is not supported\"
        #endif
        }
")


check_cxx_source_compiles("${AVX2_CXX_CODE}" AVX2_SUPPORT)
if(AVX2_SUPPORT)
        message("AVX2 is supported")
else()
        message("AVX2 is not supported")
endif()

check_cxx_source_compiles("${NEON_CXX_CODE}" NEON_SUPPORT)
if(NEON_SUPPORT)
        message("NEON is supported")
else()
        message("NEON is not supported")
endif()

check_cxx_source_compiles("${SSE42_CXX_CODE}" SSE42_SUPPORT)
if(SSE42_SUPPORT)
        message("SSE42 is supported")
else()
        message("SSE42 is not supported")
endif()

#using malloc_count for monitoring memory usage during debug mode
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        list(APPEND SOURCE_FILES external/malloc_count/malloc_count.c)
endif()

add_executable(lcg ${SOURCE_FILES})

#add DEBUG_MODE compilation flag
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        target_compile_definitions(lcg PRIVATE DEBUG_MODE)
endif()

target_compile_options(lcg PRIVATE -Wall -Wextra -Wpedantic
        -O3
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        message(STATUS "Disabling vla-extension due to ${CMAKE_CXX_COMPILER_ID}")
        target_compile_options(lcg PUBLIC -Wno-vla-extension -Wno-undefined-var-template)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        message(STATUS "Adding -march=native compiler flag due to ${CMAKE_CXX_COMPILER_ID}")
        target_compile_options(lcg PUBLIC -march=native)
endif()

target_link_libraries(lcg LINK_PUBLIC pthread)

#malloc_count in debug mode
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        target_link_libraries(lcg LINK_PUBLIC dl)
endif()

#for the file_system functionality on Linux
if(UNIX AND NOT APPLE)
        target_link_libraries(lcg LINK_PUBLIC stdc++fs)
endif()

target_include_directories(lcg PUBLIC include external/xxHash-dev external/malloc_count)
