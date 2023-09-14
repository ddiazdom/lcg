cmake_minimum_required(VERSION 3.7)  # CMake version check
project(text_parser)                 # Create project "simple_example"
set(CMAKE_CXX_STANDARD 17)           # Enable c++17 standard

# Add main.cpp file of project root directory as source file
set(SOURCE_FILES main.cpp
        external/xxHash-dev/xxhash.c
        lib/utils.cpp
        lib/cdt_common.cpp
)

# Add executable target with source files listed in SOURCE_FILES variable
add_executable(parallel_hasher_test ${SOURCE_FILES})
target_compile_options(parallel_hasher_test PRIVATE -Wall -Wextra -Wpedantic -O3
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        message(STATUS "Disabling vla-extension due to ${CMAKE_CXX_COMPILER_ID}")
        target_compile_options(parallel_hasher_test PUBLIC -Wno-vla-extension -Wno-undefined-var-template)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        message(STATUS "Adding -march=native compiler flag due to ${CMAKE_CXX_COMPILER_ID}")
        target_compile_options(parallel_hasher_test PUBLIC -march=native)
endif()

target_link_libraries(parallel_hasher_test LINK_PUBLIC pthread)
#for the file_system functionality on Linux
if(UNIX AND NOT APPLE)
        target_link_libraries(parallel_hasher_test LINK_PUBLIC stdc++fs)
endif()

target_include_directories(parallel_hasher_test PRIVATE include)

#add_executable(fix_bug fix_bug.cpp lib/utils.cpp lib/cdt_common.cpp)
#target_include_directories(fix_bug PRIVATE include)
#target_compile_options(fix_bug PRIVATE -Wall -Wextra -Wpedantic)
#
#add_executable(rd_vector scripts/generate_random_vector.cpp)
#target_compile_options(rd_vector PRIVATE -Wall -Wextra -Wpedantic -O3 )
#
#add_executable(gen_rand_str_coll scripts/create_random_str_coll.cpp)
#target_compile_options(gen_rand_str_coll PRIVATE -Wall -Wextra -Wpedantic -O3)
#
#add_executable(test_vbyte scripts/test_vbyte.cpp)
#target_compile_options(test_vbyte PRIVATE -Wall -Wextra -Wpedantic -O3)
#
#add_executable(vb_compress scripts/vbyte_compress.cpp)
#target_compile_options(vb_compress PRIVATE -Wall -Wextra -Wpedantic #-O3 )
#
#add_executable(vb_decompress scripts/vbyte_decompress.cpp)
#target_compile_options(vb_decompress PRIVATE -Wall -Wextra -Wpedantic #-O3 )


