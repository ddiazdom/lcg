cmake_minimum_required(VERSION 3.7)  # CMake version check
project(lcg)                 # Create project "simple_example"
set(CMAKE_CXX_STANDARD 17)           # Enable c++17 standard

# Add main.cpp file of project root directory as source file
set(SOURCE_FILES main.cpp
        external/xxHash-dev/xxhash.c
        lib/utils.cpp
        lib/cdt_common.cpp
        external/malloc_count/malloc_count.c
)

# Add executable target with source files listed in SOURCE_FILES variable
add_executable(lcg ${SOURCE_FILES})
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

target_link_libraries(lcg LINK_PUBLIC pthread dl)
#for the file_system functionality on Linux
if(UNIX AND NOT APPLE)
        target_link_libraries(lcg LINK_PUBLIC stdc++fs)
endif()

target_include_directories(lcg PUBLIC include external/xxHash-dev external/malloc_count)