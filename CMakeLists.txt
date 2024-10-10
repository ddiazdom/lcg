cmake_minimum_required(VERSION 3.7)
project(lcg)
set(CMAKE_CXX_STANDARD 17)

set(SOURCE_FILES main.cpp
        external/xxHash-dev/xxhash.c
        lib/utils.cpp
        lib/cdt_common.cpp
)

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
