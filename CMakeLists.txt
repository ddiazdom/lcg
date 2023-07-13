cmake_minimum_required(VERSION 3.7)

project(lcg)

set(CMAKE_CXX_STANDARD 17)
add_subdirectory(external/cdt)

add_library(lcg
            #OBJECT
            lib/build_lc_grammar.cpp
            external/malloc_count-master/malloc_count.c
        )

target_compile_options(lcg
        PRIVATE
        -Wall -Wextra -Wno-ignored-qualifiers -pedantic
        #-O3 -funroll-loops -fomit-frame-pointer -ffast-math
        INTERFACE
        -Wshadow)

#compiler-dependent flags
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message(STATUS "Disabling vla-extension due to ${CMAKE_CXX_COMPILER_ID}")
    target_compile_options(lcg PUBLIC -Wno-vla-extension -Wno-undefined-var-template)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(STATUS "Adding -march=native compiler flag due to ${CMAKE_CXX_COMPILER_ID}")
    target_compile_options(lcg PUBLIC -march=native)
endif()

# there seems to be problem with the msse4.2 compiler flag and the new Apple chips
if(NOT CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "arm64")
    target_compile_options(lcg PUBLIC -msse4.2)
else()
    message(STATUS "Disabling SSE4.2 instructions due to conflict with host system processor ${CMAKE_HOST_SYSTEM_PROCESSOR}")
endif()

#for the file_system functionality on Linux
if(UNIX AND NOT APPLE)
    target_link_libraries(lcg LINK_PUBLIC stdc++fs)
endif()

target_link_libraries(lcg LINK_PUBLIC z pthread cdt dl)

target_include_directories(lcg
        PUBLIC
        ${LIBSDSL_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/external
        )

target_include_directories(lcg SYSTEM PUBLIC ${LIBSDSL_INCLUDE_DIRS} external/malloc_count-master)
add_executable(lcg-cli main.cpp)
target_link_libraries(lcg-cli lcg)
target_include_directories(lcg-cli PUBLIC external/bioparsers/include)