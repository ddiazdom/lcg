cmake_minimum_required(VERSION 3.7)
project(lcg)
set(CMAKE_CXX_STANDARD 17)

set(SOURCE_FILES main.cpp
        external/xxHash-dev/xxhash.c
        lib/utils.cpp
        lib/cdt_common.cpp
        lib/input_reader.cpp
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

# Find htslib via pkg-config or use HTSLIB_ROOT hint
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
    pkg_check_modules(HTSLIB QUIET htslib)
endif()

if(NOT HTSLIB_FOUND)
    # Allow user to specify HTSLIB_ROOT, or fall back to $HOME/htslib-local
    if(NOT DEFINED HTSLIB_ROOT)
        set(HTSLIB_ROOT "$ENV{HOME}/htslib-local")
    endif()
    if(EXISTS "${HTSLIB_ROOT}/include/htslib/bgzf.h")
        set(HTSLIB_INCLUDE_DIRS "${HTSLIB_ROOT}/include")
        set(HTSLIB_LIBRARY_DIRS "${HTSLIB_ROOT}/lib")
        find_library(HTSLIB_LIB hts PATHS "${HTSLIB_ROOT}/lib" NO_DEFAULT_PATH)
        if(NOT HTSLIB_LIB)
            message(FATAL_ERROR "htslib library not found in ${HTSLIB_ROOT}/lib")
        endif()
        set(HTSLIB_LIBRARIES "${HTSLIB_LIB}")
        set(HTSLIB_FOUND TRUE)
    else()
        message(FATAL_ERROR "htslib not found. Install libhts-dev or set HTSLIB_ROOT.")
    endif()
endif()

target_include_directories(lcg PUBLIC include external/xxHash-dev external/malloc_count)
if(HTSLIB_INCLUDE_DIRS)
    target_include_directories(lcg PUBLIC ${HTSLIB_INCLUDE_DIRS})
endif()
if(HTSLIB_LIBRARY_DIRS)
    target_link_directories(lcg PUBLIC ${HTSLIB_LIBRARY_DIRS})
endif()
target_link_libraries(lcg LINK_PUBLIC ${HTSLIB_LIBRARIES} z)
