# - Try to find LibSDSL
# Once done this will define
#  LIBSDSL_FOUND - System has SDSL
#  LIBSDSL_INCLUDE_DIRS - The SDSL include directories
#  LIBSDSL_LIBRARIES - The libraries needed to use SDSL
#  LIBSDSL_DEFINITIONS - Compiler switches required for using SDSL

find_package(PkgConfig)

pkg_check_modules(PC_LIBSDSL QUIET libsdsl)

set(LIBSDSL_DEFINITIONS ${PC_LIBSDSL_CFLAGS_OTHER})

find_path(LIBSDSL_INCLUDE_DIR sdsl/csa_sada.hpp
          HINTS ${PC_LIBSDSL_INCLUDEDIR}
                ${PC_LIBSDSL_INCLUDE_DIRS}
          PATHS
                $ENV{HOME}/usr/include
                $ENV{HOME}/include
          PATH_SUFFIXES libsdsl)

find_library(LIBSDSL_LIBRARY NAMES sdsl
             HINTS ${PC_LIBSDSL_LIBDIR}
                   ${PC_LIBSDSL_LIBRARY_DIRS}
             PATHS
                   $ENV{HOME}/usr/lib
                   $ENV{HOME}/lib)

find_library(LIBDIVSUF_LIBRARY NAMES divsufsort
        HINTS ${PC_LIBSDSL_LIBDIR}
              ${PC_LIBSDSL_LIBRARY_DIRS}
        PATHS
              $ENV{HOME}/usr/lib
              $ENV{HOME}/lib)

find_library(LIBDIVSUF64_LIBRARY NAMES divsufsort64
        HINTS ${PC_LIBSDSL_LIBDIR}
        ${PC_LIBSDSL_LIBRARY_DIRS}
        PATHS
        $ENV{HOME}/usr/lib
        $ENV{HOME}/lib)

include(FindPackageHandleStandardArgs)

# handle the QUIETLY and REQUIRED arguments and set LIBSDSL_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(LibSDSL DEFAULT_MSG
        LIBSDSL_LIBRARY LIBDIVSUF_LIBRARY LIBDIVSUF64_LIBRARY LIBSDSL_INCLUDE_DIR)

mark_as_advanced(LIBSDSL_LIBRARY LIBDIVSUF_LIBRARY LIBDIVSUF64_LIBRARY LIBSDSL_INCLUDE_DIR)

set(LIBSDSL_LIBRARIES ${LIBSDSL_LIBRARY} ${LIBDIVSUF_LIBRARY} ${LIBDIVSUF64_LIBRARY})
set(LIBSDSL_INCLUDE_DIRS ${LIBSDSL_INCLUDE_DIR} )
