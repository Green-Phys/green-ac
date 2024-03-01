# Try to find the MPFR library
# See http://www.mpfr.org/
#
# This module supports requiring a minimum version, e.g. you can do
#   find_package(MPFR 2.3.0)
# to require version 2.3.0 to newer of MPFR.
#
# Once done this will define
#
#  MPFR_FOUND - system has MPFR lib with correct version
#  MPFR_INCLUDES - the MPFR include directory
#  MPFR_LIBRARIES - the MPFR library
#  MPFR_VERSION - MPFR version

# Copyright (c) 2006, 2007 Montel Laurent, <montel@kde.org>
# Copyright (c) 2008, 2009 Gael Guennebaud, <g.gael@free.fr>
# Copyright (c) 2010 Jitse Niesen, <jitse@maths.leeds.ac.uk>
# Copyright (c) 2015 Jack Poulson, <jack.poulson@gmail.com>
# Redistribution and use is allowed according to the terms of the BSD license.

set(gmp_ver_file "gmp.h")

if ( DEFINED ENV{GMP_DIR} )
    find_path(GMP_INCLUDES NAMES gmp.h PATHS $ENV{GMP_DIR} PATH_SUFFIXES include NO_DEFAULT_PATH)
else()
    find_path(GMP_INCLUDES NAMES gmp.h PATHS $ENV{GMP_DIR} ${INCLUDE_INSTALL_DIR} PATH_SUFFIXES include)
    find_path(GMP_X86_INCLUDES NAMES gmp-x86_64.h PATHS $ENV{GMP_DIR} ${INCLUDE_INSTALL_DIR} PATH_SUFFIXES include)
    if ( GMP_X86_INCLUDES ) # Check for redhat wrappers
        set(gmp_ver_file "gmp-x86_64.h")
    endif ()
endif()

# Set GMP_FIND_VERSION to 6.0.0 if no minimum version is specified
if(NOT GMP_FIND_VERSION)
    if(NOT GMP_FIND_VERSION_MAJOR)
        set(GMP_FIND_VERSION_MAJOR 6)
    endif()
    if(NOT GMP_FIND_VERSION_MINOR)
        set(GMP_FIND_VERSION_MINOR 0)
    endif()
    if(NOT GMP_FIND_VERSION_PATCH)
        set(GMP_FIND_VERSION_PATCH 0)
    endif()
    set(GMP_FIND_VERSION
            "${GMP_FIND_VERSION_MAJOR}.${GMP_FIND_VERSION_MINOR}.${GMP_FIND_VERSION_PATCH}")
endif()

if(GMP_INCLUDES)
    # Query GMP_VERSION
    file(READ "${GMP_INCLUDES}/${gmp_ver_file}" _gmp_version_header)

    string(REGEX MATCH "define[ \t]+__GNU_MP_VERSION[ \t]+([0-9]+)"
            _gmp_major_version_match "${_gmp_version_header}")
    set(GMP_MAJOR_VERSION "${CMAKE_MATCH_1}")
    string(REGEX MATCH "define[ \t]+__GNU_MP_VERSION_MINOR[ \t]+([0-9]+)"
            _gmp_minor_version_match "${_gmp_version_header}")
    set(GMP_MINOR_VERSION "${CMAKE_MATCH_1}")
    string(REGEX MATCH "define[ \t]+__GNU_MP_VERSION_PATCHLEVEL[ \t]+([0-9]+)"
            _gmp_patchlevel_version_match "${_gmp_version_header}")
    set(GMP_PATCHLEVEL_VERSION "${CMAKE_MATCH_1}")

    set(GMP_VERSION
            ${GMP_MAJOR_VERSION}.${GMP_MINOR_VERSION}.${GMP_PATCHLEVEL_VERSION})

    # Check whether found version exceeds minimum required
    if(${GMP_VERSION} VERSION_LESS ${GMP_FIND_VERSION})
        set(GMP_VERSION_OK FALSE)
        message(STATUS "GMP version ${GMP_VERSION} found in ${GMP_INCLUDES}, "
                "but at least version ${GMP_FIND_VERSION} is required")
    else()
        set(GMP_VERSION_OK TRUE)
    endif()
endif()

if ( DEFINED ENV{GMP_DIR} )
    find_library(GMP_LIBRARIES gmp PATHS $ENV{GMP_DIR} PATH_SUFFIXES lib NO_DEFAULT_PATH)
    find_library(GMPXX_LIBRARIES gmpxx PATHS $ENV{GMP_DIR} PATH_SUFFIXES lib NO_DEFAULT_PATH)
else ()
    find_library(GMP_LIBRARIES gmp PATHS $ENV{GMP_DIR} ${LIB_INSTALL_DIR}  PATH_SUFFIXES lib)
    find_library(GMPXX_LIBRARIES gmpxx PATHS $ENV{GMP_DIR} ${LIB_INSTALL_DIR}  PATH_SUFFIXES lib)
endif()

if( NOT GMP_LIBRARIES OR NOT GMPXX_LIBRARIES )
    message(FATAL "GMP Libraries has not been found. Please set GMPDIR environment variable to the location of the GMP library.")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP DEFAULT_MSG
        GMP_INCLUDES GMP_LIBRARIES GMPXX_LIBRARIES GMP_VERSION_OK)
mark_as_advanced(GMP_INCLUDES GMP_LIBRARIES GMPXX_LIBRARIES)

if(GMP_INCLUDES AND GMP_LIBRARIES AND NOT TARGET GMP::Library)
    add_library(GMP::Library INTERFACE IMPORTED)
    set_target_properties(GMP::Library PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES ${GMP_INCLUDES}
            )
    set_target_properties(GMP::Library PROPERTIES
            INTERFACE_LINK_LIBRARIES "${GMPXX_LIBRARIES};${GMP_LIBRARIES};"
            )
endif()
