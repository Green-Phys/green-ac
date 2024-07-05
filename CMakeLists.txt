cmake_minimum_required(VERSION 3.18)
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
project(ac)

if (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
    message(FATAL_ERROR "In-source builds are not allowed")
endif ()

################ GREEN RELEASE VERSION #####################
set(GREEN_RELEASE "v0.2.4")

set(CMAKE_CXX_STANDARD 17)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

include(GreenDeps)

add_green_dependency(green-ndarray)

add_subdirectory(src)
add_library(GREEN::AC ALIAS ac)

option(Build_Tests "Build tests" ON)
if (Build_Tests)
    enable_testing()
    add_subdirectory(test)
endif ()

if( ${GREEN_PYTHON} )
    set_property(TARGET ac PROPERTY POSITION_INDEPENDENT_CODE ON)
    set(PYBIND11_FINDPYTHON ON)
    find_package(pybind11 CONFIG REQUIRED)

    message(STATUS "Building python packaging")

    pybind11_add_module(_green_ac python/main.cpp python/nevanlinna.cpp)
    target_link_libraries(_green_ac PUBLIC GREEN::AC)

    # EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
    # define (VERSION_INFO) here.
    target_compile_definitions(_green_ac
                           PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
    set_property(TARGET _green_ac PROPERTY POSITION_INDEPENDENT_CODE ON)

    install(TARGETS _green_ac DESTINATION green_ac)
    message(STATUS ${CMAKE_INSTALL_PREFIX})
    set(CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/green_ac)
else()
    find_package(MPI REQUIRED)
    add_green_dependency(green-grids)
    add_green_dependency(green-h5pp)
    add_green_dependency(green-params)
    add_green_dependency(green-utils)
    add_executable(ac.exe main.cpp)
    target_link_libraries(ac.exe PUBLIC GREEN::UTILS GREEN::H5PP GREEN::GRIDS GREEN::PARAMS GREEN::AC MPI::MPI_CXX)
    message(STATUS ${CMAKE_INSTALL_PREFIX})
    install(TARGETS ac.exe DESTINATION bin)
    install(TARGETS ac DESTINATION lib)
    INSTALL(DIRECTORY ${ac_SOURCE_DIR}/src/green/ DESTINATION include/green)
endif()
