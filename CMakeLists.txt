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

add_green_dependency(green-grids)
add_green_dependency(green-h5pp)
add_green_dependency(green-ndarray)
add_green_dependency(green-params)
add_green_dependency(green-utils)

add_subdirectory(src)
add_library(GREEN::AC ALIAS ac)

find_package(MPI REQUIRED)

add_executable(ac.exe main.cpp)
target_link_libraries(ac.exe PUBLIC GREEN::UTILS GREEN::H5PP GREEN::GRIDS GREEN::PARAMS GREEN::AC MPI::MPI_CXX)

option(Build_Tests "Build tests" ON)
if (Build_Tests)
    enable_testing()
    add_subdirectory(test)
endif ()
install(TARGETS ac.exe DESTINATION bin)
install(TARGETS ac DESTINATION lib)
INSTALL(DIRECTORY ${ac_SOURCE_DIR}/src/green/ DESTINATION include/green)