project(ac-lib CXX)

add_library(ac nevanlinna.cpp)
target_include_directories(ac PUBLIC .)
find_package(Eigen3 REQUIRED)
find_package(GMP REQUIRED)
target_link_libraries(ac PUBLIC Eigen3::Eigen)
target_link_libraries(ac PUBLIC GMP::Library)
target_link_libraries(ac PUBLIC GREEN::NDARRAY)
