project(_test)

Include(FetchContent)

FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.4.0 # or a later release_
)

FetchContent_MakeAvailable(Catch2)
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)

add_executable(ac_test ac_test.cpp
        gmp_float_test.cpp)
target_link_libraries(ac_test
        PRIVATE
        Catch2::Catch2WithMain
        GREEN::AC)

include(CTest)
include(Catch)
catch_discover_tests(ac_test)