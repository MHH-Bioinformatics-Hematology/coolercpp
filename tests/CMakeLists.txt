# Unit tests (doctest). doctest's own CMake project declares a pre 3.5 minimum
# that CMake 4 rejects, so only its sources are fetched and the single header
# is exposed through an interface target.

include(FetchContent)
FetchContent_Declare(doctest
    GIT_REPOSITORY https://github.com/doctest/doctest.git
    GIT_TAG v2.4.11
    GIT_SHALLOW TRUE
    SOURCE_SUBDIR do-not-add)
FetchContent_MakeAvailable(doctest)

add_library(coolercpp_doctest INTERFACE)
target_include_directories(coolercpp_doctest SYSTEM INTERFACE "${doctest_SOURCE_DIR}")

add_executable(coolercpp_tests
    test_main.cpp
    test_foundation.cpp
    test_region.cpp
    test_cooler.cpp
    test_create.cpp)
target_link_libraries(coolercpp_tests PRIVATE coolercpp::coolercpp coolercpp_doctest)
# Tests read the tracked data files and write into the build tree's scratch
# directory; both paths are passed in, never derived from the working directory.
target_compile_definitions(coolercpp_tests PRIVATE
    COOLERCPP_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data"
    COOLERCPP_TEST_SCRATCH="${CMAKE_CURRENT_BINARY_DIR}/scratch")
# The tests reach into the internal numpy compatibility helpers.
target_include_directories(coolercpp_tests PRIVATE ${PROJECT_SOURCE_DIR}/src)

add_test(NAME coolercpp_unit_tests COMMAND coolercpp_tests)

# A downstream project consumes coolercpp through find_package on an installed
# tree and through FetchContent on the source tree.
add_test(NAME coolercpp_consumer
    COMMAND ${CMAKE_COMMAND}
        -DSOURCE_DIR=${PROJECT_SOURCE_DIR}
        -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/consumer
        -DPREFIX_PATH=${CMAKE_PREFIX_PATH}
        -DGENERATOR=${CMAKE_GENERATOR}
        -P ${CMAKE_CURRENT_SOURCE_DIR}/consumer/run_consumer.cmake)
set_tests_properties(coolercpp_consumer PROPERTIES TIMEOUT 900)
