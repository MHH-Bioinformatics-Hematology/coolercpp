# The C++ side of the Python-against-C++ equivalence harness. The Python side
# (oracle.py, compare.py, run.py) drives it; see harness/README.md.

add_executable(coolercpp-harness driver.cpp npy.cpp)
target_link_libraries(coolercpp-harness PRIVATE coolercpp::coolercpp)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(coolercpp-harness PRIVATE -Wall -Wextra)
endif()

# Runs a command and reports its own CPU time and peak RSS (see measure.cpp).
add_executable(coolercpp-measure measure.cpp)
