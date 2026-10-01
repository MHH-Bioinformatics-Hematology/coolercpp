# Installing and building

## Conda

```
conda install -c bioconda coolercpp
```

## From source

The library needs a C++20 compiler, CMake 3.21 or newer and the HDF5 C library.
doctest is fetched for the unit tests.

```
cmake -S . -B build -DCMAKE_PREFIX_PATH=/path/to/hdf5/prefix
cmake --build build -j
ctest --test-dir build --output-on-failure
cmake --install build --prefix /opt/coolercpp
```

## Using it from another project

```cmake
find_package(coolercpp 0.4 REQUIRED)
target_link_libraries(app PRIVATE coolercpp::coolercpp)
```

Or with FetchContent:

```cmake
include(FetchContent)
FetchContent_Declare(coolercpp GIT_REPOSITORY <url> GIT_TAG <tag>)
FetchContent_MakeAvailable(coolercpp)
target_link_libraries(app PRIVATE coolercpp::coolercpp)
```

The `coolercpp_consumer` test builds such a project both ways.
