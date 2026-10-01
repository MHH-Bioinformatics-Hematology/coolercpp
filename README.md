# coolercpp

A C++20 library that is API-compatible with the Python
[cooler](https://github.com/open2c/cooler) library (version 0.10.2) and
file-format-compatible with it: files it writes are read by cooler into the
same data and match cooler's own output structurally, and it reads
cooler-written files into the same data.

Scope: the cool format. `Cooler` with its metadata, the `chroms`, `bins`,
`pixels` and `matrix` selectors, `create_cooler` with all of its options,
multi-resolution and single-cell files, the file level operations and the
utilities for URIs, regions, chromosome sizes and bin tables.

cooler's analysis tools stay out of the library: iterative correction
(`balance_cooler`), coarsening, zoomify and merging are operations on matrices
rather than access to the format, and a C++ copy of them would be a second
implementation to keep in step with cooler's. Programs that need them call
cooler itself or implement them on top of this library.

```cpp
#include <coolercpp/coolercpp.hpp>

coolercpp::Cooler clr("matrix.mcool::/resolutions/10000");
auto bins = clr.bins().fetch("chr1:10M-12M");                   // Table
auto m = clr.matrix({.balance = true, .sparse = true}).fetch("chr1");
coolercpp::create_cooler("out.cool", bins_table, pixel_table, {.assembly = "hg38"});
```

`docs/API_MAPPING.md` maps every Python call to its C++ form;
`docs/DEVIATIONS.md` lists reproduced cooler quirks and the few deliberate
deviations; `docs/PROVENANCE.md` credits the code this library derives from.

## Build and install

Dependencies: a C++20 compiler, CMake 3.21 or newer, and the HDF5 C library.
doctest is fetched for the tests.

```sh
cmake -S . -B build -DCMAKE_PREFIX_PATH=/path/to/hdf5/prefix
cmake --build build -j
ctest --test-dir build --output-on-failure
cmake --install build --prefix /opt/coolercpp
```

Downstream projects use either

```cmake
find_package(coolercpp 0.3 REQUIRED)
target_link_libraries(app PRIVATE coolercpp::coolercpp)
```

or

```cmake
include(FetchContent)
FetchContent_Declare(coolercpp GIT_REPOSITORY <url> GIT_TAG <tag>)
FetchContent_MakeAvailable(coolercpp)
target_link_libraries(app PRIVATE coolercpp::coolercpp)
```

The `coolercpp_consumer` test builds such a project both ways.

## Equivalence harness

`harness/` runs cases through cooler 0.10.2 and through coolercpp, compares
the outputs (integers and coordinates exactly, floats by agreement class,
error types and messages, and for written files the HDF5 object tree,
datatypes, shapes, chunking, filters, fill values, attributes and decoded
data), reads each side's files with the other side, and measures CPU time and
peak RSS of both processes.

```sh
python harness/run.py --driver build/harness/coolercpp-harness \
    --extra-data /path/to/extra/test/matrices \
    --out report/
```

The Python interpreter running `run.py` must have cooler 0.10.2, h5py, numpy,
pandas and scipy. Input files are verified against
`harness/data_manifest.json`. A case passes when all comparisons agree,
coolercpp uses no more CPU time than Python, and its peak RSS stays within the
case's declared budget (Python's peak RSS when none is declared).
