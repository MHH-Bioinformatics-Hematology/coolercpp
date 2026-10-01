# coolercpp

A C++20 library for the cooler contact matrix formats, following the API of the
Python [cooler](https://github.com/open2c/cooler) library (version 0.10.2) and
its files: what it writes, cooler reads into the same data, and what cooler
writes, it reads into the same data.

The library covers the format: the `Cooler` object with its metadata, the
`chroms`, `bins`, `pixels` and `matrix` selectors, `create_cooler`,
multi-resolution and single-cell files, the file level operations, and the
utilities for URIs, regions, chromosome sizes and bin tables.

cooler's analysis tools stay out of the library. Iterative correction,
coarsening, zoomify and merging are operations on matrices rather than access to
the format, and a C++ copy of them would be a second implementation to keep in
step with cooler's. Programs that need them call cooler itself or build them on
top of this library.

```cpp
#include <iostream>
#include <coolercpp/coolercpp.hpp>

int main() {
    const coolercpp::Cooler clr("matrix.mcool::/resolutions/10000");
    const coolercpp::Table bins = clr.bins().fetch(coolercpp::Region{"chr1:10M-12M"});
    const coolercpp::MatrixResult matrix =
        clr.matrix({.balance = true, .sparse = true}).fetch(coolercpp::Region{"chr1"});
    std::cout << bins.num_rows() << " bins, " << matrix.sparse().nnz() << " pixels\n";
    return 0;
}
```

The programs under [`examples/`](https://github.com/MHH-Bioinformatics-Hematology/coolercpp/tree/main/examples)
are complete and build with the library: `read_cool` reads a file,
`create_cool` writes one, and `scool_fileops` writes a single-cell file and
operates on its groups. The pages that follow show them in full.
