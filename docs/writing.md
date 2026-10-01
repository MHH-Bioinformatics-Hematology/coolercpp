# Writing cool files

`create_cooler` writes a single resolution file from a bin table and a pixel
table, or from an iterable of pixel chunks for inputs too large to hold at once.
Its options are the keyword arguments of the Python function, including the
dtypes, the metadata, the assembly, the storage mode and the HDF5 dataset
options.

```cpp
coolercpp::CreateOptions options;
options.assembly = "hg38";
options.ordered = true;
coolercpp::create_cooler("out.cool", bins, pixels, options);
```

Two fields have no counterpart in cooler and exist to make files reproducible:
`creation_date` fixes that attribute instead of writing the current time, and
`generated_by` names the writing program, which an application built on the
library passes its own identity to.

## Multi-resolution and single-cell files

A multi-resolution file is a set of cooler groups under `/resolutions`, written
by calling `create_cooler` with those URIs. `create_scool` writes the
single-cell layout, the shared bin table with one group per cell under `/cells`,
as `cooler.create_scool_cooler` does.

## Threads

Writing filters and stores HDF5 chunks on several threads, and reading decodes
them the same way. The bytes written do not depend on the number of threads.
