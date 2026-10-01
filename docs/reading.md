# Reading cool files

`Cooler` takes a path or a URI such as `matrix.mcool::/resolutions/10000` and
answers what the Python object answers.

```cpp
coolercpp::Cooler clr("matrix.cool");
clr.binsize();                     // std::optional<std::int64_t>
clr.chromnames();
clr.info();                        // the file's attributes as JSON
clr.shape();
```

The whole of `examples/read_cool.cpp`, which prints that metadata and then a bin
table and a matrix region of the first chromosome:

```cpp title="examples/read_cool.cpp"
--8<-- "examples/read_cool.cpp"
```

```
$ read_cool matrix.cool
bin size: 10
chromosomes: 2, bins: 5
non-zero pixels: 5
bins on chrA: 3 (chrom start end )
matrix 3 x 3, first row: 5 3 0
```

## Selectors

The `chroms`, `bins` and `pixels` selectors return a `Table`, the C++
counterpart of the pandas DataFrame cooler returns; `matrix` returns a dense or
sparse matrix.

```cpp
auto bins = clr.bins().fetch("chr1:10M-12M");
auto pixels = clr.pixels(true).fetch("chr1", "chr2");        // join = true
auto dense = clr.matrix({.balance = false}).fetch("chr1");
auto sparse = clr.matrix({.balance = true, .sparse = true}).fetch("chr1");
```

`balance` takes a bool or the name of a bin table column, so the divisive
vectors of a file written from `.hic` are applied with
`{.balance = "KR", .divisive_weights = true}`.

## Regions and chromosome sizes

Region strings follow cooler's parsing, including its humanized forms
(`10M`, `1.5kb`), and `read_chromsizes`, `binnify` and `check_bins` build and
check bin tables as `cooler.util` does.
