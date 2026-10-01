# Reading cool files

`Cooler` takes a path or a URI such as `matrix.mcool::/resolutions/10000` and
answers what the Python object answers.

`binsize()` returns `std::optional<std::int64_t>`, `nullopt` where cooler
returns `None` for a variable bin size; `info()` returns the file's attributes as
JSON; `shape()` returns the number of bins on each axis.

```cpp title="examples/read_cool.cpp"
--8<-- "examples/read_cool.cpp"
```

```
$ read_cool out.cool
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

```cpp title="examples/selectors.cpp"
--8<-- "examples/selectors.cpp"
```

```
$ selectors out.cool chrA chrB
chroms: chrA(25) chrB(17)
bins of chrA: 3, columns chrom start end
pixels of chrA: 4, first chrA:0 to chrA:0 = 5
dense 3 x 3, first row: 5 3 0
sparse holds 4 entries
pixels of the pair chrA x chrB: 1, first chrA:20 to chrB:0 = 2
no weight column in this file, so nothing to balance
```

`balance` takes a bool or the name of a bin table column, so the divisive
vectors of a file written from `.hic` are applied with
`{.balance = "KR", .divisive_weights = true}`.

## Regions and chromosome sizes

Region strings follow cooler's parsing, including its humanized forms
(`10M`, `1.5kb`), and `read_chromsizes`, `binnify` and `check_bins` build and
check bin tables as `cooler.util` does.
