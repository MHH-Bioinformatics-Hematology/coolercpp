# Reproduced cooler behaviour and deliberate deviations

coolercpp reproduces what cooler 0.10.2 does, including behaviour that looks
like a bug. Each item below is pinned by a harness case (named in brackets)
that runs the same call through cooler and through coolercpp.

## cooler behaviour reproduced as is

1. `Cooler.pixels(join=True)` with a single column name raises
   `AttributeError: 'Series' object has no attribute 'columns'`, because
   `annotate` is handed a Series. [S.pixels.join.series_count_bug]
2. Divisive weights are assumed by default for balance columns named `KR`,
   `VC` and `VC_SQRT`. The docstring says `SQRT_VC`; a column of that name is
   applied multiplicatively unless `divisive_weights=True` is passed.
   [G.matrix.VC_SQRT.sparse.chr1, G.matrix.VC.multiplicative.dense.trans]
3. A slice start further below zero than the table length is not clamped for
   the index: `bins()[-100000:]` returns every row labelled from -96617.
   For the matrix selector the same slice produces a matrix of the requested,
   unclamped shape, here 1000000 by 5. [S.bins.slice.negative_beyond_start,
   S.matrix.slice.negative_beyond_start]
4. A matrix slice whose stop lies beyond the number of bins yields a matrix
   larger than the cooler, padded with empty rows and columns.
   [S.matrix.slice.beyond_end]
5. A reversed slice (`[10:5]`) raises scipy's "'shape' elements cannot be
   negative". [S.matrix.slice.reversed_rows]
6. The error for a missing balance column lacks a space:
   "No column 'bins/weight'found. ...". [S.matrix.balance_column_missing]
7. Region strings: tokens after the end coordinate are ignored
   (`"chr2L:1-2-3"` is `chr2L:1-2`); a fractional coordinate without a unit
   raises (`"1.5-2"`), with a unit it does not (`"1.5M"`); sequence names
   containing a colon cannot be addressed as strings. [S.extent,
   regions.parse.with_chromsizes, regions.parse.without_chromsizes]
8. `balance=""` disables balancing instead of naming a column.
   [S.matrix.balance_empty_string]
9. `h5opts={"maxshape": ...}` raises `TypeError` only when the pixel datasets
   are prepared, after the chromosome and bin tables have been written.
   [create.S.h5opts.maxshape_conflict]
10. `bins()` with a single field whose name merely contains "chrom" and holds
    integers converts that field to chromosome names (a substring test on a
    str). Reproduced in code; no test file carries such a column.
11. `create_cooler` accumulates the `sum` attribute from the input `count`
    column even when `count` is not among the stored `columns`.
    [create.S.sum_without_count_column]

## Deliberate deviations

1. The LZF filter (`h5opts={"compression": "lzf"}`) is h5py's own HDF5
   plugin and not part of libhdf5. coolercpp raises `ValueError` instead of
   writing LZF-compressed datasets. SZIP is available only with its default
   options.
2. Messages of errors raised inside h5py or libhdf5 (a file that cannot be
   opened, a dataset that does not exist) are coolercpp's own; the exception
   type matches. Harness cases for these compare the type only.
3. Integers beyond 64 bits (in region strings or JSON metadata) are not
   representable; coolercpp raises `ValueError` or keeps the nearest double.
4. Categorical codes are int32; pandas picks int8 or int16 for small category
   counts. The code values are identical, and the enum datasets cooler writes
   from categorical bin columns use the pandas width.
5. `create_cooler` has no `lock` argument and accepts no dask frames.
   `CreateOptions::creation_date` and `CreateOptions::generated_by` exist in
   addition, for reproducible files.
6. Groups are created without object modification times (h5py records them on
   groups), so that two runs with the same `creation_date` produce identical
   bytes. This only affects HDF5 object header bytes, not the object tree,
   datatypes, layouts, filters, attributes or data.
7. `Cooler.open()` (a raw h5py handle) has no counterpart.
