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
12. `balance_cooler` reports a run in which every bin was filtered out as
    converged: the marginal has no nonzero entry, so the loop stops with
    `var = 0.0`, which is below any tolerance, and `scale` is the mean of an
    empty array, that is NaN. Every weight is NaN. With cooler's default
    `min_nnz=10` this is what a small matrix gives.
    [balance.T.toy.default]
13. `scale` is the mean of the nonzero marginals of the last iteration before
    they were normalized, read after the loop has ended, not the scale factor
    of the weights that were returned. [balance.T.toy.unfiltered]
14. `max_iters=0` warns about the iteration limit and then raises
    `UnboundLocalError: cannot access local variable 'nzmarg' where it is not
    associated with a value`, because the statement after the loop reads a
    variable only the loop body assigns.
    [balance.T.toy.error.max_iters_zero]
15. The chunk edges run one chunk past `nnz`, so the last span reaches beyond
    the pixel table, and an `nnz` that is an exact multiple of `chunksize`
    gains a trailing empty span. `cis_only` instead uses `cooler.util.partition`
    and gets exact spans. [balance.T.toy.chunksize_equals_nnz]
16. `ignore_diags=True` zeroes the main diagonal only, because `True` compares
    as 1; `ignore_diags=False` and `ignore_diags=0` both switch the filter off
    but are stored differently in the stats (a bool and an int).
    [balance.T.toy.ignore_diags_false, balance.T.toy.ignore_diags_zero]
17. The MAD-max filter divides each chromosome's marginals by the median of its
    positive ones. A chromosome without a single positive marginal is divided
    by the median of an empty array, so its marginals become NaN, and NaN
    compares false against the cutoff: those bins keep their weight instead of
    being dropped. [balance.S.default]
18. `cis_only` returns arrays for `scale`, `converged` and `var`, one entry per
    chromosome, and stores them as array attributes of the weight column. A
    chromosome with no data is marked converged, because its variance is forced
    to 0.0. `nzmarg` outlives the chromosome loop, so a chromosome whose inner
    loop does not run reuses the previous chromosome's value.
    [balance.S.cis_only, balance.S.store.cis_only]
19. `trans_only` records `cis_only: False` in the stats, so a stored trans-only
    weight column cannot be told from a genome-wide one by its attributes.
    [balance.G.trans_only]
20. The `min_nnz`, `min_count` and MAD-max bin filters run through the base
    filters only, which never zero the cis pixels. A `trans_only` run therefore
    picks its bad bins from the whole matrix, cis data included.
    [balance.G.trans_only]

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
   addition, for reproducible files and for naming the writing program.
8. The `generated-by` attribute names coolercpp (`coolercpp-<version>`)
   instead of `cooler-0.10.2`, unless the caller passes another string.
6. Groups are created without object modification times (h5py records them on
   groups), so that two runs with the same `creation_date` produce identical
   bytes. This only affects HDF5 object header bytes, not the object tree,
   datatypes, layouts, filters, attributes or data.
7. `Cooler.open()` (a raw h5py handle) has no counterpart.
9. `balance_cooler` copies `x0` instead of taking the caller's array, filling
   its NaNs with zeros and handing the same array back as the result. An `x0`
   whose length differs from the number of bins raises `ValueError`, where
   numpy raises an `IndexError` about a boolean index of the wrong length, and
   `x0` is always float64, where a float32 `x0` makes cooler run the whole
   balancing in float32.
10. `map` becomes `threads`, and `use_lock` is accepted without effect; the
    `cooler.parallel` pipeline is internal. `docs/API_MAPPING.md` explains both.
11. The convergence warning carries cooler's message text but not its category:
    `coolercpp::warn` has no categories, so the default handler prints
    `UserWarning` where Python prints `ConvergenceWarning`, a `UserWarning`
    subclass.
12. cooler reopens the file, rereads every column of the bin table and rereads
    the pixel span for every iteration of every chromosome. coolercpp opens the
    file once and reads the chromosome ids once; the pixels are still reread per
    iteration, so the memory profile is cooler's.
