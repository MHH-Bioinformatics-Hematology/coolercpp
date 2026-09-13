# Python to C++ API mapping

coolercpp reproduces the programmatic API of cooler 0.10.2. Names, parameters
and semantics follow the Python; the spelling follows C++. This table lists
what exists so far (milestone 1) and is extended with every milestone.

Conventions:

- `import cooler` corresponds to `#include <coolercpp/coolercpp.hpp>` and the
  namespace `coolercpp`.
- Keyword arguments become fields of an options struct, used with designated
  initializers: `clr.matrix({.balance = false, .sparse = true})`.
- `None` is `std::nullopt` (or `std::monostate` inside `H5OptValue`).
- A pandas `DataFrame` is a `coolercpp::Table`; a `Series` is a one column
  `Table`; `df["x"]` is a `coolercpp::Column`; a `pd.Categorical` is a Column of
  `DType::Categorical`; `df.index` is `Table::index()`.
- A numpy array is a `Column` (one dimensional) or a `DenseMatrix` (two
  dimensional, row major); a `scipy.sparse.coo_matrix` is a `SparseMatrix`.
- Python exceptions map to classes of the same name in `coolercpp/errors.hpp`
  (`KeyError`, `ValueError`, `IndexError`, `TypeError`, `AttributeError`,
  `NotImplementedError`, `OSError`, `BadInputError`), all derived from
  `coolercpp::Error`; `python_type()` returns the Python class name.
- `warnings.warn` goes through `coolercpp::warn`, redirectable with
  `set_warning_handler`.

## Cooler

| Python (cooler 0.10.2) | C++ (coolercpp) |
|---|---|
| `cooler.Cooler(store)` | `Cooler(const std::string& store)` |
| `cooler.Cooler(store, root=...)` (deprecated) | `Cooler(store, root)` |
| `c.filename`, `c.root`, `c.uri`, `c.store` | `c.filename()`, `c.root()`, `c.uri()`, `c.store()` |
| `c._refresh()` | `c.refresh()` |
| `c.info` | `json::Value c.info()` (an object; strings JSON-decoded as cooler does; numbers carry their numpy dtype in `Value::dtype()`) |
| `c.storage_mode` | `std::string c.storage_mode()` |
| `c.binsize` | `std::optional<std::int64_t> c.binsize()` |
| `c.chromsizes` | `const ChromSizes& c.chromsizes()` (`names()`, `lengths()`, `dtype()`, `find`, `operator[]`) |
| `c.chromnames` | `std::vector<std::string> c.chromnames()` |
| `c.offset(region)` | `std::int64_t c.offset(const Region&)` |
| `c.extent(region)` | `std::pair<std::int64_t, std::int64_t> c.extent(const Region&)` |
| `c.shape` | `std::array<std::int64_t, 2> c.shape()` |
| `c.chroms(**kwargs)` | `RangeSelector1D c.chroms(bool convert_enum = true)` |
| `c.bins(**kwargs)` | `RangeSelector1D c.bins(bool convert_enum = true)` |
| `c.pixels(join=False, **kwargs)` | `RangeSelector1D c.pixels(bool join = false, bool convert_enum = true)` |
| `c.matrix(field, balance, sparse, as_pixels, join, ignore_index, divisive_weights, chunksize)` | `RangeSelector2D c.matrix(const MatrixOptions&)` with fields `field`, `balance`, `sparse`, `as_pixels`, `join`, `ignore_index`, `divisive_weights`, `chunksize` |
| `c.open(mode)` (h5py handle) | not provided: coolercpp keeps HDF5 out of its public headers |

`balance` is a `Balance`: `true`/`false`, or a column name (`.balance = "KR"`).
The option called `max_chunk` in older cooler releases is `chunksize` in 0.10.2
and in coolercpp.

A `Region` is built implicitly from a string (`"chr3:10,000,000-12M"`) or
explicitly as `Region("chr3", start, end)` with `std::optional` bounds, the C++
form of the Python tuple `("chr3", start, end)`.

## Selectors

| Python | C++ |
|---|---|
| `sel["A"]` (Series selector) | `sel["A"]` or `sel[Fields("A")]` |
| `sel[["A", "B"]]` | `sel[Fields({"A", "B"})]` |
| `sel[lo:hi]` | `sel[Slice{lo, hi}]` or `sel.slice(lo, hi)` |
| `sel[:]` | `sel.all()` |
| `sel[i]` | `sel[i]` (a one row Table) |
| `sel.fetch(region)` | `sel.fetch(region)` |
| `sel.columns`, `sel.dtypes`, `sel.keys()`, `len(sel)`, `sel.shape`, `key in sel` | `columns()`, `dtypes()`, `keys()`, `size()`, `shape()`, `contains(key)` |
| `msel[i0:i1, j0:j1]` | `msel(Slice{i0, i1}, Slice{j0, j1})` |
| `msel[i0:i1]`, `msel[i]` | `msel[Slice{i0, i1}]`, `msel[i]` |
| `msel[:]` | `msel.all()` |
| `msel.fetch(region, region2)` | `msel.fetch(region, region2)` |

A matrix selector returns a `MatrixResult` holding a `DenseMatrix`
(`is_dense()`, `dense()`), a `SparseMatrix` (`is_sparse()`, `sparse()`) or a
pixel `Table` (`is_pixels()`, `pixels()`), matching the numpy array,
`coo_matrix` or DataFrame cooler returns for the same options.

## Module functions

| Python | C++ |
|---|---|
| `cooler.annotate(pixels, bins, replace=False)` | `annotate(const Table&, const Table& or const RangeSelector1D&, bool replace = false)` |
| `cooler.create_cooler(cool_uri, bins, pixels, columns, dtypes, metadata, assembly, ordered, symmetric_upper, mode, mergebuf, delete_temp, temp_dir, max_merge, boundscheck, dupcheck, triucheck, ensure_sorted, h5opts, lock)` | `create_cooler(uri, bins, pixels, const CreateOptions&)` where `pixels` is a `Table` (DataFrame or dict) or `PixelChunks` (an iterable of chunks); `CreateOptions` has every keyword except `lock` |
| `cooler.fileops.list_coolers(filepath)` | `list_coolers(filepath)` |
| `cooler.util.parse_cooler_uri(s)` | `parse_cooler_uri(s)` |
| `cooler.util.parse_humanized(s)` | `parse_humanized(s)` |
| `cooler.util.parse_region_string(s)` | `parse_region_string(s)` returning `GenomicRange` |
| `cooler.util.parse_region(reg, chromsizes)` | `parse_region(const Region&, const ChromSizes* = nullptr)` returning `RegionTuple` |
| `cooler.util.get_binsize(bins)` | `get_binsize(const Table&)` returning a `json::Value` number with the dtype of `end - start`, or null |
| `cooler.util.get_chromsizes(bins)` | `get_chromsizes(const Table&)` |
| `cooler.util.natsorted(items)` | `natsorted(items)` |

`CreateOptions` extensions that cooler does not have: `creation_date` (a fixed
value for reproducible files) and `generated_by` (defaults to
`"cooler-0.10.2"`). `h5opts` is an ordered list of `(key, H5OptValue)` pairs;
a tuple such as `(None,)` is a vector in which `-1` stands for `None`.

## The table type

`Column` holds exactly one of: `std::vector<Bool8>`, `std::vector<int8_t>` ...
`std::vector<uint64_t>`, `std::vector<float>`, `std::vector<double>`,
`std::vector<std::string>` (a pandas object column) or `CategoricalData`
(int32 codes plus shared labels). `values<T>()` hands out the typed buffer
without copying, `as<T>()` and `astype(DType)` convert with numpy's C casts,
`result_type` implements numpy 1.26 promotion. `Table` keeps column order and
an index, and holds any number of extra columns (`weight`, `KR`, `VC`, custom
value columns) next to the standard ones. pandas picks int8 or int16 codes for
small categoricals; coolercpp always stores int32 codes with the same values.
