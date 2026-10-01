# API reference

The public headers, which carry the comments naming the cooler function each
part reproduces. The guide pages explain how they fit together.

## coolercpp/coolercpp.hpp

The umbrella header.

```cpp
--8<-- "include/coolercpp/coolercpp.hpp"
```

## coolercpp/api.hpp

The `Cooler` object, its selectors and the matrix results they return.

```cpp
--8<-- "include/coolercpp/api.hpp"
```

## coolercpp/create.hpp

`create_cooler`, `create_scool` and their options.

```cpp
--8<-- "include/coolercpp/create.hpp"
```

## coolercpp/fileops.hpp

What a file holds, and copying, moving and linking groups.

```cpp
--8<-- "include/coolercpp/fileops.hpp"
```

## coolercpp/table.hpp

`Table` and `Column`, the counterparts of the pandas DataFrame and Series.

```cpp
--8<-- "include/coolercpp/table.hpp"
```

## coolercpp/dtype.hpp

The element types, `Bool8` and the `dtype_of` traits.

```cpp
--8<-- "include/coolercpp/dtype.hpp"
```

## coolercpp/region.hpp

URIs, regions and chromosome sizes.

```cpp
--8<-- "include/coolercpp/region.hpp"
```

## coolercpp/rangequery.hpp

Range queries and the typed dataset reader.

```cpp
--8<-- "include/coolercpp/rangequery.hpp"
```

## coolercpp/util.hpp

Chromosome sizes, bin tables, natural sort and the genome segmentation helpers.

```cpp
--8<-- "include/coolercpp/util.hpp"
```

## coolercpp/json.hpp

The JSON values the file attributes carry.

```cpp
--8<-- "include/coolercpp/json.hpp"
```

## coolercpp/fastio.hpp

Threaded chunk reading and writing, and in-place editing.

```cpp
--8<-- "include/coolercpp/fastio.hpp"
```

## coolercpp/parallel.hpp

The thread pool the library runs its chunk work on.

```cpp
--8<-- "include/coolercpp/parallel.hpp"
```

## coolercpp/errors.hpp

The exception types, named after the Python exceptions they stand for.

```cpp
--8<-- "include/coolercpp/errors.hpp"
```
