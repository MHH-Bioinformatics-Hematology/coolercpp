# File level operations

The `fileops` header carries what `cooler.fileops` provides: what a file holds,
and copying, moving and linking groups.

```cpp
#include <coolercpp/fileops.hpp>

coolercpp::is_cooler("matrix.mcool::/resolutions/10000");
coolercpp::is_multires_file("matrix.mcool");
coolercpp::is_scool_file("cells.scool");
coolercpp::list_coolers("matrix.mcool");      // natural sort order, as cooler lists them
coolercpp::list_scool_cells("cells.scool");
coolercpp::cp("a.cool::/", "b.cool::/copy");
coolercpp::mv("b.cool::/copy", "b.cool::/moved");
coolercpp::ln("b.cool::/moved", "b.cool::/link", true);   // soft link
```

`pprint_data_tree` prints the tree of a file as cooler prints it.

The whole of `examples/scool_fileops.cpp`, which writes a single-cell file of
three cells, lists what it holds and then copies, renames and links a cell:

```cpp title="examples/scool_fileops.cpp"
--8<-- "examples/scool_fileops.cpp"
```

```
$ scool_fileops cells.scool
cells.scool: scool true, mcool false
  /cells/cell_1: 3 pixels, sum 15
  /cells/cell_2: 3 pixels, sum 12
  /cells/cell_3: 2 pixels, sum 10
cell_1 holds 5 bins over 2 chromosomes at 10 bp
after cp, mv and ln: /cells/cell_1 /cells/cell_2 /cells/cell_3 /cells/cell_copy /cells/cell_link
5 coolers, root holds 84 objects
cell_1
 ├── bins
 ├── chroms
 ├── indexes
 └── pixels
```
