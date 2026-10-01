# File level operations

The `fileops` header carries what `cooler.fileops` provides: what a file holds,
and copying, moving and linking groups.

`is_cooler`, `is_multires_file` and `is_scool_file` say what a URI or a file
holds. `list_coolers` and `list_scool_cells` list the cooler groups of a file in
the natural sort order cooler lists them in. `cp`, `mv` and `ln` copy, rename and
link a group, and `pprint_data_tree` prints the tree of a file as cooler prints
it.

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
