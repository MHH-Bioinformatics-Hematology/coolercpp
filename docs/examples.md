# Examples

The programs below are the files under `examples/` in the repository. They are
complete: every table, option and variable they use is defined in the code, and
a top-level build compiles them, so an example that stops compiling fails the
build.

They are built with the library and left in the build tree:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/examples/create_cool out.cool
```

A program outside this repository builds against an installed coolercpp as
[the installation page](install.md) describes; the sources below need no change.

## Quickstart

`quickstart` opens a file and reads one region out of it.

```cpp title="examples/quickstart.cpp"
--8<-- "examples/quickstart.cpp"
```

```
$ quickstart hg19.GM12878-MboI.matrix.2000kb.cool chr1:10M-60M
chr1:10M-60M covers 25 bins and holds 335 of 25 x 25 pixels
```

## Reading a file

`read_cool` opens a cool or mcool URI, prints what the file holds, and fetches a
bin table and a matrix region of the first chromosome.

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

## The selectors

`selectors` goes through the four selectors of a `Cooler` object and the forms
each one returns.

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

## Writing a file

`create_cool` builds the bin table, with the categorical chromosome column the
format stores, and the pixel table, then writes a file and reads it back.

```cpp title="examples/create_cool.cpp"
--8<-- "examples/create_cool.cpp"
```

```
$ create_cool out.cool
wrote out.cool: 2 chromosomes, 5 bins, 5 pixels, sum 18
```

## Single-cell files and file level operations

`scool_fileops` writes a single-cell file of three cells sharing one bin table,
lists what the file holds, reads one cell as a cooler, and copies, renames and
links groups.

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
