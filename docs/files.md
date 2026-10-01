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
