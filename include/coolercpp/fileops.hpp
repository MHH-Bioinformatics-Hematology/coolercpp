// cooler.fileops: the predicates that recognise cool, mcool and scool files,
// the group listings and the copy, move and link operations
// (cooler/fileops.py of cooler 0.10.2).

#ifndef COOLERCPP_FILEOPS_HPP
#define COOLERCPP_FILEOPS_HPP

#include <optional>
#include <string>
#include <vector>

namespace coolercpp {

// cooler.fileops.is_cooler: whether a URI names a cooler group (its "format"
// attribute is "HDF5::Cooler"). False for a file that is not HDF5; a group
// path that does not exist raises KeyError, as h5py's lookup does.
[[nodiscard]] bool is_cooler(const std::string& uri);

// cooler.fileops.is_multires_file: whether a file holds a multi-resolution
// cooler, either under /resolutions (format "HDF5::MCOOL") or in the legacy
// layout whose first level is the group "0". min_version gates only the
// legacy layout.
[[nodiscard]] bool is_multires_file(const std::string& filepath, int min_version = 1);

// cooler.fileops.is_scool_file: whether a file holds a single-cell cooler
// (format "HDF5::SCOOL", with chroms, bins and a non-empty cells group whose
// every member is a cooler). A file that is not HDF5 raises OSError.
[[nodiscard]] bool is_scool_file(const std::string& filepath);

// cooler.fileops.list_coolers: the group paths of every cooler in a file, in
// natural sort order.
[[nodiscard]] std::vector<std::string> list_coolers(const std::string& filepath);

// cooler.fileops.list_scool_cells: the group paths of the cells of a scool
// file, in natural sort order. A file that is not a scool raises OSError.
[[nodiscard]] std::vector<std::string> list_scool_cells(const std::string& filepath);

// cooler.fileops.ls: every group and dataset under the URI's group, in the
// order of a depth first walk, each group's members in h5py's iteration
// order.
[[nodiscard]] std::vector<std::string> ls(const std::string& uri);

// cooler.fileops.cp: copies a group or dataset within one file or into
// another file.
void cp(const std::string& src_uri, const std::string& dst_uri, bool overwrite = false);

// cooler.fileops.mv: renames a group or dataset within one file.
void mv(const std::string& src_uri, const std::string& dst_uri, bool overwrite = false);

// cooler.fileops.ln: a hard link within one file, or with soft = true a soft
// link within one file and an external link between two files.
void ln(const std::string& src_uri, const std::string& dst_uri, bool soft = false,
        bool overwrite = false);

// cooler.fileops.pprint_data_tree: the object tree under the URI's group,
// rendered with the box drawing characters asciitree's LeftAligned renderer
// uses, each dataset followed by its shape and numpy dtype. level limits the
// depth; nullopt is Python's None.
[[nodiscard]] std::string pprint_data_tree(const std::string& uri,
                                           std::optional<int> level = std::nullopt);

}  // namespace coolercpp

#endif  // COOLERCPP_FILEOPS_HPP
