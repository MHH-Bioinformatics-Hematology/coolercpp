// Port of cooler/fileops.py of cooler 0.10.2 (BSD-3-Clause): _is_cooler,
// is_cooler, is_multires_file, is_scool_file, list_coolers, list_scool_cells,
// ls, visititems, _copy with cp, mv and ln, and the data tree rendering of
// pprint_data_tree (TreeNode plus asciitree's LeftAligned and BoxStyle).

#include "coolercpp/fileops.hpp"

#include <sys/stat.h>

#include <algorithm>
#include <functional>
#include <utility>

#include "coolercpp/errors.hpp"
#include "coolercpp/region.hpp"
#include "coolercpp/util.hpp"
#include "coolercpp/version.hpp"
#include "h5.hpp"

namespace coolercpp {

namespace {

std::string join_path(const std::string& group, const std::string& name) {
    return group == "/" ? "/" + name : group + "/" + name;
}

// os.path.isfile.
bool is_regular_file(const std::string& path) {
    struct stat info {};
    return ::stat(path.c_str(), &info) == 0 && S_ISREG(info.st_mode);
}

// cooler.fileops._is_cooler. The warning for an incomplete group is cooler's.
bool is_cooler_group(const h5::File& file, const std::string& path) {
    for (const auto& [key, value] : file.attributes(path)) {
        if (key != "format") {
            continue;
        }
        if (!value.is_string() || value.as_string() != kMagic) {
            return false;
        }
        bool complete = file.is_group(path);
        for (const char* name : {"chroms", "bins", "pixels", "indexes"}) {
            complete = complete && file.exists(join_path(path, name));
        }
        if (!complete) {
            warn("Cooler path " + path + " appears to be corrupt");
        }
        return true;
    }
    return false;
}

// cooler.fileops.visititems: every descendant of a group, each group's members
// in h5py's iteration order, depth first. Unlike h5py's own visititems this
// one does not stop at an object it has already seen, so the datasets a scool
// cell hard links to the root bin table are listed under every cell.
void visititems(const h5::File& file, const std::string& group,
                const std::function<void(const std::string&)>& visit, int depth,
                std::optional<int> level) {
    if (level.has_value() && depth >= *level) {
        return;
    }
    for (const std::string& name : file.keys(group)) {
        const std::string child = join_path(group, name);
        visit(child);
        if (file.is_group(child)) {
            visititems(file, child, visit, depth + 1, level);
        }
    }
}

// What HDF5 says about a path that does not resolve: the missing last
// component by name, or that a component along the way was not found.
std::string missing_object(const h5::File& file, const std::string& group) {
    const std::size_t slash = group.find_last_of('/');
    const std::string parent = slash == 0 ? "/" : group.substr(0, slash);
    if (file.exists(parent) && slash + 1 < group.size()) {
        return "object '" + group.substr(slash + 1) + "' doesn't exist";
    }
    return "component not found";
}

// The group path h5py's own f[group_path] resolves, with the KeyError h5py
// raises for a missing one.
void require_object(const h5::File& file, const std::string& group) {
    if (!file.exists(group)) {
        throw KeyError("Unable to synchronously open object (" + missing_object(file, group) + ")");
    }
}

// H5Ocopy reports a missing source and an occupied destination itself, and
// h5py turns both into a RuntimeError.
void check_copy(const h5::File& src, const std::string& src_group, const h5::File& dst,
                const std::string& dst_group) {
    if (!src.exists(src_group)) {
        throw RuntimeError("Unable to synchronously copy object (" +
                           missing_object(src, src_group) + ")");
    }
    if (dst.exists(dst_group)) {
        throw RuntimeError(
            "Unable to synchronously copy object (destination object already exists)");
    }
}

// ---------------------------------------------------------------------------
// The data tree of pprint_data_tree

// repr() of the shape tuple h5py reports for a dataset.
std::string shape_repr(const std::vector<std::size_t>& shape) {
    std::string out = "(";
    for (std::size_t k = 0; k < shape.size(); ++k) {
        out += std::to_string(shape[k]);
        if (shape.size() == 1 || k + 1 < shape.size()) {
            out += shape.size() == 1 ? "," : ", ";
        }
    }
    return out + ")";
}

// TreeNode.get_text.
std::string node_text(const h5::File& file, const std::string& path) {
    const std::size_t slash = path.find_last_of('/');
    std::string name = slash == std::string::npos ? path : path.substr(slash + 1);
    if (name.empty()) {
        name = "/";
    }
    if (file.is_dataset(path)) {
        const h5::Dataset dataset = file.open_dataset(path);
        name += " " + shape_repr(dataset.shape()) + " " +
                h5::numpy_dtype_name(dataset.type());
    }
    return name;
}

// asciitree.LeftAligned.render with BoxStyle(gfx=BOX_LIGHT, horiz_len=2,
// label_space=1, indent=1): the lines of one subtree, its own label first.
std::vector<std::string> render_tree(const h5::File& file, const std::string& path, int depth,
                                     std::optional<int> level) {
    std::vector<std::string> lines{node_text(file, path)};
    std::vector<std::string> children;
    if (file.is_group(path) && (!level.has_value() || depth < *level)) {
        for (const std::string& name : file.keys(path)) {
            children.push_back(join_path(path, name));
        }
    }
    for (std::size_t k = 0; k < children.size(); ++k) {
        const bool last = k + 1 == children.size();
        std::vector<std::string> subtree = render_tree(file, children[k], depth + 1, level);
        lines.push_back(std::string(" ") + (last ? "└" : "├") + "── " +
                        subtree.front());
        for (std::size_t i = 1; i < subtree.size(); ++i) {
            lines.push_back(std::string(" ") + (last ? " " : "│") + "  " + subtree[i]);
        }
    }
    return lines;
}

// ---------------------------------------------------------------------------
// cooler.fileops._copy

void copy_impl(const std::string& src_uri, const std::string& dst_uri, bool overwrite, bool link,
               bool rename, bool soft_link) {
    const auto [src_path, src_group] = parse_cooler_uri(src_uri);
    const auto [dst_path, dst_group] = parse_cooler_uri(dst_uri);

    if (static_cast<int>(link) + static_cast<int>(rename) + static_cast<int>(soft_link) > 1) {
        throw ValueError("Must provide at most one of: \"link\", \"rename\", \"soft_link\"");
    }
    const h5::Mode dst_mode =
        (!is_regular_file(dst_path) || overwrite) ? h5::Mode::Truncate : h5::Mode::ReadWrite;
    const h5::Mode src_mode = src_path == dst_path ? h5::Mode::ReadWrite : h5::Mode::Read;

    h5::File src(src_path, src_mode);
    if (src_path == dst_path) {
        // Both URIs name the same file; cooler opens it twice and works
        // through the source handle only.
        h5::File dst(dst_path, dst_mode);
        if (link || rename) {
            // h5py resolves src[src_group] before linking it.
            require_object(src, src_group);
            src.hard_link(dst_group, src_group);
            if (rename) {
                src.remove(src_group);
            }
        } else if (soft_link) {
            src.soft_link(dst_group, src_group);
        } else {
            check_copy(src, src_group, dst, dst_group);
            src.copy(src_group, dst_group);
        }
        return;
    }

    h5::File dst(dst_path, dst_mode);
    if (link) {
        throw OSError("Can't hard link between two different files.");
    }
    if (soft_link) {
        dst.external_link(dst_group, src_path, src_group);
        return;
    }
    if (dst_group == "/") {
        require_object(src, src_group);
        for (const std::string& name : src.keys(src_group)) {
            src.copy(join_path(src_group, name), dst, name);
        }
        dst.copy_attributes(src, src_group, dst_group);
        return;
    }
    check_copy(src, src_group, dst, dst_group);
    src.copy(src_group, dst, dst_group);
}

}  // namespace

bool is_cooler(const std::string& uri) {
    const auto [path, group] = parse_cooler_uri(uri);
    if (!h5::is_hdf5(path)) {
        return false;
    }
    const h5::File file(path, h5::Mode::Read);
    require_object(file, group);
    return is_cooler_group(file, group);
}

bool is_multires_file(const std::string& filepath, int min_version) {
    if (!h5::is_hdf5(filepath)) {
        return false;
    }
    const h5::File file(filepath, h5::Mode::Read);
    std::string format;
    for (const auto& [key, value] : file.attributes("/")) {
        if (key == "format" && value.is_string()) {
            format = value.as_string();
        }
    }
    const std::vector<std::string> root = file.keys("/");
    const auto has = [&root](const std::string& name) {
        return std::find(root.begin(), root.end(), name) != root.end();
    };
    const std::vector<std::string> levels = has("resolutions") ? file.keys("/resolutions")
                                                               : std::vector<std::string>{};
    if (has("resolutions") && !levels.empty()) {
        // Only the first resolution is examined, and min_version does not gate
        // this layout: an mcool is multi-res whichever version is asked for.
        return format == kMagicMcool && is_cooler_group(file, "/resolutions/" + levels.front());
    }
    // The legacy layout: one group per zoom level, the finest called "0". An
    // empty /resolutions group falls through to it, as cooler's elif does.
    return has("0") && is_cooler_group(file, "/0") && min_version < 2;
}

bool is_scool_file(const std::string& filepath) {
    if (!h5::is_hdf5(filepath)) {
        // cooler raises here instead of returning False, unlike the other two
        // predicates.
        throw OSError("'" + filepath + "' is not an HDF5 file.");
    }
    const h5::File file(filepath, h5::Mode::Read);
    bool is_scool = false;
    for (const auto& [key, value] : file.attributes("/")) {
        if (key == "format" && value.is_string() && value.as_string() == kMagicScool) {
            is_scool = true;
        }
    }
    if (!is_scool) {
        return false;
    }
    const std::vector<std::string> root = file.keys("/");
    for (const char* name : {"chroms", "bins", "cells"}) {
        if (std::find(root.begin(), root.end(), name) == root.end()) {
            warn("Scool file appears to be corrupt");
            return false;
        }
    }
    const std::vector<std::string> cells = file.keys("/cells");
    if (cells.empty()) {
        return false;
    }
    for (const std::string& cell : cells) {
        if (!is_cooler_group(file, "/cells/" + cell)) {
            return false;
        }
    }
    return true;
}

std::vector<std::string> list_coolers(const std::string& filepath) {
    if (!h5::is_hdf5(filepath)) {
        throw OSError("'" + filepath + "' is not an HDF5 file.");
    }
    const h5::File file(filepath, h5::Mode::Read);
    std::vector<std::string> listing;
    const auto check = [&](const std::string& path) {
        if (is_cooler_group(file, path)) {
            listing.push_back(path);
        }
    };
    check("/");
    visititems(file, "/", check, 0, std::nullopt);
    return natsorted(std::move(listing));
}

std::vector<std::string> list_scool_cells(const std::string& filepath) {
    if (!is_scool_file(filepath)) {
        throw OSError("'" + filepath + "' is not a scool file.");
    }
    const h5::File file(filepath, h5::Mode::Read);
    std::vector<std::string> listing;
    const auto check = [&](const std::string& path) {
        if (is_cooler_group(file, path)) {
            listing.push_back(path);
        }
    };
    check("/");
    visititems(file, "/", check, 0, std::nullopt);
    // The root of a scool file is not a cooler, so this only matters for a
    // file that carries both magic strings.
    const auto root = std::find(listing.begin(), listing.end(), "/");
    if (root != listing.end()) {
        listing.erase(root);
    }
    return natsorted(std::move(listing));
}

std::vector<std::string> ls(const std::string& uri) {
    const auto [path, group] = parse_cooler_uri(uri);
    if (!h5::is_hdf5(path)) {
        throw OSError("'" + path + "' is not an HDF5 file.");
    }
    const h5::File file(path, h5::Mode::Read);
    require_object(file, group);
    std::vector<std::string> listing{group};
    visititems(file, group, [&listing](const std::string& child) { listing.push_back(child); }, 0,
               std::nullopt);
    return listing;
}

void cp(const std::string& src_uri, const std::string& dst_uri, bool overwrite) {
    copy_impl(src_uri, dst_uri, overwrite, false, false, false);
}

void mv(const std::string& src_uri, const std::string& dst_uri, bool overwrite) {
    copy_impl(src_uri, dst_uri, overwrite, false, true, false);
}

void ln(const std::string& src_uri, const std::string& dst_uri, bool soft, bool overwrite) {
    copy_impl(src_uri, dst_uri, overwrite, !soft, false, soft);
}

std::string pprint_data_tree(const std::string& uri, std::optional<int> level) {
    const auto [path, group] = parse_cooler_uri(uri);
    const h5::File file(path, h5::Mode::Read);
    require_object(file, group);
    const std::vector<std::string> lines = render_tree(file, group, 0, level);
    std::string out;
    for (std::size_t k = 0; k < lines.size(); ++k) {
        if (k != 0) {
            out += "\n";
        }
        out += lines[k];
    }
    return out;
}

}  // namespace coolercpp
