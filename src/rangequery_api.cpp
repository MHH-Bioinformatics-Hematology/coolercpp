// The public range query engines and cooler.fileops.is_cooler (cooler 0.10.2,
// BSD-3-Clause).

#include "coolercpp/rangequery.hpp"

#include "coolercpp/errors.hpp"
#include "h5.hpp"
#include "rangequery.hpp"

namespace coolercpp {

namespace {

std::string join_path(const std::string& group, const std::string& name) {
    return group == "/" ? "/" + name : group + "/" + name;
}

}  // namespace

struct RangeQuery2D::Impl {
    std::unique_ptr<h5::File> file;
    std::unique_ptr<detail::CSRReader> reader;
    std::string field;
    detail::Bbox bbox;
    std::vector<detail::Task> tasks;
    bool reflect = false;
    bool return_index = false;
};

RangeQuery2D::RangeQuery2D(const Cooler& clr, Kind kind, std::string field,
                           std::array<std::int64_t, 4> bbox, std::int64_t chunksize,
                           bool return_index) {
    auto impl = std::make_shared<Impl>();
    impl->file = std::make_unique<h5::File>(clr.filename(), h5::Mode::Read);
    const h5::Dataset offsets =
        impl->file->open_dataset(join_path(clr.root(), "indexes/bin1_offset"));
    impl->reader = std::make_unique<detail::CSRReader>(
        *impl->file, join_path(clr.root(), "pixels"),
        offsets.read<std::int64_t>(0, offsets.length()));
    impl->field = std::move(field);
    impl->bbox = detail::Bbox{bbox[0], bbox[1], bbox[2], bbox[3]};
    impl->reflect = kind == Kind::FillLower;
    impl->return_index = return_index;
    impl->tasks = impl->reflect ? detail::fill_lower_tasks(*impl->reader, impl->bbox, chunksize)
                                : detail::direct_tasks(*impl->reader, impl->bbox, chunksize);
    impl_ = std::move(impl);
}

std::size_t RangeQuery2D::n_chunks() const noexcept { return impl_->tasks.size(); }

Table RangeQuery2D::get_chunk(std::size_t i) const {
    if (i >= impl_->tasks.size()) {
        throw IndexError("");
    }
    return detail::frame_from_dict(
        detail::run_tasks(*impl_->reader, impl_->field, impl_->bbox, {impl_->tasks[i]},
                          impl_->reflect, impl_->return_index),
        impl_->field);
}

Table RangeQuery2D::to_frame() const {
    return detail::frame_from_dict(
        detail::run_tasks(*impl_->reader, impl_->field, impl_->bbox, impl_->tasks,
                          impl_->reflect, impl_->return_index),
        impl_->field);
}

SparseMatrix RangeQuery2D::to_sparse_matrix() const {
    return detail::sparse_from_dict(
        detail::run_tasks(*impl_->reader, impl_->field, impl_->bbox, impl_->tasks,
                          impl_->reflect, impl_->return_index),
        impl_->bbox);
}

DenseMatrix RangeQuery2D::to_array() const {
    return detail::dense_from_dict(
        detail::run_tasks(*impl_->reader, impl_->field, impl_->bbox, impl_->tasks,
                          impl_->reflect, impl_->return_index),
        impl_->bbox);
}

bool is_cooler(const std::string& uri) {
    const auto [path, group] = parse_cooler_uri(uri);
    if (!h5::is_hdf5(path)) {
        return false;
    }
    const h5::File file(path, h5::Mode::Read);
    if (!file.exists(group)) {
        // h5py f[path]: HDF5 names the missing object when only the last
        // component is absent.
        const std::size_t slash = group.find_last_of('/');
        const std::string parent = slash == 0 ? "/" : group.substr(0, slash);
        if (file.exists(parent) && slash + 1 < group.size()) {
            throw KeyError("Unable to synchronously open object (object '" +
                           group.substr(slash + 1) + "' doesn't exist)");
        }
        throw KeyError("Unable to synchronously open object (component not found)");
    }
    for (const auto& [key, value] : file.attributes(group)) {
        if (key == "format") {
            if (!value.is_string() || value.as_string() != "HDF5::Cooler") {
                return false;
            }
            bool complete = file.is_group(group);
            for (const char* name : {"chroms", "bins", "pixels", "indexes"}) {
                complete = complete && file.exists(join_path(group, name));
            }
            if (!complete) {
                warn("Cooler path " + group + " appears to be corrupt");
            }
            return true;
        }
    }
    return false;
}

}  // namespace coolercpp
