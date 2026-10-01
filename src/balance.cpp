// Port of cooler/balance.py of cooler 0.10.2 (BSD-3-Clause): balance_cooler
// and the split-apply-combine pipeline it runs over the pixel table.
//
// The arithmetic follows the Python operation by operation, because the weights
// are compared against cooler's: the marginal of a chunk is two separate
// bincount accumulations added together, the chunk marginals are folded into
// the total in file order, and the mean, variance and median of the nonzero
// marginals go through numpy's pairwise summation (numpy_compat.hpp).

#include "coolercpp/balance.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "coolercpp/errors.hpp"
#include "coolercpp/fastio.hpp"
#include "coolercpp/util.hpp"
#include "h5.hpp"
#include "numpy_compat.hpp"

namespace coolercpp {

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

std::string join_path(const std::string& group, const std::string& name) {
    return group == "/" ? "/" + name : group + "/" + name;
}

using Span = std::pair<std::int64_t, std::int64_t>;

// int(clr.info[key]).
std::int64_t info_int(const json::Value& info, const std::string& key) {
    const json::Value* value = info.find(key);
    if (value == nullptr) {
        throw KeyError(key);
    }
    if (value->type() == json::Type::Double) {
        return static_cast<std::int64_t>(value->as_double());
    }
    if (!value->is_integer()) {
        throw TypeError("int() argument must be a string, a bytes-like object or a real number");
    }
    return value->as_int();
}

// The pre-marginalization data transformations of cooler.balance, as flags: the
// pipeline is always a subset of _binarize, _zero_trans, _zero_diags and
// _zero_cis, applied in that order, so one pass over a chunk does all of them.
struct Filters {
    bool binarize = false;
    bool zero_trans = false;
    bool zero_diags = false;
    std::int64_t n_diags = 0;
    bool zero_cis = false;
};

// What cooler.parallel.chunkgetter hands to the pipeline: one span of the pixel
// table, plus the chromosome id of every bin. cooler reopens the file, rereads
// every column of the bin table and rereads the span for every iteration;
// coolercpp keeps the file open and the chromosome ids in memory, and rereads
// only the pixels, which changes nothing that is computed.
class PixelChunks {
  public:
    PixelChunks(const std::string& path, const std::string& root, int threads)
        : pool_(threads),
          file_(path, fastio::Mode::Read),
          bin1_(file_, join_path(root, "pixels/bin1_id")),
          bin2_(file_, join_path(root, "pixels/bin2_id")),
          count_(file_, join_path(root, "pixels/count")) {}

    // grp["pixels"][lo:hi]: a range past the end of the table is clipped, as
    // Python slicing clips it, and one entirely past it is empty.
    void load(std::int64_t lo, std::int64_t hi) {
        const auto length = static_cast<std::int64_t>(bin1_.size());
        lo = std::clamp<std::int64_t>(lo, 0, length);
        hi = std::clamp<std::int64_t>(hi, lo, length);
        rows_ = static_cast<std::size_t>(hi - lo);
        bin1_buffer_.resize(rows_);
        bin2_buffer_.resize(rows_);
        count_buffer_.resize(rows_);
        const auto from = static_cast<std::size_t>(lo);
        const auto to = static_cast<std::size_t>(hi);
        bin1_.read(pool_, from, to, bin1_buffer_.data());
        bin2_.read(pool_, from, to, bin2_buffer_.data());
        count_.read(pool_, from, to, count_buffer_.data());
    }

    [[nodiscard]] std::size_t rows() const noexcept { return rows_; }
    [[nodiscard]] const std::int64_t* bin1() const noexcept { return bin1_buffer_.data(); }
    [[nodiscard]] const std::int64_t* bin2() const noexcept { return bin2_buffer_.data(); }
    // _init: np.copy(chunk["pixels"]["count"]). Counts are read as float64;
    // every integer and float count column cooler writes converts exactly, and
    // numpy widens the counts to float64 in the same place.
    [[nodiscard]] const double* count() const noexcept { return count_buffer_.data(); }

  private:
    ThreadPool pool_;
    fastio::File file_;
    fastio::ColumnReader bin1_;
    fastio::ColumnReader bin2_;
    fastio::ColumnReader count_;
    std::size_t rows_ = 0;
    std::vector<std::int64_t> bin1_buffer_;
    std::vector<std::int64_t> bin2_buffer_;
    std::vector<double> count_buffer_;
};

// split(clr, spans).prepare(_init).pipe(filters).pipe(_timesouterproduct, vec)
// .pipe(_marginalize).reduce(add, np.zeros(n_bins)), with vec left out when it
// is null. _marginalize is two np.bincount calls that are added together, so
// the row and the column sums are accumulated apart and only then combined.
std::vector<double> marginalize(PixelChunks& pixels, const std::vector<Span>& spans,
                                const std::vector<std::int64_t>& chrom_ids,
                                const Filters& filters, const double* vec) {
    const std::size_t n = chrom_ids.size();
    std::vector<double> total(n, 0.0);
    std::vector<double> row(n);
    std::vector<double> col(n);
    for (const auto& [lo, hi] : spans) {
        pixels.load(lo, hi);
        std::fill(row.begin(), row.end(), 0.0);
        std::fill(col.begin(), col.end(), 0.0);
        const std::int64_t* bin1 = pixels.bin1();
        const std::int64_t* bin2 = pixels.bin2();
        const double* count = pixels.count();
        for (std::size_t i = 0; i < pixels.rows(); ++i) {
            const std::int64_t u = bin1[i];
            const std::int64_t v = bin2[i];
            // A guard against a corrupt pixel table, where Python would fail
            // inside np.bincount or while indexing the chrom column.
            if (u < 0 || v < 0 || u >= static_cast<std::int64_t>(n) ||
                v >= static_cast<std::int64_t>(n)) {
                throw ValueError("a pixel bin id lies outside the bin table");
            }
            double data = count[i];
            if (filters.binarize && data != 0.0) {
                data = 1.0;
            }
            if (filters.zero_trans && chrom_ids[static_cast<std::size_t>(u)] !=
                                          chrom_ids[static_cast<std::size_t>(v)]) {
                data = 0.0;
            }
            if (filters.zero_diags && (u > v ? u - v : v - u) < filters.n_diags) {
                data = 0.0;
            }
            if (filters.zero_cis && chrom_ids[static_cast<std::size_t>(u)] ==
                                        chrom_ids[static_cast<std::size_t>(v)]) {
                data = 0.0;
            }
            if (vec != nullptr) {
                // _timesouterproduct: vec[bin1] * vec[bin2] * data, left to
                // right, so a zero multiplied by a NaN weight stays NaN.
                data = (vec[static_cast<std::size_t>(u)] * vec[static_cast<std::size_t>(v)]) * data;
            }
            row[static_cast<std::size_t>(u)] += data;
            col[static_cast<std::size_t>(v)] += data;
        }
        for (std::size_t k = 0; k < n; ++k) {
            total[k] += row[k] + col[k];
        }
    }
    return total;
}

// marg[marg != 0], which keeps NaNs.
std::vector<double> nonzero(const std::span<const double> values) {
    std::vector<double> out;
    for (const double value : values) {
        if (value != 0.0) {
            out.push_back(value);
        }
    }
    return out;
}

const char* kUnboundNzmarg =
    "cannot access local variable 'nzmarg' where it is not associated with a value";

struct WholeStats {
    double scale = 1.0;
    double var = 0.0;
};

// cooler.balance._balance_genomewide, and _balance_transonly when cweights is
// given: the only differences are that the trans-only pipeline also zeroes the
// cis pixels and that it scales the bias by the per-chromosome weights before
// every outer product.
WholeStats balance_whole(std::vector<double>& bias, PixelChunks& pixels,
                         const std::vector<Span>& spans,
                         const std::vector<std::int64_t>& chrom_ids, const Filters& filters,
                         const std::vector<double>* cweights, const BalanceOptions& options) {
    const std::size_t n = bias.size();
    std::vector<double> nzmarg;
    std::vector<double> scaled;
    bool assigned = false;
    bool exhausted = true;  // the for/else branch of the Python loop
    double var = 0.0;
    for (std::int64_t iteration = 0; iteration < options.max_iters; ++iteration) {
        const double* vec = bias.data();
        if (cweights != nullptr) {
            scaled.resize(n);
            for (std::size_t i = 0; i < n; ++i) {
                scaled[i] = bias[i] * (*cweights)[i];
            }
            vec = scaled.data();
        }
        std::vector<double> marg = marginalize(pixels, spans, chrom_ids, filters, vec);
        nzmarg = nonzero(marg);
        assigned = true;
        if (nzmarg.empty()) {
            std::fill(bias.begin(), bias.end(), kNaN);
            var = 0.0;
            exhausted = false;
            break;
        }
        const double centre = npy::mean(nzmarg);
        for (std::size_t i = 0; i < n; ++i) {
            marg[i] /= centre;
        }
        for (std::size_t i = 0; i < n; ++i) {
            if (marg[i] == 0.0) {
                marg[i] = 1.0;
            }
        }
        for (std::size_t i = 0; i < n; ++i) {
            bias[i] /= marg[i];
        }
        var = npy::var(nzmarg);
        if (var < options.tol) {
            exhausted = false;
            break;
        }
    }
    if (exhausted) {
        warn("Iteration limit reached without convergence.");
    }
    if (!assigned) {
        throw UnboundLocalError(kUnboundNzmarg);
    }
    const double scale = npy::mean(nzmarg);
    for (std::size_t i = 0; i < n; ++i) {
        if (bias[i] == 0.0) {
            bias[i] = kNaN;
        }
    }
    if (options.rescale_marginals) {
        const double root = std::sqrt(scale);
        for (std::size_t i = 0; i < n; ++i) {
            bias[i] /= root;
        }
    }
    return {scale, var};
}

struct ChromStats {
    std::vector<double> scales;
    std::vector<double> variances;
};

// cooler.balance._balance_cisonly. The marginal is computed over the whole bin
// range and then cut down to the chromosome, and `nzmarg` outlives the
// chromosome loop in Python, so a chromosome whose inner loop does not run
// reuses the previous chromosome's value.
ChromStats balance_cisonly(std::vector<double>& bias, PixelChunks& pixels,
                           const std::vector<std::int64_t>& chrom_ids,
                           const std::vector<std::string>& chroms,
                           const std::vector<std::int64_t>& chrom_offsets,
                           const std::vector<std::int64_t>& bin1_offsets, const Filters& filters,
                           std::int64_t chunksize, const BalanceOptions& options) {
    const std::size_t n_chroms = chroms.size();
    ChromStats out{std::vector<double>(n_chroms, 1.0), std::vector<double>(n_chroms, kNaN)};
    std::vector<double> nzmarg;
    bool assigned = false;
    for (std::size_t cid = 0; cid < n_chroms; ++cid) {
        const std::int64_t lo = chrom_offsets[cid];
        const std::int64_t hi = chrom_offsets[cid + 1];
        const std::vector<Span> spans = partition(bin1_offsets[static_cast<std::size_t>(lo)],
                                                  bin1_offsets[static_cast<std::size_t>(hi)],
                                                  chunksize);
        double var = kNaN;
        bool exhausted = true;
        for (std::int64_t iteration = 0; iteration < options.max_iters; ++iteration) {
            const std::vector<double> full =
                marginalize(pixels, spans, chrom_ids, filters, bias.data());
            std::vector<double> marg(full.begin() + lo, full.begin() + hi);
            nzmarg = nonzero(marg);
            assigned = true;
            if (nzmarg.empty()) {
                std::fill(bias.begin() + lo, bias.begin() + hi, kNaN);
                var = 0.0;
                exhausted = false;
                break;
            }
            const double centre = npy::mean(nzmarg);
            for (double& value : marg) {
                value /= centre;
            }
            for (double& value : marg) {
                if (value == 0.0) {
                    value = 1.0;
                }
            }
            for (std::size_t k = 0; k < marg.size(); ++k) {
                bias[static_cast<std::size_t>(lo) + k] /= marg[k];
            }
            var = npy::var(nzmarg);
            if (var < options.tol) {
                exhausted = false;
                break;
            }
        }
        if (exhausted) {
            warn("Iteration limit reached without convergence on " + chroms[cid] + ".");
        }
        if (!assigned) {
            throw UnboundLocalError(kUnboundNzmarg);
        }
        const double scale = npy::mean(nzmarg);
        for (std::int64_t k = lo; k < hi; ++k) {
            if (bias[static_cast<std::size_t>(k)] == 0.0) {
                bias[static_cast<std::size_t>(k)] = kNaN;
            }
        }
        out.scales[cid] = scale;
        out.variances[cid] = var;
        if (options.rescale_marginals) {
            const double root = std::sqrt(scale);
            for (std::int64_t k = lo; k < hi; ++k) {
                bias[static_cast<std::size_t>(k)] /= root;
            }
        }
    }
    return out;
}

// The chromosome id of every bin. cooler compares the two ends of a pixel
// through pandas Categoricals built from the enum dataset, which compares the
// stored codes; a bin table whose chrom column holds plain integers or strings
// is compared by value, so distinct labels become distinct ids here.
std::vector<std::int64_t> read_chrom_ids(const h5::File& file, const std::string& root,
                                        std::int64_t n_bins) {
    const h5::Dataset dataset = file.open_dataset(join_path(root, "bins/chrom"));
    const Column column = dataset.read_column(0, static_cast<std::size_t>(n_bins));
    if (column.dtype() == DType::String) {
        const std::vector<std::string>& labels = column.values<std::string>();
        std::map<std::string, std::int64_t> ids;
        std::vector<std::int64_t> out(labels.size());
        for (std::size_t i = 0; i < labels.size(); ++i) {
            out[i] = ids.emplace(labels[i], static_cast<std::int64_t>(ids.size())).first->second;
        }
        return out;
    }
    return column.as<std::int64_t>();
}

std::vector<std::int64_t> read_index(const h5::File& file, const std::string& path) {
    const h5::Dataset dataset = file.open_dataset(path);
    return dataset.read<std::int64_t>(0, dataset.length());
}

json::Value number_array(const std::vector<double>& values) {
    json::Value out = json::Value::array({});
    for (const double value : values) {
        out.push_back(json::Value(value));
    }
    return out;
}

json::Value bool_array(const std::vector<double>& values, double tol) {
    json::Value out = json::Value::array({});
    for (const double value : values) {
        out.push_back(json::Value(value < tol));
    }
    return out;
}

}  // namespace

BalanceResult balance_cooler(const Cooler& clr, const BalanceOptions& options) {
    const json::Value info = clr.info();
    const std::int64_t nnz = info_int(info, "nnz");
    const std::int64_t n_bins = info_int(info, "nbins");

    // Divide the number of elements into non-overlapping chunks. The edges run
    // one chunk past nnz, so the last span can reach beyond the pixel table and
    // an nnz that is a multiple of the chunk size gains an empty span.
    const std::int64_t chunksize = options.chunksize.value_or(nnz);
    std::vector<Span> spans;
    if (!options.chunksize.has_value()) {
        spans.emplace_back(0, nnz);
    } else {
        if (chunksize <= 0) {
            throw ValueError("Maximum allowed size exceeded");
        }
        std::vector<std::int64_t> edges;
        for (std::int64_t edge = 0; edge < nnz + chunksize; edge += chunksize) {
            edges.push_back(edge);
        }
        for (std::size_t i = 0; i + 1 < edges.size(); ++i) {
            spans.emplace_back(edges[i], edges[i + 1]);
        }
    }

    // List of pre-marginalization data transformations.
    Filters base;
    base.zero_trans = options.cis_only;
    if (options.ignore_diags.enabled()) {
        base.zero_diags = true;
        base.n_diags = options.ignore_diags.n_diags();
    }

    // Initialize the bias weights. cooler keeps the caller's x0 array and
    // returns it; coolercpp copies it.
    std::vector<double> bias;
    if (options.x0.has_value()) {
        bias = *options.x0;
        if (static_cast<std::int64_t>(bias.size()) != n_bins) {
            throw ValueError("x0 has " + std::to_string(bias.size()) + " elements, expected " +
                             std::to_string(n_bins));
        }
        for (double& value : bias) {
            if (std::isnan(value)) {
                value = 0.0;
            }
        }
    } else {
        bias.assign(static_cast<std::size_t>(n_bins), 1.0);
    }

    const auto width = static_cast<std::size_t>(n_bins);
    double scale = 1.0;
    double var = 0.0;
    ChromStats per_chrom;
    {
        const h5::File file(clr.store(), h5::Mode::Read);
        const std::vector<std::int64_t> chrom_ids = read_chrom_ids(file, clr.root(), n_bins);
        const std::vector<std::int64_t> offsets =
            read_index(file, join_path(clr.root(), "indexes/chrom_offset"));
        PixelChunks pixels(clr.store(), clr.root(), options.threads);

        // Drop bins with too few nonzeros from bias.
        if (options.min_nnz > 0) {
            Filters filters = base;
            filters.binarize = true;
            const std::vector<double> marg_nnz =
                marginalize(pixels, spans, chrom_ids, filters, nullptr);
            for (std::size_t i = 0; i < width; ++i) {
                if (marg_nnz[i] < static_cast<double>(options.min_nnz)) {
                    bias[i] = 0.0;
                }
            }
        }

        std::vector<double> marg = marginalize(pixels, spans, chrom_ids, base, nullptr);

        // Drop bins with too few total counts from bias.
        if (options.min_count != 0) {
            for (std::size_t i = 0; i < width; ++i) {
                if (marg[i] < static_cast<double>(options.min_count)) {
                    bias[i] = 0.0;
                }
            }
        }

        // MAD-max filter on the marginals. A chromosome without a single
        // positive marginal is divided by the median of an empty array, which
        // turns its marginals into NaN; NaN compares false against the cutoff,
        // so those bins keep their weight instead of being dropped.
        if (options.mad_max > 0) {
            for (std::size_t cid = 0; cid + 1 < offsets.size(); ++cid) {
                const auto lo = static_cast<std::size_t>(offsets[cid]);
                const auto hi = static_cast<std::size_t>(offsets[cid + 1]);
                std::vector<double> positive;
                for (std::size_t i = lo; i < hi; ++i) {
                    if (marg[i] > 0.0) {
                        positive.push_back(marg[i]);
                    }
                }
                const double centre = npy::median(positive);
                for (std::size_t i = lo; i < hi; ++i) {
                    marg[i] /= centre;
                }
            }
            std::vector<double> logs;
            for (std::size_t i = 0; i < width; ++i) {
                if (marg[i] > 0.0) {
                    logs.push_back(std::log(marg[i]));
                }
            }
            const double centre = npy::median(logs);
            const double deviation = mad(logs);
            const double cutoff = std::exp(centre - static_cast<double>(options.mad_max) * deviation);
            for (std::size_t i = 0; i < width; ++i) {
                if (marg[i] < cutoff) {
                    bias[i] = 0.0;
                }
            }
        }

        // Filter out pre-determined bad bins.
        if (options.blacklist.has_value()) {
            for (const std::int64_t given : *options.blacklist) {
                const std::int64_t index = given < 0 ? given + n_bins : given;
                if (index < 0 || index >= n_bins) {
                    throw IndexError("index " + std::to_string(given) +
                                     " is out of bounds for axis 0 with size " +
                                     std::to_string(n_bins));
                }
                bias[static_cast<std::size_t>(index)] = 0.0;
            }
        }

        // Do balancing.
        if (options.cis_only) {
            const std::vector<std::int64_t> bin1_offsets =
                read_index(file, join_path(clr.root(), "indexes/bin1_offset"));
            per_chrom = balance_cisonly(bias, pixels, chrom_ids, clr.chromnames(), offsets,
                                        bin1_offsets, base, chunksize, options);
        } else if (options.trans_only) {
            // _balance_transonly's per-chromosome weights, the reciprocal of
            // the fraction of bins a chromosome does not cover.
            std::vector<double> cweights;
            cweights.reserve(width);
            for (std::size_t cid = 0; cid + 1 < offsets.size(); ++cid) {
                const std::int64_t lo = offsets[cid];
                const std::int64_t hi = offsets[cid + 1];
                const double weight =
                    1.0 - static_cast<double>(hi - lo) / static_cast<double>(n_bins);
                for (std::int64_t k = lo; k < hi; ++k) {
                    cweights.push_back(1.0 / weight);
                }
            }
            Filters filters = base;
            filters.zero_cis = true;
            const WholeStats stats =
                balance_whole(bias, pixels, spans, chrom_ids, filters, &cweights, options);
            scale = stats.scale;
            var = stats.var;
        } else {
            const WholeStats stats =
                balance_whole(bias, pixels, spans, chrom_ids, base, nullptr, options);
            scale = stats.scale;
            var = stats.var;
        }
    }

    json::Value stats = json::Value::object();
    stats["tol"] = options.tol;
    stats["min_nnz"] = options.min_nnz;
    stats["min_count"] = options.min_count;
    stats["mad_max"] = options.mad_max;
    stats["cis_only"] = options.cis_only;
    stats["ignore_diags"] = options.ignore_diags.as_stat();
    if (options.cis_only) {
        stats["scale"] = number_array(per_chrom.scales);
        stats["converged"] = bool_array(per_chrom.variances, options.tol);
        stats["var"] = number_array(per_chrom.variances);
    } else {
        stats["scale"] = scale;
        stats["converged"] = var < options.tol;
        stats["var"] = var;
    }
    stats["divisive_weights"] = false;

    if (options.store) {
        h5::File file(clr.store(), h5::Mode::ReadWrite);
        const std::string path = join_path(join_path(clr.root(), "bins"), options.store_name);
        if (file.exists(path)) {
            file.remove(path);
        }
        // grp["bins"].create_dataset(name, data=bias, compression="gzip",
        // compression_opts=6): gzip level 6 without shuffle, a fixed shape and
        // h5py's guessed chunk length.
        h5::DatasetCreate create;
        create.length = bias.size();
        create.chunk = h5::guess_chunk(bias.size(), sizeof(double));
        create.compression = h5::DatasetCreate::Compression::Gzip;
        create.gzip_level = 6;
        const h5::Handle type = h5::file_type(DType::Float64);
        h5::Dataset dataset = file.create_dataset(path, type.get(), create);
        dataset.write_column(0, Column(bias));
        for (const auto& [key, value] : stats.as_object()) {
            file.set_attribute(path, key, value);
        }
    }

    return {std::move(bias), std::move(stats)};
}

}  // namespace coolercpp
