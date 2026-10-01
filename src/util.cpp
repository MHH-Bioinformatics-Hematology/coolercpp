// Port of cooler/util.py (cooler 0.10.2, BSD-3-Clause): get_binsize,
// get_chromsizes, natsort_key, natsorted, argnatsort, partition,
// read_chromsizes, binnify, check_bins, bedslice, rlencode, mad, buffered,
// GenomeSegmentation and balanced_partition.

#include "coolercpp/util.hpp"

#include <zlib.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <regex>
#include <set>
#include <sstream>
#include <unordered_map>
#include <variant>

#include "coolercpp/errors.hpp"

namespace coolercpp {

namespace {

// Groups the rows of the chrom column by label without materialising a label
// per row: a dense group id per row, and the label of every group.
struct ChromGroups {
    std::vector<std::int32_t> ids;
    std::vector<std::string> names;
};

ChromGroups chrom_groups(const Table& bins) {
    const Column& chrom = bins["chrom"];
    ChromGroups groups;
    groups.ids.resize(chrom.size());
    if (chrom.dtype() == DType::Categorical) {
        const CategoricalData& data = chrom.categorical();
        const std::size_t n = data.categories != nullptr ? data.categories->size() : 0;
        std::vector<std::int32_t> by_code(n + 1, -1);
        for (std::size_t i = 0; i < data.codes.size(); ++i) {
            const std::int32_t code = data.codes[i];
            const std::size_t slot =
                (code < 0 || static_cast<std::size_t>(code) >= n) ? n : static_cast<std::size_t>(code);
            if (by_code[slot] < 0) {
                by_code[slot] = static_cast<std::int32_t>(groups.names.size());
                groups.names.push_back(slot == n ? std::string() : (*data.categories)[slot]);
            }
            groups.ids[i] = by_code[slot];
        }
        return groups;
    }
    std::unordered_map<std::string, std::int32_t> index;
    for (std::size_t i = 0; i < chrom.size(); ++i) {
        std::string label = chrom.dtype() == DType::String ? chrom.values<std::string>()[i]
                                                          : std::to_string(chrom.as_int64(i));
        const auto [it, inserted] =
            index.emplace(std::move(label), static_cast<std::int32_t>(groups.names.size()));
        if (inserted) {
            groups.names.push_back(it->first);
        }
        groups.ids[i] = it->second;
    }
    return groups;
}

// The last row of every group.
std::vector<std::size_t> last_rows(const ChromGroups& groups) {
    std::vector<std::size_t> last(groups.names.size(), 0);
    for (std::size_t i = 0; i < groups.ids.size(); ++i) {
        last[static_cast<std::size_t>(groups.ids[i])] = i;
    }
    return last;
}

}  // namespace

json::Value get_binsize(const Table& bins) {
    const ChromGroups groups = chrom_groups(bins);
    const std::vector<std::size_t> last_row = last_rows(groups);
    const Column& start = bins["start"];
    const Column& end = bins["end"];
    const DType dtype = result_type(end.dtype(), start.dtype());
    // (group["end"] - group["start"]).iloc[:-1].unique() accumulated into a
    // set; the result only depends on the set of widths.
    std::set<double> float_widths;
    std::set<std::int64_t> int_widths;
    std::set<std::uint64_t> uint_widths;
    const bool floating = is_float(dtype);
    const bool unsigned_result = is_unsigned_integer(dtype);
    for (std::size_t i = 0; i < groups.ids.size(); ++i) {
        if (last_row[static_cast<std::size_t>(groups.ids[i])] == i) {
            continue;
        }
        if (floating) {
            double width = end.as_double(i) - start.as_double(i);
            if (dtype == DType::Float32) {
                width = static_cast<double>(static_cast<float>(end.as_double(i)) -
                                            static_cast<float>(start.as_double(i)));
            }
            float_widths.insert(width);
        } else if (unsigned_result) {
            uint_widths.insert(static_cast<std::uint64_t>(end.as_int64(i)) -
                               static_cast<std::uint64_t>(start.as_int64(i)));
        } else {
            std::int64_t width = end.as_int64(i) - start.as_int64(i);
            switch (dtype) {
                case DType::Int8: width = static_cast<std::int8_t>(width); break;
                case DType::Int16: width = static_cast<std::int16_t>(width); break;
                case DType::Int32: width = static_cast<std::int32_t>(width); break;
                default: break;
            }
            int_widths.insert(width);
        }
        if (float_widths.size() + int_widths.size() + uint_widths.size() > 1) {
            return json::Value();
        }
    }
    if (float_widths.size() == 1) {
        return json::Value::number(*float_widths.begin(), dtype);
    }
    if (int_widths.size() == 1) {
        return json::Value::integer(*int_widths.begin(), dtype);
    }
    if (uint_widths.size() == 1) {
        return json::Value::unsigned_integer(*uint_widths.begin(), dtype);
    }
    return json::Value();
}

ChromSizes get_chromsizes(const Table& bins) {
    const ChromGroups groups = chrom_groups(bins);
    const std::vector<std::size_t> last_row = last_rows(groups);
    const Column& end = bins["end"];
    std::vector<std::string> names;
    std::vector<std::int64_t> lengths;
    for (std::size_t i = 0; i < groups.ids.size(); ++i) {
        const auto group = static_cast<std::size_t>(groups.ids[i]);
        if (last_row[group] == i) {
            names.push_back(groups.names[group]);
            lengths.push_back(end.as_int64(i));
        }
    }
    return ChromSizes(std::move(names), std::move(lengths), end.dtype());
}

namespace {

// natsort_key: re.split(r"(\d+)", s), empty pieces dropped, digit runs as int.
using Piece = std::variant<std::string, std::uint64_t>;

std::vector<Piece> natsort_key(const std::string& s) {
    std::vector<Piece> key;
    std::size_t i = 0;
    while (i < s.size()) {
        std::size_t j = i;
        const bool digits = s[i] >= '0' && s[i] <= '9';
        while (j < s.size() && ((s[j] >= '0' && s[j] <= '9') == digits)) {
            ++j;
        }
        const std::string part = s.substr(i, j - i);
        if (digits) {
            key.emplace_back(static_cast<std::uint64_t>(std::stoull(part)));
        } else {
            key.emplace_back(part);
        }
        i = j;
    }
    return key;
}

}  // namespace

std::vector<std::string> natsorted(std::vector<std::string> items) {
    std::vector<std::pair<std::vector<Piece>, std::string>> keyed;
    keyed.reserve(items.size());
    for (std::string& item : items) {
        keyed.emplace_back(natsort_key(item), std::move(item));
    }
    std::stable_sort(keyed.begin(), keyed.end(), [](const auto& a, const auto& b) {
        const auto& x = a.first;
        const auto& y = b.first;
        const std::size_t n = std::min(x.size(), y.size());
        for (std::size_t k = 0; k < n; ++k) {
            if (x[k].index() != y[k].index()) {
                throw TypeError("'<' not supported between instances of 'int' and 'str'");
            }
            if (x[k] != y[k]) {
                return x[k] < y[k];
            }
        }
        return x.size() < y.size();
    });
    std::vector<std::string> out;
    out.reserve(keyed.size());
    for (auto& entry : keyed) {
        out.push_back(std::move(entry.second));
    }
    return out;
}

std::vector<std::int64_t> argnatsort(const std::vector<std::string>& names) {
    if (names.empty()) {
        return {};
    }
    std::vector<std::vector<Piece>> keys;
    keys.reserve(names.size());
    std::size_t common = std::numeric_limits<std::size_t>::max();
    for (const std::string& name : names) {
        keys.push_back(natsort_key(name));
        common = std::min(common, keys.back().size());
    }
    std::vector<std::int64_t> order(names.size());
    std::iota(order.begin(), order.end(), std::int64_t{0});
    std::stable_sort(order.begin(), order.end(), [&](std::int64_t a, std::int64_t b) {
        const auto& x = keys[static_cast<std::size_t>(a)];
        const auto& y = keys[static_cast<std::size_t>(b)];
        for (std::size_t k = 0; k < common; ++k) {
            if (x[k].index() != y[k].index()) {
                throw TypeError("'<' not supported between instances of 'int' and 'str'");
            }
            if (x[k] != y[k]) {
                return x[k] < y[k];
            }
        }
        return false;
    });
    return order;
}

std::vector<std::pair<std::int64_t, std::int64_t>> partition(std::int64_t start, std::int64_t stop,
                                                             std::int64_t step) {
    if (step <= 0) {
        throw ValueError("range() arg 3 must not be zero");
    }
    std::vector<std::pair<std::int64_t, std::int64_t>> out;
    for (std::int64_t i = start; i < stop; i += step) {
        out.emplace_back(i, std::min(i + step, stop));
    }
    return out;
}

namespace {

// The two first tab separated fields of every non empty line, as pandas
// read_csv with usecols=[0, 1] and no header reads them.
ChromSizes parse_chromsizes(std::istream& input, const ReadChromsizesOptions& options) {
    std::vector<std::string> names;
    std::vector<std::int64_t> lengths;
    std::string line;
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty()) {
            continue;
        }
        const std::size_t first = line.find('\t');
        if (first == std::string::npos) {
            throw ValueError(
                "Usecols do not match columns, columns expected but not found: [1]");
        }
        const std::size_t second = line.find('\t', first + 1);
        const std::string value =
            line.substr(first + 1, second == std::string::npos ? std::string::npos
                                                               : second - first - 1);
        std::size_t consumed = 0;
        std::int64_t length = 0;
        try {
            length = std::stoll(value, &consumed);
        } catch (const std::exception&) {
            throw ValueError("invalid literal for int() with base 10: '" + value + "'");
        }
        if (consumed != value.size()) {
            throw ValueError("invalid literal for int() with base 10: '" + value + "'");
        }
        names.push_back(line.substr(0, first));
        lengths.push_back(length);
    }
    if (options.all_names) {
        return ChromSizes(std::move(names), std::move(lengths), DType::Int64);
    }
    // One group per pattern, each sorted in natural order and kept in the
    // order of the patterns. A name that matches two patterns is kept twice.
    std::vector<std::string> kept_names;
    std::vector<std::int64_t> kept_lengths;
    for (const std::string& pattern : options.name_patterns) {
        std::regex regex;
        try {
            regex = std::regex(pattern, std::regex::ECMAScript);
        } catch (const std::regex_error&) {
            throw ValueError("could not compile the pattern '" + pattern + "'");
        }
        std::vector<std::string> part_names;
        std::vector<std::int64_t> part_lengths;
        for (std::size_t i = 0; i < names.size(); ++i) {
            if (std::regex_search(names[i], regex)) {
                part_names.push_back(names[i]);
                part_lengths.push_back(lengths[i]);
            }
        }
        for (const std::int64_t position : argnatsort(part_names)) {
            kept_names.push_back(part_names[static_cast<std::size_t>(position)]);
            kept_lengths.push_back(part_lengths[static_cast<std::size_t>(position)]);
        }
    }
    return ChromSizes(std::move(kept_names), std::move(kept_lengths), DType::Int64);
}

}  // namespace

ChromSizes read_chromsizes(std::istream& input, const ReadChromsizesOptions& options) {
    return parse_chromsizes(input, options);
}

ChromSizes read_chromsizes(const std::string& filepath, const ReadChromsizesOptions& options) {
    std::string text;
    if (filepath.size() > 3 && filepath.compare(filepath.size() - 3, 3, ".gz") == 0) {
        gzFile file = gzopen(filepath.c_str(), "rb");
        if (file == nullptr) {
            throw OSError("[Errno 2] No such file or directory: '" + filepath + "'");
        }
        char buffer[1 << 16];
        int read = 0;
        while ((read = gzread(file, buffer, sizeof(buffer))) > 0) {
            text.append(buffer, static_cast<std::size_t>(read));
        }
        const bool failed = read < 0;
        gzclose(file);
        if (failed) {
            throw OSError("cannot decompress " + filepath);
        }
    } else {
        std::ifstream in(filepath, std::ios::binary);
        if (!in) {
            throw OSError("[Errno 2] No such file or directory: '" + filepath + "'");
        }
        std::ostringstream buffer;
        buffer << in.rdbuf();
        text = buffer.str();
    }
    std::istringstream input(text);
    return parse_chromsizes(input, options);
}

Table binnify(const ChromSizes& chromsizes, std::int64_t binsize) {
    if (binsize == 0) {
        throw ZeroDivisionError("division by zero");
    }
    std::vector<std::int32_t> codes;
    std::vector<std::int64_t> starts;
    std::vector<std::int64_t> ends;
    for (std::size_t c = 0; c < chromsizes.size(); ++c) {
        const std::int64_t clen = chromsizes.lengths()[c];
        const auto n_bins = static_cast<std::int64_t>(
            std::ceil(static_cast<double>(clen) / static_cast<double>(binsize)));
        for (std::int64_t k = 0; k < n_bins; ++k) {
            codes.push_back(static_cast<std::int32_t>(c));
            starts.push_back(k * binsize);
            // The last edge of every sequence is its length.
            ends.push_back(k + 1 == n_bins ? clen : (k + 1) * binsize);
        }
    }
    Table bins;
    bins.set("chrom", Column::categorical(std::move(codes), chromsizes.names(), true));
    bins.set("start", Column(std::move(starts)));
    bins.set("end", Column(std::move(ends)));
    return bins;
}

Table check_bins(const Table& bins, const ChromSizes& chromsizes) {
    const Column& chrom = bins["chrom"];
    Table out = bins;
    if (chrom.dtype() == DType::Categorical) {
        if (*chrom.categorical().categories != chromsizes.names()) {
            throw AssertionError("");
        }
        return out;
    }
    std::unordered_map<std::string, std::int32_t> idmap;
    for (std::size_t c = 0; c < chromsizes.size(); ++c) {
        idmap.emplace(chromsizes.names()[c], static_cast<std::int32_t>(c));
    }
    std::vector<std::int32_t> codes(chrom.size(), -1);
    if (chrom.dtype() == DType::String) {
        const std::vector<std::string>& labels = chrom.values<std::string>();
        for (std::size_t i = 0; i < labels.size(); ++i) {
            const auto found = idmap.find(labels[i]);
            if (found != idmap.end()) {
                codes[i] = found->second;
            }
        }
    }
    // A column of anything but labels matches no category, which is the -1
    // pandas fills in.
    out.set("chrom", Column::categorical(std::move(codes), chromsizes.names(), true));
    return out;
}

namespace {

// numpy searchsorted on an ordered integer column.
std::size_t searchsorted(const Column& column, std::size_t lo, std::size_t hi, std::int64_t value,
                         bool right) {
    while (lo < hi) {
        const std::size_t mid = lo + (hi - lo) / 2;
        const std::int64_t at = column.as_int64(mid);
        if (right ? at <= value : at < value) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    return lo;
}

// groupby(chrom).get_group(name): the rows of one sequence, in file order.
Table chrom_group(const Table& bed, const std::string& name) {
    const Column& chrom = bed["chrom"];
    std::vector<std::int64_t> positions;
    for (std::size_t i = 0; i < chrom.size(); ++i) {
        if (chrom.label(i) == name) {
            positions.push_back(static_cast<std::int64_t>(i));
        }
    }
    if (positions.empty()) {
        throw KeyError(name);
    }
    return bed.take(positions);
}

// The [start, end) trim bedslice and GenomeSegmentation.fetch share.
Table trim_to_region(Table result, std::int64_t start, std::int64_t end) {
    const std::size_t n = result.num_rows();
    const std::size_t lo = searchsorted(result["end"], 0, n, start, true);
    const std::size_t hi = lo + searchsorted(result["start"], lo, n, end, false) - lo;
    return result.slice(lo, hi);
}

}  // namespace

Table bedslice(const Table& bed, const ChromSizes& chromsizes, const Region& region) {
    const RegionTuple parsed = parse_region(region, &chromsizes);
    Table result = chrom_group(bed, parsed.chrom);
    if (parsed.start > 0 || parsed.end < chromsizes[parsed.chrom]) {
        result = trim_to_region(std::move(result), parsed.start, parsed.end);
    }
    return result;
}

RunLengths rlencode(const Column& array, std::optional<std::int64_t> chunksize) {
    const std::size_t n = array.size();
    RunLengths out;
    if (n == 0) {
        out.values = array.slice(0, 0);
        return out;
    }
    const bool floating = is_float(array.dtype());
    const bool labelled = array.dtype() == DType::String;
    // NaN compares unequal to itself, so a NaN always starts a new run.
    const auto same = [&](std::size_t a, std::size_t b) {
        if (labelled) {
            return array.label(a) == array.label(b);
        }
        if (floating) {
            const double x = array.as_double(a);
            const double y = array.as_double(b);
            return x == y;
        }
        return array.as_int64(a) == array.as_int64(b);
    };
    const std::size_t step = chunksize.has_value() && *chunksize > 0
                                 ? static_cast<std::size_t>(*chunksize)
                                 : n;
    bool has_previous = false;
    std::size_t previous = 0;
    for (std::size_t lo = 0; lo < n; lo += step) {
        const std::size_t hi = std::min(n, lo + step);
        // The chunk's first element continues the previous chunk's last run
        // only when it equals it; cooler starts from NaN, so the very first
        // element always opens a run.
        if (!has_previous || !same(lo, previous)) {
            out.starts.push_back(static_cast<std::int64_t>(lo));
        }
        for (std::size_t k = lo + 1; k < hi; ++k) {
            if (!same(k, k - 1)) {
                out.starts.push_back(static_cast<std::int64_t>(k));
            }
        }
        previous = hi - 1;
        has_previous = true;
    }
    out.lengths.reserve(out.starts.size());
    for (std::size_t k = 0; k < out.starts.size(); ++k) {
        const std::int64_t next =
            k + 1 < out.starts.size() ? out.starts[k + 1] : static_cast<std::int64_t>(n);
        out.lengths.push_back(next - out.starts[k]);
    }
    out.values = array.take(out.starts);
    return out;
}

namespace {

// np.median: the average of the two middle elements for an even count.
double median(std::vector<double> values) {
    if (values.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const std::size_t half = values.size() / 2;
    std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(half),
                     values.end());
    const double upper = values[half];
    if (values.size() % 2 == 1) {
        return upper;
    }
    const double lower =
        *std::max_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(half));
    return (lower + upper) / 2.0;
}

}  // namespace

double mad(const Column& data) {
    const std::vector<double> values = data.as<double>();
    const double centre = median(values);
    std::vector<double> deviations(values.size());
    for (std::size_t i = 0; i < values.size(); ++i) {
        deviations[i] = std::abs(values[i] - centre);
    }
    return median(std::move(deviations));
}

std::function<std::optional<Table>()> buffered(std::function<std::optional<Table>()> chunks,
                                              std::int64_t size) {
    struct State {
        std::function<std::optional<Table>()> chunks;
        std::int64_t size = 0;
        std::vector<Table> buffer;
        bool done = false;
    };
    auto state = std::make_shared<State>();
    state->chunks = std::move(chunks);
    state->size = size;
    return [state]() -> std::optional<Table> {
        // pd.concat(buf, axis=0): the columns of the first chunk, extended by
        // the others, and the input indexes concatenated as they are.
        const auto flush = [&state]() -> Table {
            const std::vector<std::string> names = state->buffer.front().columns();
            std::vector<Column> columns;
            for (std::size_t c = 0; c < names.size(); ++c) {
                columns.push_back(state->buffer.front().column(c));
            }
            std::vector<std::int64_t> labels = state->buffer.front().index().to_vector();
            for (std::size_t k = 1; k < state->buffer.size(); ++k) {
                const Table& more = state->buffer[k];
                if (more.columns() != names) {
                    throw ValueError("buffered() needs chunks with the same columns");
                }
                for (std::size_t c = 0; c < names.size(); ++c) {
                    columns[c].append(more.column(c));
                }
                const std::vector<std::int64_t> more_labels = more.index().to_vector();
                labels.insert(labels.end(), more_labels.begin(), more_labels.end());
            }
            Table out;
            for (std::size_t c = 0; c < names.size(); ++c) {
                out.set(names[c], std::move(columns[c]));
            }
            out.set_index(Index::values(std::move(labels)));
            state->buffer.clear();
            return out;
        };
        std::int64_t rows = 0;
        for (const Table& chunk : state->buffer) {
            rows += static_cast<std::int64_t>(chunk.num_rows());
        }
        while (!state->done) {
            std::optional<Table> chunk = state->chunks();
            if (!chunk.has_value()) {
                state->done = true;
                break;
            }
            rows += static_cast<std::int64_t>(chunk->num_rows());
            state->buffer.push_back(std::move(*chunk));
            if (rows > state->size) {
                return flush();
            }
        }
        if (state->buffer.empty()) {
            return std::nullopt;
        }
        return flush();
    };
}

GenomeSegmentation::GenomeSegmentation(ChromSizes chromsizes, const Table& bins)
    : chromsizes_(std::move(chromsizes)), bins_(check_bins(bins, chromsizes_)) {
    binsize_ = get_binsize(bins_);
    // groupby("chrom", observed=True, sort=False): the sequences the bin table
    // carries, in the order they appear, each a contiguous run.
    const Column& chrom = bins_["chrom"];
    const std::size_t n = chrom.size();
    for (std::size_t i = 0; i < n;) {
        const std::string name = chrom.label(i);
        std::size_t j = i;
        while (j < n && chrom.label(j) == name) {
            ++j;
        }
        observed_.push_back(name);
        extents_.emplace_back(i, j);
        i = j;
    }
    chrom_binoffset_.push_back(0);
    for (const auto& [lo, hi] : extents_) {
        chrom_binoffset_.push_back(chrom_binoffset_.back() +
                                   static_cast<std::int64_t>(hi - lo));
    }
    chrom_abspos_.push_back(0);
    for (const std::int64_t length : chromsizes_.lengths()) {
        chrom_abspos_.push_back(chrom_abspos_.back() + length);
    }
    const Column& start = bins_["start"];
    start_abspos_.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        const std::int32_t code = chrom.categorical().codes[i];
        start_abspos_[i] = chrom_abspos_[static_cast<std::size_t>(code)] + start.as_int64(i);
    }
}

std::int64_t GenomeSegmentation::idmap(const std::string& contig) const {
    for (std::size_t c = 0; c < chromsizes_.size(); ++c) {
        if (chromsizes_.names()[c] == contig) {
            return static_cast<std::int64_t>(c);
        }
    }
    throw KeyError(contig);
}

Table GenomeSegmentation::group(const std::string& contig) const {
    for (std::size_t k = 0; k < observed_.size(); ++k) {
        if (observed_[k] == contig) {
            return bins_.slice(extents_[k].first, extents_[k].second);
        }
    }
    throw KeyError(contig);
}

Table GenomeSegmentation::fetch(const Region& region) const {
    const RegionTuple parsed = parse_region(region, &chromsizes_);
    Table result = group(parsed.chrom);
    if (parsed.start > 0 || parsed.end < chromsizes_[parsed.chrom]) {
        result = trim_to_region(std::move(result), parsed.start, parsed.end);
    }
    return result;
}

std::vector<RegionTuple> balanced_partition(const GenomeSegmentation& gs,
                                            std::int64_t n_chunk_max,
                                            const std::vector<std::string>& file_contigs,
                                            const std::vector<std::pair<std::string, double>>&
                                                loadings) {
    const std::vector<std::string>& observed = gs.observed();
    std::vector<double> weights;
    if (loadings.empty()) {
        for (std::size_t k = 0; k < observed.size(); ++k) {
            weights.push_back(static_cast<double>(gs.chrom_binoffset()[k + 1] -
                                                  gs.chrom_binoffset()[k]));
        }
    } else {
        for (const std::string& name : observed) {
            const auto found = std::find_if(
                loadings.begin(), loadings.end(),
                [&name](const auto& entry) { return entry.first == name; });
            if (found == loadings.end()) {
                throw KeyError(name);
            }
            weights.push_back(found->second);
        }
    }
    if (weights.empty()) {
        return {};
    }
    // idxmax(): the first sequence reaching the largest weight.
    std::size_t heaviest = 0;
    for (std::size_t k = 1; k < weights.size(); ++k) {
        if (weights[k] > weights[heaviest]) {
            heaviest = k;
        }
    }
    const double reference = weights[heaviest];
    const double constant =
        static_cast<double>(gs.chrom_binoffset()[heaviest + 1] -
                            gs.chrom_binoffset()[heaviest]) /
        static_cast<double>(n_chunk_max);
    std::vector<RegionTuple> granges;
    for (std::size_t k = 0; k < observed.size(); ++k) {
        const std::string& contig = observed[k];
        if (std::find(file_contigs.begin(), file_contigs.end(), contig) == file_contigs.end()) {
            continue;
        }
        const std::int64_t clen = gs.chromsizes()[contig];
        const auto step =
            static_cast<std::int64_t>(std::ceil(constant / (weights[k] / reference)));
        if (step <= 0) {
            throw ValueError("slice step cannot be zero");
        }
        const Table rows = gs.group(contig);
        const Column& start = rows["start"];
        std::vector<std::int64_t> anchors;
        for (std::size_t i = 0; i < start.size(); i += static_cast<std::size_t>(step)) {
            anchors.push_back(start.as_int64(i));
        }
        if (anchors.empty()) {
            throw IndexError("index -1 is out of bounds for axis 0 with size 0");
        }
        if (anchors.back() != clen) {
            anchors.push_back(clen);
        }
        for (std::size_t i = 0; i + 1 < anchors.size(); ++i) {
            granges.push_back(RegionTuple{contig, anchors[i], anchors[i + 1]});
        }
    }
    return granges;
}

}  // namespace coolercpp

#include <cstdio>
#include <mutex>

namespace coolercpp {

namespace {

std::mutex& warning_mutex() {
    static std::mutex mutex;
    return mutex;
}

WarningHandler& warning_handler() {
    static WarningHandler handler = [](const std::string& message) {
        std::fprintf(stderr, "UserWarning: %s\n", message.c_str());
    };
    return handler;
}

}  // namespace

void set_warning_handler(WarningHandler handler) {
    const std::lock_guard<std::mutex> lock(warning_mutex());
    warning_handler() = handler ? std::move(handler) : [](const std::string&) {};
}

void warn(const std::string& message) {
    WarningHandler handler;
    {
        const std::lock_guard<std::mutex> lock(warning_mutex());
        handler = warning_handler();
    }
    handler(message);
}

}  // namespace coolercpp
