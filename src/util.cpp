// Port of cooler/util.py get_binsize, get_chromsizes and natsorted
// (cooler 0.10.2, BSD-3-Clause).

#include "coolercpp/util.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <unordered_map>
#include <variant>

#include "coolercpp/errors.hpp"

namespace coolercpp {

namespace {

std::vector<std::string> chrom_labels(const Table& bins) {
    const Column& chrom = bins["chrom"];
    if (chrom.dtype() == DType::String || chrom.dtype() == DType::Categorical) {
        return chrom.labels();
    }
    // A numeric chrom column groups by value.
    std::vector<std::string> out(chrom.size());
    for (std::size_t i = 0; i < out.size(); ++i) {
        out[i] = std::to_string(chrom.as_int64(i));
    }
    return out;
}

}  // namespace

json::Value get_binsize(const Table& bins) {
    const std::vector<std::string> labels = chrom_labels(bins);
    const Column& start = bins["start"];
    const Column& end = bins["end"];
    const DType dtype = result_type(end.dtype(), start.dtype());
    std::unordered_map<std::string, std::size_t> last_row;
    for (std::size_t i = 0; i < labels.size(); ++i) {
        last_row[labels[i]] = i;
    }
    // (group["end"] - group["start"]).iloc[:-1].unique() accumulated into a
    // set; the result only depends on the set of widths.
    std::set<double> float_widths;
    std::set<std::int64_t> int_widths;
    std::set<std::uint64_t> uint_widths;
    const bool floating = is_float(dtype);
    const bool unsigned_result = is_unsigned_integer(dtype);
    for (std::size_t i = 0; i < labels.size(); ++i) {
        if (last_row[labels[i]] == i) {
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
    const std::vector<std::string> labels = chrom_labels(bins);
    const Column& end = bins["end"];
    std::unordered_map<std::string, std::size_t> last_row;
    for (std::size_t i = 0; i < labels.size(); ++i) {
        last_row[labels[i]] = i;
    }
    std::vector<std::string> names;
    std::vector<std::int64_t> lengths;
    for (std::size_t i = 0; i < labels.size(); ++i) {
        if (last_row[labels[i]] == i) {
            names.push_back(labels[i]);
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
