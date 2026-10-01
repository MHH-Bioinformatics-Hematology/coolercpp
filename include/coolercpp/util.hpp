// cooler.util helpers that operate on bin tables.

#ifndef COOLERCPP_UTIL_HPP
#define COOLERCPP_UTIL_HPP

#include <cstdint>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "coolercpp/json.hpp"
#include "coolercpp/region.hpp"
#include "coolercpp/table.hpp"

namespace coolercpp {

// cooler.util.get_binsize: the single width shared by every bin except the
// last one of each chromosome, as a number carrying the dtype of end - start;
// null when the widths differ or no chromosome has more than one bin.
[[nodiscard]] json::Value get_binsize(const Table& bins);

// cooler.util.get_chromsizes: the end of the last bin of each chromosome, in
// the order of those last bins.
[[nodiscard]] ChromSizes get_chromsizes(const Table& bins);

// cooler.util.natsorted.
[[nodiscard]] std::vector<std::string> natsorted(std::vector<std::string> items);

// cooler.util.partition: [start, stop) cut into pieces of `step`, the last one
// short. Like range(), an empty interval yields nothing.
[[nodiscard]] std::vector<std::pair<std::int64_t, std::int64_t>> partition(std::int64_t start,
                                                                          std::int64_t stop,
                                                                          std::int64_t step);

// cooler.util.mad: the median absolute deviation from the median, over the
// whole array (the axis=None case). NaN for an empty array, as np.median is.
[[nodiscard]] double mad(std::span<const double> data);

}  // namespace coolercpp

#endif  // COOLERCPP_UTIL_HPP
