// Minimal .npy reading and writing for the harness: one or two dimensional,
// little endian, C order, the numeric dtypes of coolercpp::Column.

#ifndef COOLERCPP_HARNESS_NPY_HPP
#define COOLERCPP_HARNESS_NPY_HPP

#include <coolercpp/table.hpp>

#include <cstdint>
#include <string>
#include <vector>

namespace harness {

void write_npy(const std::string& path, const coolercpp::Column& column,
               const std::vector<std::int64_t>& shape);
void write_npy(const std::string& path, const coolercpp::Column& column);

// Reads a one dimensional array into a Column of the matching dtype.
[[nodiscard]] coolercpp::Column read_npy(const std::string& path);

}  // namespace harness

#endif  // COOLERCPP_HARNESS_NPY_HPP
