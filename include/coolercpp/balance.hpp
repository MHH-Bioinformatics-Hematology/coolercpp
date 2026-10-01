// Matrix balancing: cooler.balance_cooler and the marginal machinery it runs
// on (cooler/balance.py of cooler 0.10.2).

#ifndef COOLERCPP_BALANCE_HPP
#define COOLERCPP_BALANCE_HPP

#include <concepts>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "coolercpp/api.hpp"
#include "coolercpp/json.hpp"

namespace coolercpp {

// The ignore_diags argument of balance_cooler: a number of diagonals, or
// Python's False. The two are not interchangeable, because balance_cooler
// copies the argument into its stats dictionary, where cooler stores False as
// a bool and 0 as an int, and both end up as an attribute of the stored
// weight column.
class IgnoreDiags {
  public:
    IgnoreDiags(bool enabled) : value_(enabled ? 1 : 0), is_bool_(true) {}  // NOLINT
    template <std::integral I>
        requires(!std::same_as<I, bool>)
    IgnoreDiags(I n) : value_(static_cast<std::int64_t>(n)) {}  // NOLINT

    // Python truthiness: `if ignore_diags:` decides whether the filter runs.
    [[nodiscard]] bool enabled() const noexcept { return value_ != 0; }
    // The n_diags of cooler.balance._zero_diags.
    [[nodiscard]] std::int64_t n_diags() const noexcept { return value_; }
    // The value as it appears in the stats dictionary.
    [[nodiscard]] json::Value as_stat() const {
        return is_bool_ ? json::Value(value_ != 0) : json::Value(value_);
    }

  private:
    std::int64_t value_;
    bool is_bool_ = false;
};

// The keyword arguments of cooler.balance_cooler. `map` and `use_lock` serve
// Python's process pools; see the notes on `threads` and `use_lock` below.
struct BalanceOptions {
    bool cis_only = false;
    bool trans_only = false;
    IgnoreDiags ignore_diags = 2;
    std::int64_t mad_max = 5;
    std::int64_t min_nnz = 10;
    std::int64_t min_count = 0;
    // The bad bins to zero out. Negative indices count from the end, as they
    // do in the numpy assignment cooler performs.
    std::optional<std::vector<std::int64_t>> blacklist;
    bool rescale_marginals = true;
    // The initial weight vector. cooler replaces the NaNs in it with zeros and
    // keeps the array itself as the result; coolercpp copies it.
    std::optional<std::vector<double>> x0;
    double tol = 1e-5;
    std::int64_t max_iters = 200;
    // nullopt is Python's chunksize=None: all pixels in one chunk.
    std::optional<std::int64_t> chunksize = 10000000;
    // Where cooler takes a `map`, coolercpp takes a thread count. Threads
    // decode the HDF5 chunks of one pixel span; the marginals are still
    // accumulated span by span in file order, so the weights are the same
    // whatever the thread count is. The default is one thread, as cooler's
    // default builtin `map` is sequential.
    int threads = 1;
    // cooler's lock guards an HDF5 file shared with forked workers. coolercpp
    // reads the pixels from a single process and accepts the flag without
    // changing what it does.
    bool use_lock = false;
    bool store = false;
    std::string store_name = "weight";
};

// What balance_cooler returns: the bias vector and the stats dictionary, with
// the keys in the order cooler builds them. `scale`, `converged` and `var` are
// numbers for a genome-wide or trans-only run and arrays with one entry per
// chromosome for cis_only, as in Python.
struct BalanceResult {
    std::vector<double> bias;
    json::Value stats;
};

// cooler.balance_cooler (also exported as cooler.iterative_correction).
[[nodiscard]] BalanceResult balance_cooler(const Cooler& clr, const BalanceOptions& options = {});

}  // namespace coolercpp

#endif  // COOLERCPP_BALANCE_HPP
