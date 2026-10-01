// balance_cooler and the helpers it needs. The expected values were taken from
// cooler 0.10.2 on 2026-10-01; harness/cases/balance.json compares every
// option against cooler itself on real data.

#include <doctest/doctest.h>

#include <coolercpp/coolercpp.hpp>

#include <sys/stat.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>

#include "numpy_compat.hpp"

using namespace coolercpp;

namespace {

std::string data(const std::string& name) { return std::string(COOLERCPP_TEST_DATA) + "/" + name; }

std::string scratch(const std::string& name) {
    ::mkdir(COOLERCPP_TEST_SCRATCH, 0755);
    const std::string path = std::string(COOLERCPP_TEST_SCRATCH) + "/" + name;
    std::remove(path.c_str());
    return path;
}

void copy_file(const std::string& from, const std::string& to) {
    std::ifstream in(from, std::ios::binary);
    std::ofstream out(to, std::ios::binary);
    out << in.rdbuf();
}

double stat_double(const BalanceResult& result, const char* key) {
    return result.stats.find(key)->as_double();
}

}  // namespace

TEST_CASE("cooler.util partition and mad") {
    CHECK(partition(0, 9, 2) == std::vector<std::pair<std::int64_t, std::int64_t>>{
                                    {0, 2}, {2, 4}, {4, 6}, {6, 8}, {8, 9}});
    CHECK(partition(5, 5, 3).empty());
    CHECK(partition(10, 20, 100) ==
          std::vector<std::pair<std::int64_t, std::int64_t>>{{10, 20}});
    CHECK_THROWS_AS(partition(0, 10, 0), ValueError);

    // np.median(np.abs(x - np.median(x))) for x = [1, 1, 2, 2, 4, 6, 9]
    const std::vector<double> x{1, 1, 2, 2, 4, 6, 9};
    CHECK(npy::median(x) == 2.0);
    CHECK(mad(x) == 1.0);
    // An even count averages the two middle order statistics.
    const std::vector<double> even{4, 1, 3, 2};
    CHECK(npy::median(even) == 2.5);
    CHECK(std::isnan(npy::median(std::vector<double>{})));
    CHECK(std::isnan(npy::median(std::vector<double>{1.0, std::nan("")})));
    CHECK(npy::mean(x) == doctest::Approx(25.0 / 7.0));
    CHECK(npy::var(even) == doctest::Approx(1.25));
    CHECK(std::isnan(npy::mean(std::vector<double>{})));
}

TEST_CASE("ignore_diags keeps Python's int and False apart") {
    CHECK(IgnoreDiags(2).enabled());
    CHECK(IgnoreDiags(2).n_diags() == 2);
    CHECK(IgnoreDiags(2).as_stat().is_integer());
    CHECK_FALSE(IgnoreDiags(0).enabled());
    CHECK(IgnoreDiags(0).as_stat().is_integer());
    CHECK_FALSE(IgnoreDiags(false).enabled());
    CHECK(IgnoreDiags(false).as_stat().is_bool());
    CHECK(IgnoreDiags(true).enabled());
    CHECK(IgnoreDiags(true).n_diags() == 1);
}

TEST_CASE("balance_cooler on a toy cooler") {
    const Cooler c(data("toy.symm.upper.2.cool"));
    const BalanceResult result =
        balance_cooler(c, {.ignore_diags = 0, .mad_max = 0, .min_nnz = 0});
    REQUIRE(result.bias.size() == 32);
    CHECK(stat_double(result, "scale") == 6.0);
    CHECK(stat_double(result, "var") == 0.0);
    CHECK(result.stats.find("converged")->as_bool());
    CHECK(result.stats.find("cis_only")->as_bool() == false);
    CHECK(result.stats.find("divisive_weights")->as_bool() == false);
    CHECK(result.stats.find("min_nnz")->as_int() == 0);
    // The balanced matrix has rows summing to one, which is what
    // rescale_marginals promises.
    const DenseMatrix m = c.matrix({.balance = false}).all().dense();
    for (std::int64_t i = 0; i < 32; ++i) {
        double sum = 0.0;
        for (std::int64_t j = 0; j < 32; ++j) {
            sum += m.value(i, j) * result.bias[static_cast<std::size_t>(i)] *
                   result.bias[static_cast<std::size_t>(j)];
        }
        CHECK(sum == doctest::Approx(1.0).epsilon(1e-9));
    }

    // The default filters drop every bin of this 48 pixel cooler, which leaves
    // an empty marginal: cooler reports a NaN scale and claims convergence.
    const BalanceResult empty = balance_cooler(c);
    CHECK(std::isnan(stat_double(empty, "scale")));
    CHECK(stat_double(empty, "var") == 0.0);
    CHECK(result.stats.find("converged")->as_bool());
    for (const double weight : empty.bias) {
        CHECK(std::isnan(weight));
    }
}

TEST_CASE("balance_cooler does not depend on the thread count") {
    const Cooler c(data("yeast.10kb.cool"));
    const BalanceOptions base{.max_iters = 12, .chunksize = 40000};
    const BalanceResult one = balance_cooler(c, base);
    BalanceOptions many = base;
    many.threads = 4;
    const BalanceResult four = balance_cooler(c, many);
    REQUIRE(one.bias.size() == four.bias.size());
    for (std::size_t i = 0; i < one.bias.size(); ++i) {
        if (std::isnan(one.bias[i])) {
            CHECK(std::isnan(four.bias[i]));
        } else {
            // Bit for bit: the chunk marginals are folded in file order.
            CHECK(one.bias[i] == four.bias[i]);
        }
    }
    CHECK(stat_double(one, "scale") == stat_double(four, "scale"));
}

TEST_CASE("balance_cooler cis_only and trans_only") {
    const Cooler c(data("odd.1.cool"));
    const BalanceResult cis = balance_cooler(c, {.cis_only = true, .mad_max = 0, .min_nnz = 0});
    const json::Value* scales = cis.stats.find("scale");
    REQUIRE(scales->is_array());
    CHECK(scales->as_array().size() == 3);
    CHECK(cis.stats.find("converged")->as_array().size() == 3);
    CHECK(cis.stats.find("cis_only")->as_bool());

    const BalanceResult trans =
        balance_cooler(c, {.trans_only = true, .mad_max = 0, .min_nnz = 0, .max_iters = 20});
    CHECK(trans.stats.find("scale")->type() == json::Type::Double);
    // trans_only leaves cis_only False in the stats, as cooler does.
    CHECK(trans.stats.find("cis_only")->as_bool() == false);
}

TEST_CASE("balance_cooler option errors") {
    const Cooler c(data("toy.symm.upper.2.cool"));
    CHECK_THROWS_AS(balance_cooler(c, {.blacklist = std::vector<std::int64_t>{40}}), IndexError);
    CHECK_THROWS_AS(balance_cooler(c, {.mad_max = 0, .min_nnz = 0, .max_iters = 0}),
                    UnboundLocalError);
    CHECK_THROWS_AS(balance_cooler(c, {.x0 = std::vector<double>{1.0, 1.0}}), ValueError);
}

TEST_CASE("balance_cooler stores the weights") {
    const std::string path = scratch("balanced.cool");
    copy_file(data("toy.symm.upper.2.cool"), path);
    const Cooler c(path);
    const BalanceResult result = balance_cooler(
        c, {.ignore_diags = 0, .mad_max = 0, .min_nnz = 0, .store = true, .store_name = "ICE"});
    const Table bins = Cooler(path).bins().all();
    REQUIRE(bins.contains("ICE"));
    const std::vector<double> stored = bins["ICE"].values<double>();
    REQUIRE(stored.size() == result.bias.size());
    for (std::size_t i = 0; i < stored.size(); ++i) {
        CHECK(stored[i] == result.bias[i]);
    }
    // Storing again replaces the column instead of failing.
    CHECK_NOTHROW(balance_cooler(
        c, {.ignore_diags = 0, .mad_max = 0, .min_nnz = 0, .store = true, .store_name = "ICE"}));
    CHECK(Cooler(path).bins().all().contains("ICE"));
    // The weights the matrix selector applies are the stored ones.
    const DenseMatrix m = Cooler(path).matrix({.balance = "ICE"}).fetch("chr1").dense();
    CHECK(m.rows == 16);
}
