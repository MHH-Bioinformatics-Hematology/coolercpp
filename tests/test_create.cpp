// create_cooler: round trips, validation, multi-resolution writes and
// byte-identical output between runs.

#include <doctest/doctest.h>

#include <coolercpp/coolercpp.hpp>

#include <sys/stat.h>

#include <cstdio>
#include <fstream>
#include <iterator>
#include <string>

using namespace coolercpp;

namespace {

std::string scratch(const std::string& name) {
    ::mkdir(COOLERCPP_TEST_SCRATCH, 0755);
    const std::string path = std::string(COOLERCPP_TEST_SCRATCH) + "/" + name;
    std::remove(path.c_str());
    return path;
}

Table toy_bins() {
    return Table{{"chrom", Column::categorical({0, 0, 0, 1, 1}, {"chrA", "chrB"})},
                 {"start", Column{0, 10, 20, 0, 10}},
                 {"end", Column{10, 20, 25, 10, 17}},
                 {"weight", Column{1.0, 0.5, 2.0, 1.0, 0.25}}};
}

Table toy_pixels() {
    return Table{{"bin1_id", Column(std::vector<std::int64_t>{0, 0, 1, 3, 2})},
                 {"bin2_id", Column(std::vector<std::int64_t>{0, 2, 4, 4, 3})},
                 {"count", Column{5, 1, 2, 7, 3}}};
}

std::string slurp(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}

}  // namespace

TEST_CASE("a created cooler reads back") {
    const std::string path = scratch("roundtrip.cool");
    create_cooler(path, toy_bins(), toy_pixels(), {.assembly = "toy"});
    const Cooler c(path);
    CHECK(c.chromnames() == std::vector<std::string>{"chrA", "chrB"});
    CHECK(c.chromsizes().lengths() == std::vector<std::int64_t>{25, 17});
    CHECK(c.binsize() == 10);
    const json::Value info = c.info();
    CHECK(info.find("nnz")->as_int() == 5);
    CHECK(info.find("sum")->as_int() == 18);
    CHECK(info.find("genome-assembly")->as_string() == "toy");
    CHECK(info.find("format-version")->as_int() == 3);
    const Table pixels = c.pixels().all();
    // The input was sorted by (bin1_id, bin2_id) on the way in.
    CHECK(pixels["bin1_id"].values<std::int64_t>() == std::vector<std::int64_t>{0, 0, 1, 2, 3});
    CHECK(pixels["count"].values<std::int32_t>() == std::vector<std::int32_t>{5, 1, 2, 3, 7});
    const Table bins = c.bins().all();
    CHECK(bins.columns() == std::vector<std::string>{"chrom", "start", "end", "weight"});
    CHECK(bins["weight"].values<double>()[4] == 0.25);
    CHECK(c.matrix().all().dense().value(0, 2) == 1 * 1.0 * 2.0);
}

TEST_CASE("variable bins and float counts") {
    const std::string path = scratch("variable.cool");
    Table bins = toy_bins();
    bins.set("end", Column{10, 18, 25, 10, 17});
    Table pixels = toy_pixels();
    pixels.set("count", Column{0.1, 0.2, 0.3, 0.4, 0.5});
    create_cooler(path, bins, pixels, {.dtypes = {{"count", DType::Float64}}});
    const Cooler c(path);
    CHECK_FALSE(c.binsize().has_value());
    CHECK(c.info().find("bin-type")->as_string() == "variable");
    CHECK(c.info().find("sum")->as_double() == 1.5);
}

TEST_CASE("pixel validation raises BadInputError") {
    const std::string path = scratch("invalid.cool");
    Table duplicate = toy_pixels();
    duplicate.set("bin2_id", Column(std::vector<std::int64_t>{0, 0, 4, 4, 3}));
    duplicate.set("bin1_id", Column(std::vector<std::int64_t>{0, 0, 1, 3, 2}));
    CHECK_THROWS_AS(create_cooler(path, toy_bins(), duplicate), BadInputError);

    Table lower = toy_pixels();
    lower.set("bin1_id", Column(std::vector<std::int64_t>{0, 3, 1, 3, 2}));
    CHECK_THROWS_WITH_AS(create_cooler(path, toy_bins(), lower), "Found bin1_id greater than bin2_id",
                         BadInputError);

    Table bounds = toy_pixels();
    bounds.set("bin2_id", Column(std::vector<std::int64_t>{0, 2, 4, 4, 5}));
    CHECK_THROWS_AS(create_cooler(path, toy_bins(), bounds), BadInputError);

    Table no_chrom = toy_bins();
    no_chrom.drop("chrom");
    CHECK_THROWS_WITH_AS(create_cooler(path, no_chrom, toy_pixels()),
                         "Missing column from bin table: 'chrom'.", ValueError);
    CHECK_THROWS_WITH_AS(
        create_cooler(path, toy_bins(), toy_pixels(), {.h5opts = H5Opts{{"compresion", std::string("gzip")}}}),
        "Unknown storage option 'compresion'.", ValueError);
}

TEST_CASE("chunks, unordered input and resolutions in one file") {
    const std::string path = scratch("multi.mcool");
    const Table all = toy_pixels();
    const std::vector<std::int64_t> first{3, 4};
    const std::vector<std::int64_t> second{0, 1, 2};
    // Unordered: the external sort merges the chunks.
    create_cooler(path + "::/resolutions/10", toy_bins(),
                  iterate_chunks({all.take(first), all.take(second)}), {.mode = "a"});
    create_cooler(path + "::/resolutions/20", toy_bins(), toy_pixels(), {.mode = "a"});
    CHECK(list_coolers(path) == std::vector<std::string>{"/resolutions/10", "/resolutions/20"});
    const Cooler c(path + "::/resolutions/10");
    CHECK(c.pixels().all()["count"].values<std::int32_t>() == std::vector<std::int32_t>{5, 1, 2, 3, 7});
}

TEST_CASE("output is byte-identical from run to run") {
    const CreateOptions options{.metadata = json::Value::object(), .creation_date = "2026-01-01T00:00:00"};
    const std::string a = scratch("same_a.cool");
    const std::string b = scratch("same_b.cool");
    create_cooler(a, toy_bins(), toy_pixels(), options);
    create_cooler(b, toy_bins(), toy_pixels(), options);
    const std::string bytes_a = slurp(a);
    CHECK(!bytes_a.empty());
    CHECK(bytes_a == slurp(b));
}
