// Reading coolers written by Python cooler (cooler's own test files). The
// expected values were taken from cooler 0.10.2 on 2026-09-13; the harness
// compares every selector against cooler itself on real data.

#include <doctest/doctest.h>

#include <coolercpp/coolercpp.hpp>

#include <string>

using namespace coolercpp;

namespace {

std::string data(const std::string& name) { return std::string(COOLERCPP_TEST_DATA) + "/" + name; }

}  // namespace

TEST_CASE("Cooler metadata") {
    const Cooler c(data("toy.symm.upper.2.cool"));
    CHECK(c.chromnames() == std::vector<std::string>{"chr1", "chr2"});
    CHECK(c.chromsizes().lengths() == std::vector<std::int64_t>{32, 32});
    CHECK(c.chromsizes().dtype() == DType::Int32);
    CHECK(c.binsize() == 2);
    CHECK(c.storage_mode() == "symmetric-upper");
    CHECK(c.shape() == std::array<std::int64_t, 2>{32, 32});
    CHECK(c.info().find("nnz")->as_int() == 48);
    CHECK(c.info().find("sum")->as_int() == 96);
    CHECK(c.offset("chr2") == 16);
    CHECK(c.extent("chr1:3-9") == std::pair<std::int64_t, std::int64_t>{1, 5});
    CHECK(c.root() == "/");
    CHECK(Cooler(data("toy.asymm.2.cool")).storage_mode() == "square");
}

TEST_CASE("table selectors") {
    const Cooler c(data("toy.symm.upper.2.cool"));
    const Table pixels = c.pixels().slice(0, 3);
    CHECK(pixels.columns() == std::vector<std::string>{"bin1_id", "bin2_id", "count"});
    CHECK(pixels["bin2_id"].values<std::int64_t>() == std::vector<std::int64_t>{8, 16, 31});
    const Table chr2 = c.bins().fetch("chr2");
    CHECK(chr2.index().start() == 16);
    CHECK(chr2["chrom"].dtype() == DType::Categorical);
    CHECK(chr2["chrom"].label(0) == "chr2");
    CHECK(chr2["end"].values<std::int32_t>()[1] == 4);
    CHECK(c.bins()["start"].slice(0, 2).num_columns() == 1);
    CHECK_THROWS_AS((void)c.bins()["start"].columns(), AttributeError);
    CHECK_THROWS_AS((void)c.chroms().fetch("chr1"), NotImplementedError);
    CHECK_THROWS_AS((void)c.bins()[32], IndexError);

    const Table joined = c.pixels(true).slice(0, 2);
    CHECK(joined.columns() ==
          std::vector<std::string>{"chrom1", "start1", "end1", "chrom2", "start2", "end2", "count"});
}

TEST_CASE("matrix selector fills the lower triangle") {
    const Cooler c(data("toy.symm.upper.2.cool"));
    const MatrixResult dense = c.matrix({.balance = false}).all();
    REQUIRE(dense.is_dense());
    const DenseMatrix& m = dense.dense();
    std::int64_t total = 0;
    for (std::int64_t i = 0; i < m.rows; ++i) {
        for (std::int64_t j = 0; j < m.cols; ++j) {
            CHECK(m.value(i, j) == m.value(j, i));
            total += static_cast<std::int64_t>(m.value(i, j));
        }
    }
    CHECK(total == 192);
    CHECK(c.matrix({.balance = false, .sparse = true}).all().sparse().nnz() == 96);
    CHECK_THROWS_AS((void)c.matrix().all(), ValueError);
}

TEST_CASE("balanced values are count * w1 * w2") {
    const Cooler c(data("yeast.10kb.cool"));
    const DenseMatrix m = c.matrix().fetch(Region("chrMito", 0, 30000)).dense();
    CHECK(m.value(0, 2) == 0.03162385247691608);
    CHECK(m.value(2, 2) == 0.047286293435556334);
    const std::vector<double> w = c.bins()["weight"].slice(0, 3)["weight"].values<double>();
    CHECK(m.value(0, 2) == 4665 * (w[0] * w[2]));
}

TEST_CASE("multi-resolution files and errors") {
    CHECK(list_coolers(data("toy.symm.upper.2.mcool")) ==
          std::vector<std::string>{"/resolutions/2", "/resolutions/4", "/resolutions/8",
                                   "/resolutions/16", "/resolutions/32"});
    const Cooler res(data("toy.symm.upper.2.mcool") + "::resolutions/4");
    CHECK(res.binsize() == 4);
    try {
        const Cooler whole(data("toy.symm.upper.2.mcool"));
        FAIL("opening a multi-resolution file without a group must raise");
    } catch (const KeyError& e) {
        CHECK(std::string(e.what()).find("Coolers found in ['/resolutions/2', '/resolutions/4'") !=
              std::string::npos);
    }
    CHECK_THROWS_AS(Cooler(data("does_not_exist.cool")), OSError);
    CHECK_THROWS_AS((void)list_coolers(data("toy.chrom.sizes")), OSError);
}

TEST_CASE("chromosome IDs stored without an enum header") {
    const Cooler c(data("manycontigs.1.cool"));
    CHECK(c.chromnames().size() == 10000);
    const Table bins = c.bins().slice(0, 2);
    CHECK(bins["chrom"].dtype() == DType::Categorical);
    CHECK(bins["chrom"].label(1) == "contig0000");
    // convert_enum=False leaves the plain integer IDs alone.
    CHECK(c.bins(false).slice(0, 2)["chrom"].dtype() == DType::Int32);
}
