// create_scool and cooler.fileops: the single-cell layout, the file
// predicates, the group listings and the copy, move and link operations.

#include <doctest/doctest.h>

#include <coolercpp/coolercpp.hpp>

#include <sys/stat.h>

#include <algorithm>
#include <cstdio>
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

Table pixels(std::size_t n) {
    const std::vector<std::int64_t> bin1{0, 0, 1, 2, 3};
    const std::vector<std::int64_t> bin2{0, 2, 4, 3, 4};
    const std::vector<std::int32_t> count{5, 1, 2, 3, 7};
    return Table{{"bin1_id", Column(std::vector<std::int64_t>(bin1.begin(), bin1.begin() + n))},
                 {"bin2_id", Column(std::vector<std::int64_t>(bin2.begin(), bin2.begin() + n))},
                 {"count", Column(std::vector<std::int32_t>(count.begin(), count.begin() + n))}};
}

std::string write_scool(const std::string& name) {
    const std::string path = scratch(name);
    create_scool(path, toy_bins(),
                 {{"cell10", pixels(5)}, {"cell2", pixels(3)}, {"cell1", pixels(1)}},
                 {.assembly = "toy", .creation_date = "2026-01-01T00:00:00"});
    return path;
}

}  // namespace

TEST_CASE("create_scool writes the single-cell layout") {
    const std::string path = write_scool("layout.scool");

    CHECK(is_scool_file(path));
    CHECK_FALSE(is_cooler(path));
    CHECK_FALSE(is_multires_file(path));
    // The cells are listed in natural order, not in the order they were given
    // and not in the lexicographic order create_scool writes them in.
    CHECK(list_scool_cells(path) ==
          std::vector<std::string>{"/cells/cell1", "/cells/cell2", "/cells/cell10"});
    CHECK(list_coolers(path) == list_scool_cells(path));

    // The root carries the shared tables and the scool magic, and no pixels.
    const std::vector<std::string> entries = ls(path);
    CHECK(std::find(entries.begin(), entries.end(), "/bins/chrom") != entries.end());
    CHECK(std::find(entries.begin(), entries.end(), "/chroms/name") != entries.end());
    CHECK(std::find(entries.begin(), entries.end(), "/pixels") == entries.end());
    // The cell shares the three main bin columns and keeps its own weight.
    CHECK(std::find(entries.begin(), entries.end(), "/cells/cell1/bins/start") != entries.end());
    CHECK(std::find(entries.begin(), entries.end(), "/cells/cell1/bins/weight") != entries.end());
}

TEST_CASE("a scool cell reads back as a cooler") {
    const std::string path = write_scool("read.scool");
    const Cooler cell(path + "::/cells/cell2");
    CHECK(cell.chromnames() == std::vector<std::string>{"chrA", "chrB"});
    CHECK(cell.binsize() == 10);
    CHECK(cell.shape() == std::array<std::int64_t, 2>{5, 5});
    CHECK(cell.info().find("nnz")->as_int() == 3);
    CHECK(cell.info().find("sum")->as_int() == 8);
    CHECK(cell.info().find("format")->as_string() == "HDF5::Cooler");
    CHECK(cell.info().find("genome-assembly")->as_string() == "toy");
    const Table bins = cell.bins().all();
    CHECK(bins.num_rows() == 5);
    CHECK(bins["weight"].values<double>()[1] == doctest::Approx(0.5));
    CHECK(cell.pixels().all()["count"].values<std::int32_t>() ==
          std::vector<std::int32_t>{5, 1, 2});
}

TEST_CASE("the scool root info counts the cells") {
    const std::string path = write_scool("info.scool");
    const json::Value info = Cooler(path + "::/cells/cell1").info();
    CHECK(info.find("nbins")->as_int() == 5);
    // The root is not a cooler, so its attributes are read through fileops.
    CHECK(is_scool_file(path));
    const std::string tree = pprint_data_tree(path, 1);
    CHECK(tree.find("/\n") == 0);
    CHECK(tree.find("bins") != std::string::npos);
    // Depth 1 stops before the contents of the groups.
    CHECK(tree.find("chrom (5,) int32") == std::string::npos);
}

TEST_CASE("create_scool with one bin table per cell") {
    const std::string path = scratch("percell.scool");
    Table other = toy_bins();
    other["weight"] = Column{2.0, 1.0, 4.0, 2.0, 0.5};
    create_scool(path, CellTables{{"a", toy_bins()}, {"b", other}},
                 {{"a", pixels(5)}, {"b", pixels(2)}});
    CHECK(list_scool_cells(path) == std::vector<std::string>{"/cells/a", "/cells/b"});
    // The shared root table keeps only chrom, start and end.
    const std::vector<std::string> entries = ls(path);
    CHECK(std::find(entries.begin(), entries.end(), "/bins/weight") == entries.end());
    CHECK(Cooler(path + "::/cells/b").bins().all()["weight"].values<double>()[0] ==
          doctest::Approx(2.0));
    CHECK(Cooler(path + "::/cells/a").bins().all()["weight"].values<double>()[0] ==
          doctest::Approx(1.0));
}

TEST_CASE("create_scool rejects bad input") {
    const std::string path = scratch("bad.scool");
    CHECK_THROWS_AS(create_scool(path, CellTables{}, {{"a", pixels(1)}}), ValueError);
    CHECK_THROWS_AS(create_scool(path, CellTables{{"a", toy_bins()}, {"c", toy_bins()}},
                                 {{"a", pixels(1)}, {"b", pixels(1)}}),
                    ValueError);
    Table without_end = toy_bins();
    without_end.drop("end");
    CHECK_THROWS_AS(create_scool(path, without_end, {{"a", pixels(1)}}), ValueError);
    // The cell URI is built by appending "::/cells/<name>", so a URI that
    // already names a group leaves three parts behind.
    CHECK_THROWS_AS(create_scool(path + "::/experiment", toy_bins(), {{"a", pixels(1)}}),
                    ValueError);
}

TEST_CASE("fileops predicates recognise the file kinds") {
    const std::string cool = scratch("plain.cool");
    create_cooler(cool, toy_bins(), pixels(5));
    CHECK(is_cooler(cool));
    CHECK_FALSE(is_multires_file(cool));
    CHECK_FALSE(is_scool_file(cool));
    CHECK_THROWS_AS((void)list_scool_cells(cool), OSError);

    const std::string mcool = scratch("two.mcool");
    create_cooler(mcool + "::/resolutions/10", toy_bins(), pixels(5), {.mode = "a"});
    create_cooler(mcool + "::/resolutions/20", toy_bins(), pixels(3), {.mode = "a"});
    // Without the mcool magic on the root the file is not multi-res.
    CHECK_FALSE(is_multires_file(mcool));
    CHECK(list_coolers(mcool) ==
          std::vector<std::string>{"/resolutions/10", "/resolutions/20"});
    CHECK_FALSE(is_cooler(mcool));
    CHECK(is_cooler(mcool + "::/resolutions/20"));
    CHECK_THROWS_AS((void)is_cooler(mcool + "::/resolutions/30"), KeyError);

    const std::string text = scratch("not.hdf5");
    { std::FILE* out = std::fopen(text.c_str(), "w"); std::fputs("chrA\t25\n", out); std::fclose(out); }
    CHECK_FALSE(is_cooler(text));
    CHECK_FALSE(is_multires_file(text));
    CHECK_THROWS_AS((void)is_scool_file(text), OSError);
    CHECK_THROWS_AS((void)list_coolers(text), OSError);
    CHECK_THROWS_AS((void)ls(text), OSError);
}

TEST_CASE("cp, mv and ln move groups") {
    const std::string path = write_scool("links.scool");
    cp(path + "::/cells/cell1", path + "::/copied");
    ln(path + "::/cells/cell2", path + "::/hard");
    ln(path + "::/cells/cell2", path + "::/soft", true);
    mv(path + "::/cells/cell10", path + "::/cells/renamed");

    CHECK(is_cooler(path + "::/copied"));
    CHECK(is_cooler(path + "::/hard"));
    CHECK(is_cooler(path + "::/soft"));
    CHECK(Cooler(path + "::/cells/renamed").info().find("nnz")->as_int() == 5);
    CHECK_THROWS_AS(Cooler(path + "::/cells/cell10"), KeyError);
    // The listing resolves links, so the hard and the soft one both appear.
    const std::vector<std::string> listed = list_coolers(path);
    CHECK(std::find(listed.begin(), listed.end(), "/hard") != listed.end());
    CHECK(std::find(listed.begin(), listed.end(), "/soft") != listed.end());

    // Copying the same name twice and copying what is not there.
    CHECK_THROWS_AS(cp(path + "::/cells/cell1", path + "::/copied"), RuntimeError);
    CHECK_THROWS_AS(mv(path + "::/cells/absent", path + "::/gone"), KeyError);

    const std::string other = scratch("other.cool");
    CHECK_THROWS_AS(ln(path + "::/cells/cell1", other + "::/h"), OSError);
    cp(path + "::/cells/cell1", other);
    CHECK(is_cooler(other));
    CHECK(Cooler(other).info().find("nnz")->as_int() == 1);
    const std::string grouped = scratch("grouped.cool");
    cp(path + "::/cells/cell1", grouped + "::/x");
    CHECK(list_coolers(grouped) == std::vector<std::string>{"/x"});
}

TEST_CASE("util: partition, argnatsort, binnify and read_chromsizes") {
    CHECK(partition(0, 9, 2) ==
          std::vector<std::pair<std::int64_t, std::int64_t>>{{0, 2}, {2, 4}, {4, 6}, {6, 8}, {8, 9}});
    CHECK(argnatsort({"chr10", "chr2", "chr1"}) == std::vector<std::int64_t>{2, 1, 0});
    CHECK(natsorted({"chr10", "chr2", "chr1"}) ==
          std::vector<std::string>{"chr1", "chr2", "chr10"});

    const ChromSizes sizes({"chrA", "chrB"}, {25, 17});
    const Table bins = binnify(sizes, 10);
    CHECK(bins.num_rows() == 5);
    CHECK(bins["start"].values<std::int64_t>() == std::vector<std::int64_t>{0, 10, 20, 0, 10});
    CHECK(bins["end"].values<std::int64_t>() == std::vector<std::int64_t>{10, 20, 25, 10, 17});
    CHECK(*bins["chrom"].categorical().categories == std::vector<std::string>{"chrA", "chrB"});
    CHECK(get_binsize(bins).as_int() == 10);
    CHECK_THROWS_AS(binnify(sizes, 0), ZeroDivisionError);

    const std::string path = scratch("sizes.txt");
    {
        std::FILE* out = std::fopen(path.c_str(), "w");
        std::fputs("chr10\t10\nchr2\t20\nchrUn_x\t5\n\nchrM\t3\nchr1\t100\textra\n", out);
        std::fclose(out);
    }
    const ChromSizes filtered = read_chromsizes(path);
    CHECK(filtered.names() == std::vector<std::string>{"chr1", "chr2", "chr10", "chrM"});
    CHECK(filtered.lengths() == std::vector<std::int64_t>{100, 20, 10, 3});
    const ChromSizes all = read_chromsizes(path, {.all_names = true});
    CHECK(all.names() == std::vector<std::string>{"chr10", "chr2", "chrUn_x", "chrM", "chr1"});
}

TEST_CASE("util: rlencode, mad, bedslice and buffered") {
    const Column values(std::vector<std::int32_t>{1, 1, 2, 2, 2, 3});
    const RunLengths runs = rlencode(values);
    CHECK(runs.starts == std::vector<std::int64_t>{0, 2, 5});
    CHECK(runs.lengths == std::vector<std::int64_t>{2, 3, 1});
    CHECK(runs.values.values<std::int32_t>() == std::vector<std::int32_t>{1, 2, 3});
    // Chunking does not change the result; a run spanning two chunks stays one.
    const RunLengths chunked = rlencode(values, 2);
    CHECK(chunked.starts == runs.starts);
    CHECK(chunked.lengths == runs.lengths);

    CHECK(mad(Column{1.0, 2.0, 3.0, 4.0}) == doctest::Approx(1.0));

    const ChromSizes sizes({"chrA", "chrB"}, {25, 17});
    const Table bins = binnify(sizes, 10);
    const Table slice = bedslice(bins, sizes, "chrA:5-20");
    CHECK(slice.num_rows() == 2);
    CHECK(slice["start"].values<std::int64_t>() == std::vector<std::int64_t>{0, 10});
    CHECK(bedslice(bins, sizes, "chrB").num_rows() == 2);
    CHECK_THROWS_AS((void)bedslice(bins, sizes, Region("chrC", std::nullopt, std::nullopt)),
                    ValueError);
    // A sequence of the chromosome table with no bins has no group.
    const ChromSizes extra({"chrA", "chrB", "chrC"}, {25, 17, 9});
    CHECK_THROWS_AS((void)bedslice(bins, extra, "chrC"), KeyError);

    auto chunks = buffered(iterate_chunks({pixels(2), pixels(2), pixels(1)}), 3);
    const std::optional<Table> first = chunks();
    REQUIRE(first.has_value());
    CHECK(first->num_rows() == 4);
    // pd.concat keeps the indexes of the chunks as they are.
    CHECK(first->index().to_vector() == std::vector<std::int64_t>{0, 1, 0, 1});
    const std::optional<Table> second = chunks();
    REQUIRE(second.has_value());
    CHECK(second->num_rows() == 1);
    CHECK_FALSE(chunks().has_value());
}

TEST_CASE("util: GenomeSegmentation and balanced_partition") {
    const ChromSizes sizes({"chrA", "chrB"}, {25, 17});
    const GenomeSegmentation gs(sizes, binnify(sizes, 10));
    CHECK(gs.contigs() == std::vector<std::string>{"chrA", "chrB"});
    CHECK(gs.observed() == std::vector<std::string>{"chrA", "chrB"});
    CHECK(gs.binsize().as_int() == 10);
    CHECK(gs.chrom_binoffset() == std::vector<std::int64_t>{0, 3, 5});
    CHECK(gs.chrom_abspos() == std::vector<std::int64_t>{0, 25, 42});
    CHECK(gs.start_abspos() == std::vector<std::int64_t>{0, 10, 20, 25, 35});
    CHECK(gs.idmap("chrB") == 1);
    CHECK_THROWS_AS((void)gs.idmap("chrC"), KeyError);
    CHECK(gs.fetch("chrA:5-20").num_rows() == 2);

    const std::vector<RegionTuple> granges = balanced_partition(gs, 3, {"chrA", "chrB"});
    // chrB carries fewer bins than chrA, so its step is larger and it stays
    // one range.
    CHECK(granges == std::vector<RegionTuple>{
                         {"chrA", 0, 10}, {"chrA", 10, 20}, {"chrA", 20, 25}, {"chrB", 0, 17}});
    // A contig the file does not carry is skipped.
    CHECK(balanced_partition(gs, 1, {"chrB"}) ==
          std::vector<RegionTuple>{{"chrB", 0, 17}});

    // check_bins turns labels into the categories of the chromosome table and
    // refuses a categorical column with other categories.
    Table labelled = Table{{"chrom", Column(std::vector<std::string>{"chrA", "chrB"})},
                           {"start", Column(std::vector<std::int64_t>{0, 0})},
                           {"end", Column(std::vector<std::int64_t>{25, 17})}};
    CHECK(check_bins(labelled, sizes)["chrom"].categorical().codes ==
          std::vector<std::int32_t>{0, 1});
    CHECK_THROWS_AS((void)check_bins(binnify(ChromSizes({"chrA"}, {25}), 10), sizes), AssertionError);
}
