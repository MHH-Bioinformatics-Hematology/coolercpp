// A single-cell cooler file: writing cells, listing what a file holds, and the
// file-level operations that copy, rename and link groups.
//
//     scool_fileops cells.scool

#include <cstdint>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <coolercpp/coolercpp.hpp>

namespace {

// The bin table the cells share: two chromosomes at 10 bp bins. The chrom
// column is categorical, as cooler stores it: codes into a list of names.
coolercpp::Table makeBins() {
    const std::vector<std::int32_t> codes{0, 0, 0, 1, 1};
    const std::vector<std::string> names{"chrA", "chrB"};
    const std::vector<std::int64_t> start{0, 10, 20, 0, 10};
    const std::vector<std::int64_t> end{10, 20, 30, 10, 20};
    return coolercpp::Table{{"chrom", coolercpp::Column::categorical(codes, names)},
                            {"start", coolercpp::Column{start}},
                            {"end", coolercpp::Column{end}}};
}

// The pixels of one cell: the upper triangle, as bin ids into the table above.
coolercpp::Table makePixels(const std::vector<std::int64_t>& bin1,
                            const std::vector<std::int64_t>& bin2,
                            const std::vector<std::int64_t>& count) {
    return coolercpp::Table{{"bin1_id", coolercpp::Column{bin1}},
                            {"bin2_id", coolercpp::Column{bin2}},
                            {"count", coolercpp::Column{count}}};
}

}  // namespace

int main(int argc, char** argv) {
    const std::string path = argc > 1 ? argv[1] : "cells.scool";

    const coolercpp::Table bins = makeBins();

    // One entry per cell. The names given here become the groups under /cells,
    // sorted; every cell hard links the bin table rather than repeating it.
    const coolercpp::CellPixelTables cells{
        {"cell_1", makePixels({0, 0, 1}, {0, 1, 1}, {5, 3, 7})},
        {"cell_2", makePixels({1, 2, 3}, {1, 2, 4}, {4, 6, 2})},
        {"cell_3", makePixels({0, 3}, {2, 3}, {1, 9})},
    };

    coolercpp::CreateOptions options;
    options.assembly = "toy";
    options.ordered = true;   // the pixels are sorted by (bin1_id, bin2_id)
    coolercpp::create_scool(path, bins, cells, options);

    // What the file holds. list_coolers returns every cooler group in a file,
    // which for a scool file is its cells.
    std::cout << path << ": scool " << std::boolalpha << coolercpp::is_scool_file(path) << ", mcool "
              << coolercpp::is_multires_file(path) << '\n';
    for (const std::string& cell : coolercpp::list_scool_cells(path)) {
        const coolercpp::Cooler clr(path + "::" + cell);
        std::cout << "  " << cell << ": " << clr.info().find("nnz")->as_int() << " pixels, sum "
                  << clr.info().find("sum")->as_int() << '\n';
    }

    // A cell is a cooler group like any other, so the Cooler object reads it.
    const coolercpp::Cooler cell(path + "::/cells/cell_1");
    std::cout << "cell_1 holds " << cell.shape()[0] << " bins over " << cell.chromnames().size()
              << " chromosomes at " << *cell.binsize() << " bp\n";

    // The file-level operations. cp copies a group within one file or into
    // another, mv renames one, and ln links one without copying its data.
    coolercpp::cp(path + "::/cells/cell_1", path + "::/cells/cell_4");
    coolercpp::mv(path + "::/cells/cell_4", path + "::/cells/cell_copy");
    coolercpp::ln(path + "::/cells/cell_2", path + "::/cells/cell_link", /*soft=*/true);
    std::cout << "after cp, mv and ln:";
    for (const std::string& group : coolercpp::list_coolers(path)) {
        std::cout << ' ' << group;
    }
    std::cout << '\n';

    // ls lists the groups and datasets under a URI's group, and
    // pprint_data_tree renders the same tree as cooler's `cooler tree`.
    std::cout << coolercpp::list_coolers(path).size() << " coolers, root holds "
              << coolercpp::ls(path + "::/").size() << " objects\n";
    std::cout << coolercpp::pprint_data_tree(path + "::/cells/cell_1", 1) << '\n';
    return 0;
}
