// The four selectors of a Cooler object: chroms, bins, pixels and matrix, with
// the forms each one returns.
//
//     selectors matrix.cool chrA chrB

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>

#include <coolercpp/coolercpp.hpp>

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "usage: selectors <cool uri> <region> <region2>\n";
        return 2;
    }
    const coolercpp::Cooler clr(argv[1]);
    const coolercpp::Region region{argv[2]};
    const coolercpp::Region region2{argv[3]};

    // chroms and bins return a Table, the counterpart of the DataFrame cooler
    // returns. A Table names its columns and hands out one at a time.
    const coolercpp::Table chroms = clr.chroms().all();
    std::cout << "chroms:";
    for (std::size_t i = 0; i < chroms.num_rows(); ++i) {
        std::cout << ' ' << chroms["name"].label(i) << '(' << chroms["length"].as_int64(i) << ')';
    }
    std::cout << '\n';

    const coolercpp::Table bins = clr.bins().fetch(region);
    std::cout << "bins of " << argv[2] << ": " << bins.num_rows() << ", columns";
    for (const std::string& name : bins.columns()) {
        std::cout << ' ' << name;
    }
    std::cout << '\n';

    // pixels(join = true) replaces the two bin ids with the coordinates of both
    // bins, as cooler's join option does. The 1D selectors take one region; a
    // pair of regions is a matrix query, below.
    const coolercpp::Table pixels = clr.pixels(/*join=*/true).fetch(region);
    std::cout << "pixels of " << argv[2] << ": " << pixels.num_rows();
    if (pixels.num_rows() > 0) {
        std::cout << ", first " << pixels["chrom1"].label(0) << ':'
                  << pixels["start1"].as_int64(0) << " to " << pixels["chrom2"].label(0) << ':'
                  << pixels["start2"].as_int64(0) << " = " << pixels["count"].as_double(0);
    }
    std::cout << '\n';

    // matrix returns a dense array by default, a sparse matrix with
    // sparse = true, and a pixel table with as_pixels = true.
    // The result owns the matrix, so it is held in a variable: a reference
    // bound straight to .dense() of the returned value would dangle.
    const coolercpp::MatrixResult denseResult = clr.matrix({.balance = false}).fetch(region);
    const coolercpp::DenseMatrix& dense = denseResult.dense();
    std::cout << "dense " << dense.rows << " x " << dense.cols << ", first row:";
    for (std::int64_t j = 0; j < std::min<std::int64_t>(dense.cols, 5); ++j) {
        std::cout << ' ' << dense.value(0, j);
    }
    std::cout << '\n';

    const coolercpp::MatrixResult sparseResult =
        clr.matrix({.balance = false, .sparse = true}).fetch(region);
    std::cout << "sparse holds " << sparseResult.sparse().nnz() << " entries\n";

    // as_pixels returns the pixels of the block instead of a matrix, and a
    // second region makes it the block of a chromosome pair.
    const coolercpp::MatrixResult pairResult =
        clr.matrix({.balance = false, .as_pixels = true, .join = true}).fetch(region, region2);
    const coolercpp::Table& pairPixels = pairResult.pixels();
    std::cout << "pixels of the pair " << argv[2] << " x " << argv[3] << ": "
              << pairPixels.num_rows();
    if (pairPixels.num_rows() > 0) {
        std::cout << ", first " << pairPixels["chrom1"].label(0) << ':'
                  << pairPixels["start1"].as_int64(0) << " to "
                  << pairPixels["chrom2"].label(0) << ':' << pairPixels["start2"].as_int64(0)
                  << " = " << pairPixels["count"].as_double(0);
    }
    std::cout << '\n';

    // balance takes a bool or the name of a bin table column. A file converted
    // from .hic carries divisive vectors such as KR, which are applied with
    // divisive_weights = true.
    if (bins.contains("weight")) {
        const coolercpp::MatrixResult balanced = clr.matrix({.balance = true}).fetch(region);
        std::cout << "balanced first value: " << balanced.dense().value(0, 0) << '\n';
    } else {
        std::cout << "no weight column in this file, so nothing to balance\n";
    }
    return 0;
}
