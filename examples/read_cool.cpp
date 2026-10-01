// Reading a cool file: metadata, the bin table of a region and a matrix.
//
//     read_cool matrix.cool
//     read_cool matrix.mcool::/resolutions/10000

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <optional>
#include <string>

#include <coolercpp/coolercpp.hpp>

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: read_cool <cool uri>\n";
        return 2;
    }
    const coolercpp::Cooler clr(argv[1]);

    std::cout << "bin size: ";
    if (const std::optional<std::int64_t> binsize = clr.binsize()) {
        std::cout << *binsize << '\n';
    } else {
        std::cout << "variable\n";
    }
    std::cout << "chromosomes: " << clr.chromnames().size() << ", bins: " << clr.shape()[0] << '\n';
    std::cout << "non-zero pixels: " << clr.info().find("nnz")->as_int() << '\n';

    // The bin table of one region. Columns are chrom, start, end, plus any
    // weight or normalization column the file carries.
    const std::string region = clr.chromnames().front();
    const coolercpp::Table bins = clr.bins().fetch(coolercpp::Region{region});
    std::cout << "bins on " << region << ": " << bins.num_rows() << " (";
    for (const std::string& name : bins.columns()) {
        std::cout << name << ' ';
    }
    std::cout << ")\n";

    // The matrix of that region, as raw counts in a dense array.
    const coolercpp::MatrixResult matrix = clr.matrix({.balance = false}).fetch(coolercpp::Region{region});
    const coolercpp::DenseMatrix& dense = matrix.dense();
    std::cout << "matrix " << dense.rows << " x " << dense.cols << ", first row:";
    for (std::int64_t j = 0; j < std::min<std::int64_t>(dense.cols, 5); ++j) {
        std::cout << ' ' << dense.value(0, j);
    }
    std::cout << '\n';
    return 0;
}
