// The shortest program that opens a cool file and reads a region out of it.
//
//     quickstart matrix.mcool::/resolutions/10000 chr1:10M-12M

#include <cstdint>
#include <iostream>
#include <string>

#include <coolercpp/coolercpp.hpp>

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: quickstart <cool uri> <region>\n";
        return 2;
    }
    const coolercpp::Cooler clr(argv[1]);
    const coolercpp::Region region{argv[2]};

    const coolercpp::Table bins = clr.bins().fetch(region);
    // balance = false reads the raw counts. The default, as in cooler, applies
    // the "weight" column, which a file carries only once it has been balanced.
    const coolercpp::MatrixResult matrix =
        clr.matrix({.balance = false, .sparse = true}).fetch(region);
    const coolercpp::SparseMatrix& sparse = matrix.sparse();

    std::cout << argv[2] << " covers " << bins.num_rows() << " bins and holds " << sparse.nnz()
              << " of " << sparse.rows << " x " << sparse.cols << " pixels\n";
    return 0;
}
