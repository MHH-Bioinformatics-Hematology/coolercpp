// Writing a cool file from tables built in this program, then reading it back.
//
//     create_cool out.cool

#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include <coolercpp/coolercpp.hpp>

namespace {

// The bin table: one row per bin, with the chromosome as a categorical column,
// as cooler stores it. Two chromosomes of 25 bp and 17 bp at a bin size of 10.
coolercpp::Table make_bins() {
    const std::vector<std::int32_t> chrom_codes{0, 0, 0, 1, 1};
    const std::vector<std::string> chrom_labels{"chrA", "chrB"};
    return coolercpp::Table{
        {"chrom", coolercpp::Column::categorical(chrom_codes, chrom_labels)},
        {"start", coolercpp::Column{std::vector<std::int64_t>{0, 10, 20, 0, 10}}},
        {"end", coolercpp::Column{std::vector<std::int64_t>{10, 20, 25, 10, 17}}},
    };
}

// The pixels: the upper triangle, as bin indexes into the table above.
coolercpp::Table make_pixels() {
    return coolercpp::Table{
        {"bin1_id", coolercpp::Column{std::vector<std::int64_t>{0, 0, 1, 2, 3}}},
        {"bin2_id", coolercpp::Column{std::vector<std::int64_t>{0, 1, 1, 3, 4}}},
        {"count", coolercpp::Column{std::vector<std::int32_t>{5, 3, 4, 2, 4}}},
    };
}

}  // namespace

int main(int argc, char** argv) {
    const std::string path = argc > 1 ? argv[1] : "out.cool";

    coolercpp::CreateOptions options;
    options.assembly = "toy";
    options.ordered = true;                  // the pixels are already sorted
    options.dtypes = {{"count", coolercpp::DType::Int32}};

    coolercpp::create_cooler(path, make_bins(), make_pixels(), options);

    const coolercpp::Cooler clr(path);
    std::cout << "wrote " << path << ": " << clr.chromnames().size() << " chromosomes, "
              << clr.shape()[0] << " bins, " << clr.info().find("nnz")->as_int() << " pixels, sum "
              << clr.info().find("sum")->as_int() << '\n';
    return 0;
}
