// A downstream program: opens a cooler, prints its shape and one matrix sum.
#include <coolercpp/coolercpp.hpp>

#include <cstdio>

int main(int argc, char** argv) {
    if (argc != 2) {
        std::fprintf(stderr, "usage: consumer FILE.cool\n");
        return 2;
    }
    const coolercpp::Cooler clr(argv[1]);
    const auto shape = clr.shape();
    const auto matrix = clr.matrix({.balance = false}).all();
    long long total = 0;
    for (std::size_t k = 0; k < matrix.dense().values.size(); ++k) {
        total += matrix.dense().values.as_int64(k);
    }
    std::printf("coolercpp %s: %lld x %lld, sum %lld\n", coolercpp::kVersion,
                static_cast<long long>(shape[0]), static_cast<long long>(shape[1]), total);
    return total > 0 ? 0 : 1;
}
