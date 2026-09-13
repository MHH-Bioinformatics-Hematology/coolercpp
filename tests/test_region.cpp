// cooler.util URI and region parsing.

#include <doctest/doctest.h>

#include <coolercpp/coolercpp.hpp>

using namespace coolercpp;

namespace {

std::string error_of(const std::function<void()>& body) {
    try {
        body();
    } catch (const Error& e) {
        return std::string(e.python_type()) + ": " + e.what();
    }
    return "no error";
}

}  // namespace

TEST_CASE("parse_cooler_uri") {
    CHECK(parse_cooler_uri("a.cool") == std::pair<std::string, std::string>{"a.cool", "/"});
    CHECK(parse_cooler_uri("a.mcool::resolutions/5") ==
          std::pair<std::string, std::string>{"a.mcool", "/resolutions/5"});
    CHECK(error_of([] { (void)parse_cooler_uri("a::b::c"); }) == "ValueError: Invalid Cooler URI string");
}

TEST_CASE("parse_humanized") {
    CHECK(parse_humanized("1,000") == 1000);
    CHECK(parse_humanized("1.5M") == 1500000);
    CHECK(parse_humanized("2kb") == 2000);
    CHECK(parse_humanized("1.1Mb") == 1100000);
    CHECK(error_of([] { (void)parse_humanized("10x"); }) == "ValueError: Unknown unit 'X'");
    CHECK(error_of([] { (void)parse_humanized("1.5"); }) ==
          "ValueError: invalid literal for int() with base 10: '1.5'");
}

TEST_CASE("parse_region_string") {
    CHECK(parse_region_string("chr1:1-2") == GenomicRange{"chr1", 1, 2});
    CHECK(parse_region_string(" chr1 ") == GenomicRange{"chr1", std::nullopt, std::nullopt});
    CHECK(parse_region_string("chr1:5-") == GenomicRange{"chr1", 5, std::nullopt});
    CHECK(parse_region_string("chr1:1-2-3") == GenomicRange{"chr1", 1, 2});
    CHECK(parse_region_string("chr1:  7  -  9  ") == GenomicRange{"chr1", 7, 9});
    CHECK(error_of([] { (void)parse_region_string("chr1:10"); }) ==
          "ValueError: Expected HYPHEN token missing");
    CHECK(error_of([] { (void)parse_region_string(":1-2"); }) ==
          "ValueError: Chromosome name cannot be empty");
    CHECK(error_of([] { (void)parse_region_string("chr1:5-1"); }) ==
          "ValueError: End coordinate less than start");
    CHECK(error_of([] { (void)parse_region_string("chr1:abc"); }) ==
          "ValueError: Unexpected token \"abc\"");
}

TEST_CASE("parse_region") {
    const ChromSizes sizes({"chr1", "chr2"}, {1000, 500}, DType::Int32);
    CHECK(parse_region("chr2", &sizes) == RegionTuple{"chr2", 0, 500});
    CHECK(parse_region(Region("chr1", 10, std::nullopt), &sizes) == RegionTuple{"chr1", 10, 1000});
    CHECK(error_of([&] { (void)parse_region("chr3", &sizes); }) ==
          "ValueError: Unknown sequence label: chr3");
    CHECK(error_of([&] { (void)parse_region("chr2:0-501", &sizes); }) ==
          "ValueError: Genomic region out of bounds: [0, 501)");
    CHECK(error_of([] { (void)parse_region("chr2"); }) ==
          "ValueError: Cannot determine end coordinate.");
    CHECK(error_of([&] { (void)parse_region(Region("chr1", 20, 10), &sizes); }) ==
          "ValueError: End cannot be less than start");
}
