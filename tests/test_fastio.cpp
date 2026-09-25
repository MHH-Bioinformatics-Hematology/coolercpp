// The threaded chunk writer and reader and the in-place file editing: round
// trips, files that do not depend on the number of threads, chunk shapes.

#include <doctest/doctest.h>

#include <coolercpp/fastio.hpp>

#include <sys/stat.h>

#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <numeric>
#include <string>
#include <vector>

using namespace coolercpp;

namespace {

std::string scratch(const std::string& name) {
    ::mkdir(COOLERCPP_TEST_SCRATCH, 0755);
    const std::string path = std::string(COOLERCPP_TEST_SCRATCH) + "/" + name;
    std::remove(path.c_str());
    return path;
}

std::string slurp(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}

// Writes a pixel-like triple of columns with `threads` threads.
void write_pixels(const std::string& path, int threads, std::size_t rows) {
    ThreadPool pool(threads);
    fastio::File file(path, fastio::Mode::Create);
    file.create_group("/pixels");
    fastio::ChunkedWriter bin1(file.create_dataset("/pixels/bin1_id", fastio::ColumnType::int64(), 0, true), 8,
                               true);
    fastio::ChunkedWriter count(file.create_dataset("/pixels/count", fastio::ColumnType::int32(), 0, true), 4,
                                true);
    std::vector<std::int64_t> ids(rows);
    std::vector<std::int32_t> counts(rows);
    for (std::size_t i = 0; i < rows; ++i) {
        ids[i] = static_cast<std::int64_t>(i / 3);
        counts[i] = static_cast<std::int32_t>(i % 97);
    }
    const std::size_t step = 100000;
    for (std::size_t lo = 0; lo < rows; lo += step) {
        const std::size_t n = std::min(step, rows - lo);
        bin1.append(ids.data() + lo, n);
        count.append(counts.data() + lo, n);
        fastio::flush(pool, {&bin1, &count}, false);
    }
    fastio::flush(pool, {&bin1, &count}, true);
    file.set_attribute("/", "nnz", std::int64_t(rows));
    file.set_attribute("/", "format", std::string("HDF5::Cooler"));
}

}  // namespace

TEST_CASE("fastio: written columns read back through the file and the chunk reader") {
    const std::string path = scratch("fastio_roundtrip.h5");
    const std::size_t rows = 350001;
    write_pixels(path, 4, rows);

    fastio::File file(path, fastio::Mode::Read);
    CHECK(file.length("/pixels/bin1_id") == rows);
    const auto ids = file.read_int64("/pixels/bin1_id");
    REQUIRE(ids.size() == rows);
    CHECK(ids[0] == 0);
    CHECK(ids[rows - 1] == static_cast<std::int64_t>((rows - 1) / 3));

    ThreadPool pool(4);
    fastio::ColumnReader reader(file, "/pixels/count");
    REQUIRE(reader.size() == rows);
    std::vector<std::int64_t> counts(200000);
    reader.read(pool, 12345, 12345 + counts.size(), counts.data());
    for (std::size_t i = 0; i < counts.size(); ++i) {
        REQUIRE(counts[i] == static_cast<std::int64_t>((12345 + i) % 97));
    }
    const auto nnz = file.attribute("/", "nnz");
    REQUIRE(nnz.has_value());
    CHECK(std::get<std::int64_t>(*nnz) == static_cast<std::int64_t>(rows));
}

TEST_CASE("fastio: the file does not depend on the number of threads") {
    const std::string one = scratch("fastio_t1.h5");
    const std::string eight = scratch("fastio_t8.h5");
    write_pixels(one, 1, 420000);
    write_pixels(eight, 8, 420000);
    CHECK(slurp(one) == slurp(eight));
}

TEST_CASE("fastio: fixed string, enum and float columns, attributes and unlinking") {
    const std::string path = scratch("fastio_types.h5");
    ThreadPool pool(2);
    {
        fastio::File file(path, fastio::Mode::Create);
        file.create_group("/chroms");
        const char names[] = {'c', 'h', 'r', 'A', 0, 0, 'c', 'h', 'r', 'B', 0, 0};
        fastio::write_dataset(pool, file, "/chroms/name", fastio::ColumnType::fixed_string(6), names, 2);
        const std::vector<std::int32_t> chrom = {0, 0, 1};
        fastio::write_dataset(pool, file, "/chroms/chrom", fastio::ColumnType::enumeration({"chrA", "chrB"}),
                              chrom.data(), chrom.size());
        const std::vector<double> weight = {1.5, 0.25, 2.0};
        fastio::write_dataset(pool, file, "/chroms/weight", fastio::ColumnType::float64(), weight.data(),
                              weight.size());
    }
    fastio::File file(path, fastio::Mode::ReadWrite);
    CHECK(file.is_dataset("/chroms/name"));
    CHECK(file.read_strings("/chroms/name") == std::vector<std::string>{"chrA", "chrB"});
    CHECK(file.read_doubles("/chroms/weight") == std::vector<double>{1.5, 0.25, 2.0});
    file.write_doubles("/chroms/weight", {3.0, 4.0, 5.0});
    CHECK(file.read_doubles("/chroms/weight") == std::vector<double>{3.0, 4.0, 5.0});
    file.set_attribute("/chroms", "note", std::string("x"));
    file.set_attribute("/chroms", "note", std::string("y"));
    CHECK(std::get<std::string>(*file.attribute("/chroms", "note")) == "y");
    file.unlink("/chroms/weight");
    CHECK_FALSE(file.exists("/chroms/weight"));
}
