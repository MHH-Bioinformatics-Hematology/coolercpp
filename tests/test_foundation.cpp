// dtype promotion, the columnar table, JSON and the numpy compatibility helpers.

#include <doctest/doctest.h>

#include <coolercpp/coolercpp.hpp>

#include <cmath>
#include <limits>

#include "numpy_compat.hpp"

using namespace coolercpp;

TEST_CASE("numpy 1.26 array promotion") {
    CHECK(result_type(DType::Int32, DType::Float32) == DType::Float64);
    CHECK(result_type(DType::Int16, DType::Float32) == DType::Float32);
    CHECK(result_type(DType::UInt8, DType::Int8) == DType::Int16);
    CHECK(result_type(DType::UInt32, DType::Int32) == DType::Int64);
    CHECK(result_type(DType::UInt64, DType::Int64) == DType::Float64);
    CHECK(result_type(DType::Bool, DType::Int32) == DType::Int32);
    CHECK(result_type(DType::Float32, DType::Float32) == DType::Float32);
    CHECK(true_divide_result(DType::Int32) == DType::Float64);
    CHECK(true_divide_result(DType::Float32) == DType::Float32);
    CHECK_THROWS_AS((void)result_type(DType::String, DType::Int32), TypeError);
}

TEST_CASE("dtype names") {
    CHECK(dtype_from_name("i4") == DType::Int32);
    CHECK(dtype_from_name("<f8") == DType::Float64);
    CHECK(dtype_from_name("float") == DType::Float64);
    CHECK(dtype_from_name("int") == DType::Int64);
    CHECK(dtype_name(DType::Categorical) == "category");
    CHECK_THROWS_AS((void)dtype_from_name("complex128"), TypeError);
}

TEST_CASE("columns convert like numpy astype and keep categories on slicing") {
    const Column floats{1.7, -2.6, 3.0};
    CHECK(floats.as<std::int32_t>() == std::vector<std::int32_t>{1, -2, 3});
    CHECK(floats.astype(DType::Bool).values<Bool8>()[0]);

    Column cats = Column::categorical({2, 0, 1}, {"a", "b", "c"});
    CHECK(cats.label(0) == "c");
    const Column part = cats.slice(1, 3);
    CHECK(part.categorical().categories == cats.categorical().categories);
    cats.append(part);
    CHECK(cats.size() == 5);
    CHECK_THROWS_AS(cats.append(Column::categorical({0}, {"z"})), ValueError);
    CHECK_THROWS_AS(cats.append(floats), TypeError);

    const std::vector<std::int64_t> positions{2, 0};
    CHECK(floats.take(positions).values<double>() == std::vector<double>{3.0, 1.7});
    CHECK_THROWS_AS((void)floats.values<float>(), TypeError);
}

TEST_CASE("tables keep column order, their index and pandas' errors") {
    Table t{{"b", Column{1, 2, 3}}, {"a", Column{0.5, 1.5, 2.5}}};
    CHECK(t.columns() == std::vector<std::string>{"b", "a"});
    CHECK(t.num_rows() == 3);
    CHECK_THROWS_AS(t.set("c", Column{1, 2}), ValueError);
    CHECK_THROWS_AS((void)t.select({"missing"}), KeyError);
    t.set("b", Column{7, 8, 9});
    CHECK(t.columns().front() == "b");

    const Table tail = t.slice(1, 3);
    CHECK(tail.index().start() == 1);
    const std::vector<std::int64_t> order{2, 0};
    const Table picked = t.take(order);
    CHECK(picked.index()[0] == 2);
    CHECK(picked["a"].values<double>()[1] == 0.5);
}

TEST_CASE("simplejson.loads semantics") {
    CHECK(json::parse("null")->is_null());
    CHECK_FALSE(json::parse("HDF5::Cooler").has_value());
    CHECK_FALSE(json::parse("2018-04-26T16:00:58.951917").has_value());
    CHECK(json::parse("123")->as_int() == 123);
    CHECK(json::parse(" 1.5e3 ")->as_double() == 1500.0);
    CHECK_FALSE(json::parse("01").has_value());
    CHECK_FALSE(json::parse("[1,]").has_value());
    CHECK(std::isnan(json::parse("NaN")->as_double()));
    CHECK(json::parse("-Infinity")->as_double() == -std::numeric_limits<double>::infinity());
    CHECK(json::parse("\"\\u00e9\"")->as_string() == "\xc3\xa9");
    CHECK_FALSE(json::parse("\"tab\there\"").has_value());

    const auto object = json::parse(R"({"b": 1, "a": 2, "b": 3})");
    REQUIRE(object.has_value());
    REQUIRE(object->as_object().size() == 2);
    CHECK(object->as_object()[0].first == "b");
    CHECK(object->as_object()[0].second.as_int() == 3);
}

TEST_CASE("simplejson.dumps output") {
    json::Value value = json::Value::object();
    value["a"] = json::Value::array({json::Value(1.5), json::Value(), json::Value(true)});
    value["\xc3\xa9"] = "x\"y";
    CHECK(json::dumps(value) == R"({"a": [1.5, null, true], "\u00e9": "x\"y"})");
    CHECK(json::dumps(json::Value(1e16)) == "1e+16");
    CHECK(json::dumps(json::Value(0.1)) == "0.1");
    CHECK(json::dumps(json::Value(std::nan(""))) == "NaN");
    CHECK(json::dumps(json::Value::object()) == "{}");
}

TEST_CASE("numpy reductions and linspace") {
    std::vector<double> x(10000);
    for (std::size_t i = 0; i < x.size(); ++i) {
        x[i] = static_cast<double>(i) * 0.1;
    }
    CHECK(npy::pairwise_sum(x.data(), x.size()) == 4999500.0);
    // np.sum((np.sin(np.arange(20000)) * 1e3).astype(np.float32))
    std::vector<float> y(20000);
    for (std::size_t i = 0; i < y.size(); ++i) {
        y[i] = static_cast<float>(std::sin(static_cast<double>(i)) * 1e3);
    }
    CHECK(npy::pairwise_sum(y.data(), y.size()) == -120.02685546875F);

    CHECK(npy::linspace_int(0, 24108, 3) == std::vector<std::int64_t>{0, 12054, 24108});
    CHECK(npy::linspace_int(7, 1000, 13) ==
          std::vector<std::int64_t>{7, 89, 172, 255, 338, 420, 503, 586, 669, 751, 834, 917, 1000});
    CHECK(npy::linspace_int(0, 5, 1) == std::vector<std::int64_t>{0});
}

TEST_CASE("Python slice and index resolution") {
    CHECK(npy::slice_indices(-100000, std::nullopt, 10).start == 0);
    CHECK(npy::slice_indices(3, 100, 10).stop == 10);
    CHECK(npy::slice_indices(8, 2, 10).size() == 0);
    CHECK(npy::wrap_index(-1, 10) == 9);
    CHECK_THROWS_AS((void)npy::wrap_index(10, 10), IndexError);
    CHECK(npy::float_repr(1e-5) == "1e-05");
    CHECK(npy::float_repr(123456789.0) == "123456789.0");
    CHECK(npy::float32_repr(0.1F) == "0.1");
}
