// coolercpp-harness: runs one harness case through coolercpp and writes the
// outputs in the neutral format compare.py reads (result.json plus .npy
// files). The Python counterpart is oracle.py; both take the same case file.
//
//   coolercpp-harness CASE.json OUT_DIR
//
// Library errors are results (kind "error"), not failures of the driver. The
// exit status is non-zero only when result.json could not be written.

#include <coolercpp/coolercpp.hpp>

#include <sys/stat.h>

#include <chrono>
#include <cstdio>
#include <fstream>
#include <functional>
#include <sstream>
#include <stdexcept>

#include "npy.hpp"

using coolercpp::AxisKey;
using coolercpp::Column;
using coolercpp::Cooler;
using coolercpp::DType;
using coolercpp::Fields;
using coolercpp::MatrixResult;
using coolercpp::RangeSelector1D;
using coolercpp::Region;
using coolercpp::Slice;
using coolercpp::Table;
using coolercpp::json::Value;

namespace {

std::string read_text(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("cannot read " + path);
    }
    std::stringstream buffer;
    buffer << in.rdbuf();
    return buffer.str();
}

void write_text(const std::string& path, const std::string& text) {
    std::ofstream out(path, std::ios::binary);
    out << text;
    if (!out) {
        throw std::runtime_error("cannot write " + path);
    }
}

const Value& member(const Value& v, const char* key) {
    const Value* found = v.find(key);
    if (found == nullptr) {
        throw std::runtime_error(std::string("case lacks '") + key + "'");
    }
    return *found;
}

bool flag(const Value& v, const char* key, bool fallback) {
    const Value* found = v.find(key);
    return (found == nullptr || found->is_null()) ? fallback : found->as_bool();
}

std::optional<std::int64_t> optional_int(const Value& v) {
    if (v.is_null()) {
        return std::nullopt;
    }
    return v.as_int();
}

Value error_value(const char* type, const std::string& message) {
    Value out = Value::object();
    out["kind"] = "error";
    out["type"] = type;
    out["message"] = message;
    return out;
}

Value capture(const std::function<Value()>& body) {
    try {
        return body();
    } catch (const coolercpp::Error& e) {
        return error_value(e.python_type(), e.what());
    }
}

Value value_result(Value v) {
    Value out = Value::object();
    out["kind"] = "value";
    out["value"] = std::move(v);
    return out;
}

// A JSON value with the type information of the Python object it stands for.
Value tagged(const Value& v) {
    Value out = Value::object();
    switch (v.type()) {
        case coolercpp::json::Type::Null: out["t"] = "null"; break;
        case coolercpp::json::Type::Bool:
            out["t"] = "bool";
            out["v"] = v.as_bool();
            break;
        case coolercpp::json::Type::Int:
            out["t"] = "int";
            out["dtype"] = std::string(coolercpp::dtype_name(v.dtype()));
            out["v"] = v.as_int();
            break;
        case coolercpp::json::Type::UInt:
            out["t"] = "int";
            out["dtype"] = std::string(coolercpp::dtype_name(v.dtype()));
            out["v"] = v.as_uint();
            break;
        case coolercpp::json::Type::Double:
            out["t"] = "float";
            out["dtype"] = std::string(coolercpp::dtype_name(v.dtype()));
            out["v"] = v.as_double();
            break;
        case coolercpp::json::Type::String:
            out["t"] = "str";
            out["v"] = v.as_string();
            break;
        case coolercpp::json::Type::Bytes:
            out["t"] = "bytes";
            out["v"] = v.as_string();
            break;
        case coolercpp::json::Type::Array: {
            out["t"] = "list";
            Value items = Value::array({});
            for (const Value& item : v.as_array()) {
                items.push_back(tagged(item));
            }
            out["v"] = items;
            break;
        }
        case coolercpp::json::Type::Object: {
            out["t"] = "dict";
            Value items = Value::array({});
            for (const auto& [key, item] : v.as_object()) {
                items.push_back(Value::array({Value(key), tagged(item)}));
            }
            out["v"] = items;
            break;
        }
    }
    return out;
}

Region region_from(const Value& v) {
    if (v.is_string()) {
        return Region(v.as_string());
    }
    const auto& items = v.as_array();
    return Region(items.at(0).as_string(), optional_int(items.at(1)), optional_int(items.at(2)));
}

Value string_array(const std::vector<std::string>& items) {
    Value out = Value::array({});
    for (const std::string& item : items) {
        out.push_back(Value(item));
    }
    return out;
}

// ---- neutral output ----

class Output {
  public:
    explicit Output(std::string dir) : dir_(std::move(dir)) {}

    std::string next_file() { return "a" + std::to_string(counter_++) + ".npy"; }

    Value table(const Table& t) {
        Value columns = Value::array({});
        for (std::size_t c = 0; c < t.num_columns(); ++c) {
            const Column& column = t.column(c);
            Value entry = Value::object();
            entry["name"] = t.columns()[c];
            entry["dtype"] = std::string(coolercpp::dtype_name(column.dtype()));
            if (column.dtype() == DType::String) {
                entry["values"] = string_array(column.values<std::string>());
            } else if (column.dtype() == DType::Categorical) {
                const auto& data = column.categorical();
                const std::string file = next_file();
                harness::write_npy(dir_ + "/" + file, Column(data.codes));
                entry["codes"] = file;
                entry["categories"] = string_array(*data.categories);
                entry["ordered"] = data.ordered;
            } else {
                const std::string file = next_file();
                harness::write_npy(dir_ + "/" + file, column);
                entry["file"] = file;
            }
            columns.push_back(entry);
        }
        Value out = Value::object();
        out["kind"] = "table";
        out["columns"] = columns;
        Value index = Value::object();
        if (t.index().is_range()) {
            index["kind"] = "range";
            index["start"] = t.index().start();
            index["size"] = static_cast<std::int64_t>(t.index().size());
        } else {
            const std::string file = next_file();
            harness::write_npy(dir_ + "/" + file, Column(t.index().to_vector()));
            index["kind"] = "values";
            index["file"] = file;
        }
        out["index"] = index;
        return out;
    }

    Value matrix(const MatrixResult& result) {
        if (result.is_pixels()) {
            return table(result.pixels());
        }
        Value out = Value::object();
        if (result.is_dense()) {
            const auto& dense = result.dense();
            const std::string file = next_file();
            harness::write_npy(dir_ + "/" + file, dense.values, {dense.rows, dense.cols});
            out["kind"] = "dense";
            out["shape"] = Value::array({Value(dense.rows), Value(dense.cols)});
            out["dtype"] = std::string(coolercpp::dtype_name(dense.dtype()));
            out["file"] = file;
            return out;
        }
        const auto& sparse = result.sparse();
        out["kind"] = "sparse";
        out["shape"] = Value::array({Value(sparse.rows), Value(sparse.cols)});
        for (const auto& [name, column] :
             {std::pair<const char*, const Column*>{"row", &sparse.row},
              {"col", &sparse.col},
              {"data", &sparse.data}}) {
            const std::string file = next_file();
            harness::write_npy(dir_ + "/" + file, *column);
            Value entry = Value::object();
            entry["file"] = file;
            entry["dtype"] = std::string(coolercpp::dtype_name(column->dtype()));
            out[name] = entry;
        }
        return out;
    }

  private:
    std::string dir_;
    int counter_ = 0;
};

// ---- neutral input ----

Table read_table(const std::string& dir) {
    const std::optional<Value> meta = coolercpp::json::parse(read_text(dir + "/meta.json"));
    if (!meta.has_value()) {
        throw std::runtime_error("bad table meta in " + dir);
    }
    Table table;
    for (const Value& entry : member(*meta, "columns").as_array()) {
        const std::string name = member(entry, "name").as_string();
        const std::string dtype = member(entry, "dtype").as_string();
        if (dtype == "object") {
            std::vector<std::string> values;
            for (const Value& item : member(entry, "values").as_array()) {
                values.push_back(item.as_string());
            }
            table.set(name, Column(std::move(values)));
        } else if (dtype == "category") {
            const Column codes = harness::read_npy(dir + "/" + member(entry, "codes").as_string());
            std::vector<std::string> categories;
            for (const Value& item : member(entry, "categories").as_array()) {
                categories.push_back(item.as_string());
            }
            table.set(name, Column::categorical(codes.as<std::int32_t>(), std::move(categories),
                                                flag(entry, "ordered", true)));
        } else {
            table.set(name, harness::read_npy(dir + "/" + member(entry, "file").as_string()));
        }
    }
    return table;
}

coolercpp::H5OptValue h5opt_value(const Value& v) {
    switch (v.type()) {
        case coolercpp::json::Type::Null: return std::monostate{};
        case coolercpp::json::Type::Bool: return v.as_bool();
        case coolercpp::json::Type::Int:
        case coolercpp::json::Type::UInt: return v.as_int();
        case coolercpp::json::Type::Double: return v.as_double();
        case coolercpp::json::Type::String: return v.as_string();
        case coolercpp::json::Type::Array: {
            std::vector<std::int64_t> items;
            for (const Value& item : v.as_array()) {
                items.push_back(item.is_null() ? -1 : item.as_int());
            }
            return items;
        }
        default: throw std::runtime_error("unsupported h5opts value");
    }
}

coolercpp::CreateOptions create_options(const Value& o) {
    coolercpp::CreateOptions options;
    if (const Value* v = o.find("columns"); v != nullptr && !v->is_null()) {
        options.columns.emplace();
        for (const Value& item : v->as_array()) {
            options.columns->push_back(item.as_string());
        }
    }
    if (const Value* v = o.find("dtypes"); v != nullptr && !v->is_null()) {
        for (const auto& [key, item] : v->as_object()) {
            options.dtypes.emplace_back(key, coolercpp::dtype_from_name(item.as_string()));
        }
    }
    if (const Value* v = o.find("metadata"); v != nullptr && !v->is_null()) {
        options.metadata = *v;
    }
    if (const Value* v = o.find("assembly"); v != nullptr && !v->is_null()) {
        options.assembly = v->as_string();
    }
    options.ordered = flag(o, "ordered", false);
    options.symmetric_upper = flag(o, "symmetric_upper", true);
    if (const Value* v = o.find("mode"); v != nullptr && !v->is_null()) {
        options.mode = v->as_string();
    }
    if (const Value* v = o.find("mergebuf"); v != nullptr && !v->is_null()) {
        options.mergebuf = v->as_int();
    }
    if (const Value* v = o.find("max_merge"); v != nullptr && !v->is_null()) {
        options.max_merge = v->as_int();
    }
    if (const Value* v = o.find("temp_dir"); v != nullptr && !v->is_null()) {
        options.temp_dir = v->as_string();
    }
    options.delete_temp = flag(o, "delete_temp", true);
    options.boundscheck = flag(o, "boundscheck", true);
    options.dupcheck = flag(o, "dupcheck", true);
    options.triucheck = flag(o, "triucheck", true);
    options.ensure_sorted = flag(o, "ensure_sorted", false);
    if (const Value* v = o.find("h5opts"); v != nullptr && !v->is_null()) {
        options.h5opts.emplace();
        for (const auto& [key, item] : v->as_object()) {
            options.h5opts->emplace_back(key, h5opt_value(item));
        }
    }
    return options;
}

coolercpp::MatrixOptions matrix_options(const Value& o) {
    coolercpp::MatrixOptions options;
    if (const Value* v = o.find("field"); v != nullptr && !v->is_null()) {
        options.field = v->as_string();
    }
    if (const Value* v = o.find("balance"); v != nullptr && !v->is_null()) {
        options.balance = v->is_bool() ? coolercpp::Balance(v->as_bool())
                                       : coolercpp::Balance(v->as_string());
    }
    options.sparse = flag(o, "sparse", false);
    options.as_pixels = flag(o, "as_pixels", false);
    options.join = flag(o, "join", false);
    options.ignore_index = flag(o, "ignore_index", true);
    if (const Value* v = o.find("divisive_weights"); v != nullptr && !v->is_null()) {
        options.divisive_weights = v->as_bool();
    }
    if (const Value* v = o.find("chunksize"); v != nullptr && !v->is_null()) {
        options.chunksize = v->as_int();
    }
    return options;
}

Slice slice_from(const Value& v) {
    const auto& items = v.as_array();
    return Slice{optional_int(items.at(0)), optional_int(items.at(1)),
                 items.size() > 2 ? optional_int(items[2]) : std::nullopt};
}

AxisKey axis_from(const Value& v) {
    if (v.is_array()) {
        return slice_from(v);
    }
    return AxisKey(v.as_int());
}

// ---- operations ----

Value op_info(const Value& spec) {
    const Cooler c(member(spec, "uri").as_string());
    Value items = Value::object();
    items["info"] = capture([&] { return value_result(tagged(c.info())); });
    items["chromnames"] = capture([&] { return value_result(string_array(c.chromnames())); });
    items["chromsizes"] = capture([&] {
        Value v = Value::object();
        v["names"] = string_array(c.chromsizes().names());
        Value lengths = Value::array({});
        for (const std::int64_t l : c.chromsizes().lengths()) {
            lengths.push_back(Value(l));
        }
        v["lengths"] = lengths;
        v["dtype"] = std::string(coolercpp::dtype_name(c.chromsizes().dtype()));
        return value_result(v);
    });
    items["binsize"] = capture([&] {
        const auto bs = c.binsize();
        return value_result(bs.has_value() ? Value(*bs) : Value());
    });
    items["storage_mode"] = capture([&] { return value_result(Value(c.storage_mode())); });
    items["shape"] = capture([&] {
        const auto s = c.shape();
        return value_result(Value::array({Value(s[0]), Value(s[1])}));
    });
    items["root"] = value_result(Value(c.root()));
    Value out = Value::object();
    out["kind"] = "multi";
    out["items"] = items;
    return out;
}

Value op_extent(const Value& spec) {
    const Cooler c(member(spec, "uri").as_string());
    const bool offset = member(spec, "op").as_string() == "offset";
    Value items = Value::array({});
    for (const Value& region : member(spec, "regions").as_array()) {
        items.push_back(capture([&] {
            if (offset) {
                return value_result(Value(c.offset(region_from(region))));
            }
            const auto [a, b] = c.extent(region_from(region));
            return value_result(Value::array({Value(a), Value(b)}));
        }));
    }
    Value out = Value::object();
    out["kind"] = "list";
    out["items"] = items;
    return out;
}

Value op_table(const Value& spec, Output& output) {
    const Cooler c(member(spec, "uri").as_string());
    const std::string which = member(spec, "selector").as_string();
    const Value kwargs = spec.find("selector_kwargs") ? member(spec, "selector_kwargs") : Value::object();
    const bool convert_enum = flag(kwargs, "convert_enum", true);
    RangeSelector1D sel = which == "chroms" ? c.chroms(convert_enum)
                          : which == "bins" ? c.bins(convert_enum)
                                            : c.pixels(flag(kwargs, "join", false), convert_enum);
    if (const Value* fields = spec.find("fields"); fields != nullptr && !fields->is_null()) {
        if (fields->is_string()) {
            sel = sel[Fields(fields->as_string())];
        } else {
            std::vector<std::string> names;
            for (const Value& item : fields->as_array()) {
                names.push_back(item.as_string());
            }
            sel = sel[Fields(names)];
        }
    }
    const Value& access = member(spec, "access");
    if (const Value* v = access.find("slice")) {
        return output.table(sel[slice_from(*v)]);
    }
    if (const Value* v = access.find("row")) {
        return output.table(sel[v->as_int()]);
    }
    if (const Value* v = access.find("fetch")) {
        return output.table(sel.fetch(region_from(*v)));
    }
    if (access.find("columns") != nullptr) {
        return value_result(string_array(sel.columns()));
    }
    if (access.find("dtypes") != nullptr) {
        Value list = Value::array({});
        for (const auto& [name, dtype] : sel.dtypes()) {
            list.push_back(Value::array({Value(name), Value(std::string(coolercpp::dtype_name(dtype)))}));
        }
        return value_result(list);
    }
    if (access.find("len") != nullptr) {
        return value_result(Value(static_cast<std::int64_t>(sel.size())));
    }
    throw std::runtime_error("unknown table access");
}

Value op_matrix(const Value& spec, Output& output) {
    const Cooler c(member(spec, "uri").as_string());
    const Value options = spec.find("options") ? member(spec, "options") : Value::object();
    const auto sel = c.matrix(matrix_options(options));
    const Value& access = member(spec, "access");
    if (const Value* v = access.find("slice")) {
        const auto& axes = v->as_array();
        if (axes.size() == 1) {
            return output.matrix(sel[axis_from(axes[0])]);
        }
        return output.matrix(sel(axis_from(axes.at(0)), axis_from(axes.at(1))));
    }
    if (const Value* v = access.find("fetch")) {
        const auto& regions = v->as_array();
        std::optional<Region> second;
        if (regions.size() > 1 && !regions[1].is_null()) {
            second = region_from(regions[1]);
        }
        return output.matrix(sel.fetch(region_from(regions.at(0)), second));
    }
    throw std::runtime_error("unknown matrix access");
}

Value op_parse_region(const Value& spec) {
    std::optional<coolercpp::ChromSizes> sizes;
    if (const Value* uri = spec.find("chromsizes_from"); uri != nullptr && !uri->is_null()) {
        sizes = Cooler(uri->as_string()).chromsizes();
    }
    Value items = Value::array({});
    for (const Value& region : member(spec, "regions").as_array()) {
        items.push_back(capture([&] {
            const auto r = coolercpp::parse_region(region_from(region), sizes ? &*sizes : nullptr);
            return value_result(Value::array({Value(r.chrom), Value(r.start), Value(r.end)}));
        }));
    }
    Value out = Value::object();
    out["kind"] = "list";
    out["items"] = items;
    return out;
}

Value op_dump(const Value& spec, Output& output) {
    const Cooler c(member(spec, "uri").as_string());
    Value items = Value::object();
    items["info"] = capture([&] { return value_result(tagged(c.info())); });
    items["chroms"] = capture([&] { return output.table(c.chroms().all()); });
    items["bins"] = capture([&] { return output.table(c.bins().all()); });
    items["pixels"] = capture([&] { return output.table(c.pixels().all()); });
    Value out = Value::object();
    out["kind"] = "multi";
    out["items"] = items;
    return out;
}

Value op_range_query(const Value& spec, Output& output) {
    const Cooler c(member(spec, "uri").as_string());
    const auto& box = member(spec, "bbox").as_array();
    const coolercpp::RangeQuery2D query(
        c,
        member(spec, "kind").as_string() == "direct" ? coolercpp::RangeQuery2D::Kind::Direct
                                                     : coolercpp::RangeQuery2D::Kind::FillLower,
        member(spec, "field").as_string(),
        {box.at(0).as_int(), box.at(1).as_int(), box.at(2).as_int(), box.at(3).as_int()},
        member(spec, "chunksize").as_int(), flag(spec, "return_index", false));
    const Value& access = member(spec, "access");
    if (access.find("n_chunks") != nullptr) {
        return value_result(Value(static_cast<std::int64_t>(query.n_chunks())));
    }
    if (const Value* v = access.find("chunk")) {
        return output.table(query.get_chunk(static_cast<std::size_t>(v->as_int())));
    }
    if (access.find("chunks") != nullptr) {
        Value items = Value::array({});
        for (std::size_t i = 0; i < query.n_chunks(); ++i) {
            items.push_back(output.table(query.get_chunk(i)));
        }
        Value out = Value::object();
        out["kind"] = "list";
        out["items"] = items;
        return out;
    }
    if (access.find("frame") != nullptr) {
        return output.table(query.to_frame());
    }
    if (access.find("sparse") != nullptr) {
        return output.matrix(MatrixResult(query.to_sparse_matrix()));
    }
    if (access.find("array") != nullptr) {
        return output.matrix(MatrixResult(query.to_array()));
    }
    throw std::runtime_error("unknown range query access");
}

Value op_is_cooler(const Value& spec) {
    Value items = Value::array({});
    for (const Value& uri : member(spec, "uris").as_array()) {
        items.push_back(capture([&] { return value_result(Value(coolercpp::is_cooler(uri.as_string()))); }));
    }
    Value out = Value::object();
    out["kind"] = "list";
    out["items"] = items;
    return out;
}

Value op_create(const Value& spec) {
    for (const Value& step : member(spec, "steps").as_array()) {
        const Table bins = read_table(member(step, "bins").as_string());
        const coolercpp::CreateOptions options =
            create_options(step.find("options") ? member(step, "options") : Value::object());
        const std::string uri = member(step, "uri").as_string();
        if (const Value* pixels = step.find("pixels")) {
            const Table table = read_table(pixels->as_string());
            coolercpp::create_cooler(uri, bins, table, options);
        } else {
            std::vector<std::string> dirs;
            for (const Value& item : member(step, "pixel_chunks").as_array()) {
                dirs.push_back(item.as_string());
            }
            auto position = std::make_shared<std::size_t>(0);
            coolercpp::PixelChunks chunks = [dirs, position]() -> std::optional<Table> {
                if (*position >= dirs.size()) {
                    return std::nullopt;
                }
                return read_table(dirs[(*position)++]);
            };
            coolercpp::create_cooler(uri, bins, chunks, options);
        }
    }
    Value out = Value::object();
    out["kind"] = "created";
    return out;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::fprintf(stderr, "usage: coolercpp-harness CASE.json OUT_DIR\n");
        return 2;
    }
    const std::string out_dir = argv[2];
    ::mkdir(out_dir.c_str(), 0755);
    Value result;
    const auto start = std::chrono::steady_clock::now();
    try {
        const std::optional<Value> spec = coolercpp::json::parse(read_text(argv[1]));
        if (!spec.has_value()) {
            throw std::runtime_error("case file is not valid JSON");
        }
        const std::string op = member(*spec, "op").as_string();
        Output output(out_dir);
        result = capture([&]() -> Value {
            if (op == "info") return op_info(*spec);
            if (op == "extent" || op == "offset") return op_extent(*spec);
            if (op == "table") return op_table(*spec, output);
            if (op == "matrix") return op_matrix(*spec, output);
            if (op == "parse_region") return op_parse_region(*spec);
            if (op == "list_coolers") {
                return value_result(string_array(coolercpp::list_coolers(member(*spec, "path").as_string())));
            }
            if (op == "dump") return op_dump(*spec, output);
            if (op == "create") return op_create(*spec);
            if (op == "range_query") return op_range_query(*spec, output);
            if (op == "is_cooler") return op_is_cooler(*spec);
            throw std::runtime_error("unknown op " + op);
        });
    } catch (const std::exception& e) {
        result = error_value("HarnessError", e.what());
    }
    const auto seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    Value document = Value::object();
    document["result"] = result;
    document["op_seconds"] = seconds;
    try {
        write_text(out_dir + "/result.json", coolercpp::json::dumps(document));
    } catch (const std::exception& e) {
        std::fprintf(stderr, "%s\n", e.what());
        return 1;
    }
    return 0;
}
