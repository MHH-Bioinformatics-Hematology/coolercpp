#include "npy.hpp"

#include <cstdio>
#include <cstring>
#include <fstream>
#include <stdexcept>

namespace harness {

using coolercpp::Column;
using coolercpp::DType;

namespace {

const char* descr_of(DType dtype) {
    switch (dtype) {
        case DType::Bool: return "|b1";
        case DType::Int8: return "|i1";
        case DType::Int16: return "<i2";
        case DType::Int32: return "<i4";
        case DType::Int64: return "<i8";
        case DType::UInt8: return "|u1";
        case DType::UInt16: return "<u2";
        case DType::UInt32: return "<u4";
        case DType::UInt64: return "<u8";
        case DType::Float32: return "<f4";
        case DType::Float64: return "<f8";
        default: throw std::runtime_error("npy: no descr for a non-numeric column");
    }
}

DType dtype_of(const std::string& descr) {
    std::string d = descr;
    if (!d.empty() && (d[0] == '<' || d[0] == '|' || d[0] == '=')) {
        d.erase(0, 1);
    }
    if (d == "b1") return DType::Bool;
    if (d == "i1") return DType::Int8;
    if (d == "i2") return DType::Int16;
    if (d == "i4") return DType::Int32;
    if (d == "i8") return DType::Int64;
    if (d == "u1") return DType::UInt8;
    if (d == "u2") return DType::UInt16;
    if (d == "u4") return DType::UInt32;
    if (d == "u8") return DType::UInt64;
    if (d == "f4") return DType::Float32;
    if (d == "f8") return DType::Float64;
    throw std::runtime_error("npy: unsupported descr " + descr);
}

const void* raw_data(const Column& column) {
    return std::visit(
        [](const auto& v) -> const void* {
            using V = std::decay_t<decltype(v)>;
            if constexpr (std::is_same_v<V, coolercpp::CategoricalData> ||
                          std::is_same_v<V, std::vector<std::string>>) {
                return nullptr;
            } else {
                return v.data();
            }
        },
        column.data());
}

void* raw_data(Column& column) {
    return std::visit(
        [](auto& v) -> void* {
            using V = std::decay_t<decltype(v)>;
            if constexpr (std::is_same_v<V, coolercpp::CategoricalData> ||
                          std::is_same_v<V, std::vector<std::string>>) {
                return nullptr;
            } else {
                return v.data();
            }
        },
        column.data());
}

}  // namespace

void write_npy(const std::string& path, const Column& column,
               const std::vector<std::int64_t>& shape) {
    std::string shape_text = "(";
    for (std::size_t i = 0; i < shape.size(); ++i) {
        shape_text += std::to_string(shape[i]);
        shape_text += (shape.size() == 1 || i + 1 < shape.size()) ? "," : "";
        if (i + 1 < shape.size()) {
            shape_text += " ";
        }
    }
    shape_text += ")";
    std::string header = std::string("{'descr': '") + descr_of(column.dtype()) +
                         "', 'fortran_order': False, 'shape': " + shape_text + ", }";
    const std::size_t total = 10 + header.size() + 1;
    header.append((64 - total % 64) % 64, ' ');
    header.push_back('\n');
    std::ofstream out(path, std::ios::binary);
    if (!out) {
        throw std::runtime_error("npy: cannot write " + path);
    }
    out.write("\x93NUMPY\x01\x00", 8);
    const auto length = static_cast<std::uint16_t>(header.size());
    out.write(reinterpret_cast<const char*>(&length), 2);
    out.write(header.data(), static_cast<std::streamsize>(header.size()));
    out.write(static_cast<const char*>(raw_data(column)),
              static_cast<std::streamsize>(column.size() * coolercpp::itemsize(column.dtype())));
    if (!out) {
        throw std::runtime_error("npy: short write to " + path);
    }
}

void write_npy(const std::string& path, const Column& column) {
    write_npy(path, column, {static_cast<std::int64_t>(column.size())});
}

Column read_npy(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("npy: cannot open " + path);
    }
    char magic[8];
    in.read(magic, 8);
    if (std::memcmp(magic, "\x93NUMPY", 6) != 0) {
        throw std::runtime_error("npy: not an npy file: " + path);
    }
    std::size_t header_length = 0;
    if (magic[6] == 1) {
        std::uint16_t n = 0;
        in.read(reinterpret_cast<char*>(&n), 2);
        header_length = n;
    } else {
        std::uint32_t n = 0;
        in.read(reinterpret_cast<char*>(&n), 4);
        header_length = n;
    }
    std::string header(header_length, '\0');
    in.read(header.data(), static_cast<std::streamsize>(header_length));
    const auto field = [&header](const std::string& key) {
        const std::size_t k = header.find("'" + key + "'");
        if (k == std::string::npos) {
            throw std::runtime_error("npy: header lacks " + key);
        }
        return header.substr(header.find(':', k) + 1);
    };
    std::string descr = field("descr");
    descr = descr.substr(descr.find('\'') + 1);
    descr = descr.substr(0, descr.find('\''));
    if (field("fortran_order").find("True") < field("fortran_order").find(',')) {
        throw std::runtime_error("npy: Fortran order is not supported");
    }
    std::string shape = field("shape");
    shape = shape.substr(shape.find('(') + 1);
    shape = shape.substr(0, shape.find(')'));
    std::int64_t count = 1;
    std::size_t pos = 0;
    bool any = false;
    while (pos < shape.size()) {
        const std::size_t comma = shape.find(',', pos);
        const std::string part = shape.substr(pos, comma == std::string::npos ? std::string::npos : comma - pos);
        if (part.find_first_of("0123456789") != std::string::npos) {
            count *= std::stoll(part);
            any = true;
        }
        if (comma == std::string::npos) {
            break;
        }
        pos = comma + 1;
    }
    if (!any) {
        count = 1;
    }
    Column column = Column::empty(dtype_of(descr), static_cast<std::size_t>(count));
    in.read(static_cast<char*>(raw_data(column)),
            static_cast<std::streamsize>(column.size() * coolercpp::itemsize(column.dtype())));
    if (!in) {
        throw std::runtime_error("npy: short read from " + path);
    }
    return column;
}

}  // namespace harness
