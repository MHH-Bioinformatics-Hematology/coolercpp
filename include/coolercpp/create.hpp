// cooler.create_cooler (cooler/create/_create.py of cooler 0.10.2).

#ifndef COOLERCPP_CREATE_HPP
#define COOLERCPP_CREATE_HPP

#include <cstdint>
#include <functional>
#include <optional>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "coolercpp/json.hpp"
#include "coolercpp/table.hpp"
#include "coolercpp/version.hpp"

namespace coolercpp {

// One value of the h5opts dictionary (h5py dataset creation keywords).
// std::monostate is Python's None. A tuple such as maxshape=(None,) is a
// vector in which -1 stands for None.
using H5OptValue =
    std::variant<std::monostate, bool, std::int64_t, double, std::string, std::vector<std::int64_t>>;
// The h5opts dict, in insertion order. Recognised keys: chunks, maxshape,
// compression, compression_opts, scaleoffset, shuffle, fletcher32, fillvalue,
// track_times.
using H5Opts = std::vector<std::pair<std::string, H5OptValue>>;

// An iterable of pixel chunks: every call yields the next chunk, nullopt at
// the end. Used for inputs too large for one table.
using PixelChunks = std::function<std::optional<Table>()>;

// Wraps a list of chunks as PixelChunks.
[[nodiscard]] PixelChunks iterate_chunks(std::vector<Table> chunks);

// The keyword arguments of cooler.create_cooler.
struct CreateOptions {
    // Pixel value columns to store. nullopt: only "count".
    std::optional<std::vector<std::string>> columns;
    // Column name -> dtype overrides (a Python dict; empty is None).
    std::vector<std::pair<std::string, DType>> dtypes;
    std::optional<json::Value> metadata;
    std::optional<std::string> assembly;
    bool ordered = false;
    bool symmetric_upper = true;
    std::string mode = "w";
    std::int64_t mergebuf = 20000000;
    bool delete_temp = true;
    std::optional<std::string> temp_dir;
    std::int64_t max_merge = 200;
    bool boundscheck = true;
    bool dupcheck = true;
    bool triucheck = true;
    bool ensure_sorted = false;
    std::optional<H5Opts> h5opts;

    // coolercpp extensions, not present in cooler:
    // A fixed value for the "creation-date" attribute, for reproducible
    // files. nullopt writes datetime.now().isoformat() as cooler does.
    std::optional<std::string> creation_date;
    // The "generated-by" attribute: the program that wrote the file. cooler
    // writes "cooler-<version>"; coolercpp names itself unless the caller
    // passes another string (an application writing through coolercpp may
    // pass its own identity).
    std::string generated_by = std::string("coolercpp-") + kVersion;
};

// create_cooler with a DataFrame or dict of columns: the table is sorted by
// (bin1_id, bin2_id) and written in one step (ordered is implied).
void create_cooler(const std::string& cool_uri, const Table& bins, const Table& pixels,
                   const CreateOptions& options = {});

// create_cooler with an iterable of chunks: written in one pass when
// options.ordered, otherwise through cooler's external sort (temporary
// coolers merged with CoolerMerger).
void create_cooler(const std::string& cool_uri, const Table& bins, PixelChunks pixels,
                   const CreateOptions& options = {});

// The pixel input of one cell of a scool file: a table, or an iterable of
// chunks for an input too large for one table.
using CellPixels = std::variant<Table, PixelChunks>;
// A Python dict of cell name -> table. The order given does not reach the
// file: create_scool sorts the names.
using CellTables = std::vector<std::pair<std::string, Table>>;
using CellPixelTables = std::vector<std::pair<std::string, CellPixels>>;

// cooler.create_scool with one bin table shared by every cell: the root of the
// file holds the chromosome and bin tables and one cooler per cell under
// /cells, each hard linking the shared tables.
void create_scool(const std::string& cool_uri, const Table& bins,
                  const CellPixelTables& cell_name_pixels, const CreateOptions& options = {});

// cooler.create_scool with one bin table per cell (the dict form of the bins
// argument): the shared root bin table is the chrom, start and end columns of
// the first cell's bins, and every cell keeps its own extra bin columns.
void create_scool(const std::string& cool_uri, const CellTables& cell_name_bins,
                  const CellPixelTables& cell_name_pixels, const CreateOptions& options = {});

}  // namespace coolercpp

#endif  // COOLERCPP_CREATE_HPP
