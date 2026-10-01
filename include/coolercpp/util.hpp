// cooler.util helpers that operate on bin tables.

#ifndef COOLERCPP_UTIL_HPP
#define COOLERCPP_UTIL_HPP

#include <cstdint>
#include <functional>
#include <iosfwd>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "coolercpp/json.hpp"
#include "coolercpp/region.hpp"
#include "coolercpp/table.hpp"

namespace coolercpp {

// cooler.util.get_binsize: the single width shared by every bin except the
// last one of each chromosome, as a number carrying the dtype of end - start;
// null when the widths differ or no chromosome has more than one bin.
[[nodiscard]] json::Value get_binsize(const Table& bins);

// cooler.util.get_chromsizes: the end of the last bin of each chromosome, in
// the order of those last bins.
[[nodiscard]] ChromSizes get_chromsizes(const Table& bins);

// cooler.util.natsorted.
[[nodiscard]] std::vector<std::string> natsorted(std::vector<std::string> items);

// cooler.util.argnatsort: the positions that put the names in natural order.
// Only the first min(len(key)) pieces of the sort keys are compared, because
// cooler zips the keys into columns and zip stops at the shortest one; names
// that agree on that prefix keep their original order.
[[nodiscard]] std::vector<std::int64_t> argnatsort(const std::vector<std::string>& names);

// cooler.util.partition: [start, stop) in steps of step, as end point pairs.
[[nodiscard]] std::vector<std::pair<std::int64_t, std::int64_t>> partition(std::int64_t start,
                                                                          std::int64_t stop,
                                                                          std::int64_t step);

// The keyword arguments of cooler.util.read_chromsizes.
struct ReadChromsizesOptions {
    // Each pattern's matching names are kept as one group, sorted naturally;
    // a name matching two patterns is kept twice, as in cooler.
    std::vector<std::string> name_patterns = {"^chr[0-9]+$", "^chr[XY]$", "^chrM$"};
    bool all_names = false;
};

// cooler.util.read_chromsizes: the first two tab separated columns of a
// <db>.chrom.sizes or <db>.chromInfo.txt file, as a name to length mapping. A
// path ending in .gz is decompressed. cooler also accepts a URL and a file
// object; a C++ caller opens those itself and uses the stream overload.
[[nodiscard]] ChromSizes read_chromsizes(const std::string& filepath,
                                         const ReadChromsizesOptions& options = {});
[[nodiscard]] ChromSizes read_chromsizes(std::istream& input,
                                         const ReadChromsizesOptions& options = {});

// cooler.util.binnify: a genome divided into bins of binsize, the last bin of
// each sequence truncated to its length. chrom is a categorical column whose
// categories are every name of chromsizes, in its order.
[[nodiscard]] Table binnify(const ChromSizes& chromsizes, std::int64_t binsize);

// cooler.util.check_bins: the bin table with chrom as a categorical column
// over the names of chromsizes. A categorical column whose categories are not
// exactly those names raises (cooler asserts).
[[nodiscard]] Table check_bins(const Table& bins, const ChromSizes& chromsizes);

// cooler.util.bedslice: the rows of a BED-like table with non-overlapping
// intervals that fall inside a region. cooler takes a pandas groupby object
// and looks the group up by name; here the table is grouped by its chrom
// column on the spot, and a chrom with no rows raises KeyError as
// get_group does.
[[nodiscard]] Table bedslice(const Table& bed, const ChromSizes& chromsizes,
                             const Region& region);

// cooler.util.rlencode: run length encoding of a column.
struct RunLengths {
    std::vector<std::int64_t> starts;
    std::vector<std::int64_t> lengths;
    Column values;
};
[[nodiscard]] RunLengths rlencode(const Column& array,
                                  std::optional<std::int64_t> chunksize = std::nullopt);

// cooler.util.mad for a one dimensional array: the median absolute deviation
// from the median. cooler's axis argument only applies to a two dimensional
// numpy array, which has no Column counterpart.
[[nodiscard]] double mad(const Column& data);

// cooler.util.buffered: an iterable of chunks regrouped into chunks of at
// least size rows.
[[nodiscard]] std::function<std::optional<Table>()> buffered(
    std::function<std::optional<Table>()> chunks, std::int64_t size = 10000000);

// cooler.util.GenomeSegmentation: a bin table together with the positions its
// sequences occupy in the genome, grouped by sequence. The grouping observes
// only the sequences the bin table carries, in the order they appear, which is
// what chrom_binoffset counts over.
class GenomeSegmentation {
  public:
    GenomeSegmentation(ChromSizes chromsizes, const Table& bins);

    [[nodiscard]] const ChromSizes& chromsizes() const noexcept { return chromsizes_; }
    // get_binsize(bins): a number with the dtype of end - start, or null.
    [[nodiscard]] const json::Value& binsize() const noexcept { return binsize_; }
    // list(chromsizes.keys())
    [[nodiscard]] const std::vector<std::string>& contigs() const noexcept {
        return chromsizes_.names();
    }
    [[nodiscard]] const Table& bins() const noexcept { return bins_; }
    // idmap[contig]: the position of a contig in chromsizes; KeyError when
    // absent.
    [[nodiscard]] std::int64_t idmap(const std::string& contig) const;
    // The sequences the bin table observes, in the order they appear.
    [[nodiscard]] const std::vector<std::string>& observed() const noexcept { return observed_; }
    // np.r_[0, cumsum(bins per observed sequence)]
    [[nodiscard]] const std::vector<std::int64_t>& chrom_binoffset() const noexcept {
        return chrom_binoffset_;
    }
    // np.r_[0, cumsum(chromsizes.values)]
    [[nodiscard]] const std::vector<std::int64_t>& chrom_abspos() const noexcept {
        return chrom_abspos_;
    }
    // The absolute genomic start position of every bin.
    [[nodiscard]] const std::vector<std::int64_t>& start_abspos() const noexcept {
        return start_abspos_;
    }
    // The bins of one region, trimmed as bedslice trims them.
    [[nodiscard]] Table fetch(const Region& region) const;
    // The rows of one observed sequence (groupby.get_group).
    [[nodiscard]] Table group(const std::string& contig) const;

  private:
    ChromSizes chromsizes_;
    Table bins_;
    json::Value binsize_;
    std::vector<std::string> observed_;
    // [lo, hi) rows of every observed sequence, in the order observed.
    std::vector<std::pair<std::size_t, std::size_t>> extents_;
    std::vector<std::int64_t> chrom_binoffset_;
    std::vector<std::int64_t> chrom_abspos_;
    std::vector<std::int64_t> start_abspos_;
};

// cooler.util.balanced_partition: regions of roughly equal load, at most
// n_chunk_max of them for the most loaded sequence. loadings is cooler's
// Series of per sequence weights; empty is Python's None, which weighs every
// sequence by its number of bins.
[[nodiscard]] std::vector<RegionTuple> balanced_partition(
    const GenomeSegmentation& gs, std::int64_t n_chunk_max,
    const std::vector<std::string>& file_contigs,
    const std::vector<std::pair<std::string, double>>& loadings = {});

}  // namespace coolercpp

#endif  // COOLERCPP_UTIL_HPP
