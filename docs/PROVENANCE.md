# Provenance

## Licence

coolercpp is BSD-3-Clause (`LICENSE`), the licence of cooler, whose API, file
format and ported code it follows. The notices the projects below require are
reproduced at the end of this document and summarised in `NOTICE`.

## Code taken from HiCExplorer v4

The HiCExplorer v4 C++ port (`HiCExplorer-v4/cpp/core`) had a validated cool
reader and writer. The parts below were copied and adapted; coolercpp does not
include any HiCExplorer header or link any HiCExplorer target.

| coolercpp file | origin |
|---|---|
| `src/h5.hpp`, `src/h5.cpp` | `hdf5_util.hpp/.cpp`: the `Handle` class, `guess_chunk`, the fixed string and enum type builders and the attribute reader. The blosc filter was left out; h5py file modes and dataset creation were added. |
| `src/json.cpp`, `include/coolercpp/json.hpp` | `json_lite.hpp/.cpp`, reworked for ordered objects, strict simplejson parsing and `simplejson.dumps` output. |
| `src/numpy_compat.hpp`, `src/numpy_compat.cpp` | `numpy_compat.hpp/.cpp`: `pairwise_sum` and `float_repr`. |
| `src/create.cpp` | The conventions validated in `cool_file.cpp` (`write_cool`): h5py's chunk rule, gzip-6 without shuffle for extra bin columns, the ENUM chrom column, numpy's accumulation order for the `sum` attribute. |

## Code derived from BSD-3-Clause projects

Ports of cooler 0.10.2 (https://github.com/open2c/cooler, BSD-3-Clause):

| coolercpp file | cooler source |
|---|---|
| `include/coolercpp/api.hpp`, `src/api.cpp` | `cooler/api.py`, `cooler/core/_selectors.py`, `cooler/core/_tableops.py` (`get`), `cooler/fileops.py` (`list_coolers`) |
| `src/rangequery.hpp`, `src/rangequery.cpp` | `cooler/core/_rangequery.py` |
| `include/coolercpp/create.hpp`, `src/create.cpp` | `cooler/create/_create.py`, `cooler/create/_ingest.py` (`validate_pixels`), `cooler/core/_tableops.py` (`put`), `cooler/reduce.py` (`merge_breakpoints`, `CoolerMerger`), `cooler/util.py` (`rlencode`) |
| `include/coolercpp/region.hpp`, `src/region.cpp` | `cooler/util.py` (`parse_cooler_uri`, `parse_humanized`, `parse_region_string`, `parse_region`) |
| `include/coolercpp/util.hpp`, `src/util.cpp` | `cooler/util.py` (`get_binsize`, `get_chromsizes`, `natsort_key`, `natsorted`, `partition`, `mad`) |
| `include/coolercpp/balance.hpp`, `src/balance.cpp` | `cooler/balance.py`, `cooler/parallel.py` (`chunkgetter`, `split`, `MultiplexDataPipe`) |
| `tests/data/*.cool`, `tests/data/*.mcool`, `tests/data/toy.chrom.sizes` | `tests/data/` of cooler 0.10.2 |

Behaviour reproduced from other BSD-3-Clause projects, implemented from their
documented or observed semantics:

- h5py 3.12 (`h5py/_hl/filters.py` `guess_chunk` and `fill_dcpl`,
  `h5py/_hl/dataset.py` `make_new_dset`): `src/h5.cpp`, `src/create.cpp`.
- numpy 1.26 (pairwise summation, `mean`, `var`, `median`, `linspace`, dtype
  promotion, slice resolution): `src/numpy_compat.cpp`, `src/dtype.cpp`.
- pandas 2.2 (`group_sum` Kahan summation, `Categorical.from_codes` checks,
  label slicing in `annotate`): `src/create.cpp`, `src/api.cpp`.
- scipy 1.14 (`coo_matrix` index dtype selection and bounds checks):
  `src/rangequery.cpp`.

## cooler's licence notice

    BSD 3-Clause License

    Copyright (c) 2015-2023, Cooler developers
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice, this
      list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    * Neither the name of the copyright holder nor the names of its
      contributors may be used to endorse or promote products derived from
      this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
    FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
