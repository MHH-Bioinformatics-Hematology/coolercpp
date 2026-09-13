"""Python side of the coolercpp equivalence harness.

    python oracle.py run CASE.json OUT_DIR     run a case with cooler 0.10.2
    python oracle.py prep SPEC.json OUT_DIR    build the input tables of a
                                               create case (not measured)

The outputs use the same neutral format as the C++ driver: result.json plus
.npy files.
"""

import json
import os
import sys
import time

import numpy as np
import pandas as pd

HARNESS_ERROR = "HarnessError"


# ---------------------------------------------------------------------------
# neutral tables


class Output:
    def __init__(self, out_dir):
        self.dir = out_dir
        self.counter = 0

    def next_file(self):
        name = f"a{self.counter}.npy"
        self.counter += 1
        return name

    def save(self, array):
        name = self.next_file()
        np.save(os.path.join(self.dir, name), np.ascontiguousarray(array))
        return name

    def table(self, obj):
        if isinstance(obj, pd.Series):
            frame = obj.to_frame()
        else:
            frame = obj
        columns = []
        for position in range(frame.shape[1]):
            name = frame.columns[position]
            series = frame.iloc[:, position]
            dtype = series.dtype
            entry = {"name": str(name)}
            if isinstance(dtype, pd.CategoricalDtype):
                entry["dtype"] = "category"
                entry["codes"] = self.save(series.cat.codes.to_numpy().astype(np.int32))
                entry["categories"] = [str(x) for x in dtype.categories]
                entry["ordered"] = bool(dtype.ordered)
            elif dtype == object:
                entry["dtype"] = "object"
                entry["values"] = [str(x) for x in series.tolist()]
            else:
                entry["dtype"] = dtype.name
                entry["file"] = self.save(series.to_numpy())
            columns.append(entry)
        index = frame.index
        if isinstance(index, pd.RangeIndex) and index.step == 1:
            index_entry = {"kind": "range", "start": int(index.start), "size": len(index)}
        else:
            index_entry = {"kind": "values", "file": self.save(np.asarray(index, dtype=np.int64))}
        return {"kind": "table", "columns": columns, "index": index_entry}

    def matrix(self, result):
        from scipy.sparse import coo_matrix

        if isinstance(result, (pd.DataFrame, pd.Series)):
            return self.table(result)
        if isinstance(result, coo_matrix):
            out = {"kind": "sparse", "shape": [int(result.shape[0]), int(result.shape[1])]}
            for name, array in (("row", result.row), ("col", result.col), ("data", result.data)):
                out[name] = {"file": self.save(array), "dtype": array.dtype.name}
            return out
        array = np.asarray(result)
        return {
            "kind": "dense",
            "shape": [int(array.shape[0]), int(array.shape[1])],
            "dtype": array.dtype.name,
            "file": self.save(array),
        }


def read_table(directory):
    with open(os.path.join(directory, "meta.json")) as handle:
        meta = json.load(handle)
    data = {}
    for entry in meta["columns"]:
        name = entry["name"]
        if entry["dtype"] == "object":
            data[name] = pd.Series(entry["values"], dtype=object)
        elif entry["dtype"] == "category":
            codes = np.load(os.path.join(directory, entry["codes"]))
            data[name] = pd.Categorical.from_codes(
                codes, entry["categories"], ordered=entry.get("ordered", True)
            )
        else:
            data[name] = np.load(os.path.join(directory, entry["file"]))
    return pd.DataFrame(data, columns=[e["name"] for e in meta["columns"]])


def write_input_table(frame, directory):
    os.makedirs(directory, exist_ok=True)
    out = Output(directory)
    meta = out.table(frame)
    with open(os.path.join(directory, "meta.json"), "w") as handle:
        json.dump(meta, handle)


# ---------------------------------------------------------------------------
# results


def error_value(exc):
    import cooler.create

    kind = type(exc).__name__
    if isinstance(exc, cooler.create.BadInputError):
        kind = "BadInputError"
    elif isinstance(exc, OSError):
        kind = "OSError"
    message = exc.args[0] if exc.args and isinstance(exc.args[0], str) else str(exc)
    return {"kind": "error", "type": kind, "message": message}


def capture(body):
    try:
        return body()
    except Exception as exc:  # noqa: BLE001 - every library error is a result
        return error_value(exc)


def value_result(value):
    return {"kind": "value", "value": value}


def tagged(value):
    if value is None:
        return {"t": "null"}
    if isinstance(value, (bool, np.bool_)):
        return {"t": "bool", "v": bool(value)}
    if isinstance(value, np.integer):
        return {"t": "int", "dtype": value.dtype.name, "v": int(value)}
    if isinstance(value, int):
        return {"t": "int", "dtype": "int64", "v": value}
    if isinstance(value, np.floating):
        return {"t": "float", "dtype": value.dtype.name, "v": float(value)}
    if isinstance(value, float):
        return {"t": "float", "dtype": "float64", "v": value}
    if isinstance(value, str):
        return {"t": "str", "v": value}
    if isinstance(value, (bytes, np.bytes_)):
        return {"t": "bytes", "v": bytes(value).decode("latin-1")}
    if isinstance(value, dict):
        return {"t": "dict", "v": [[str(k), tagged(v)] for k, v in value.items()]}
    if isinstance(value, (list, tuple, np.ndarray)):
        return {"t": "list", "v": [tagged(v) for v in list(value)]}
    raise TypeError(f"cannot tag {type(value)}")


def region_from(value):
    if isinstance(value, str):
        return value
    return tuple(value)


def slice_from(items):
    return slice(*items)


def axis_from(value):
    if isinstance(value, list):
        return slice_from(value)
    return value


# ---------------------------------------------------------------------------
# operations


def op_info(spec, output):
    import cooler

    c = cooler.Cooler(spec["uri"])
    items = {
        "info": capture(lambda: value_result(tagged(c.info))),
        "chromnames": capture(lambda: value_result([str(x) for x in c.chromnames])),
        "chromsizes": capture(
            lambda: value_result(
                {
                    "names": [str(x) for x in c.chromsizes.index],
                    "lengths": [int(x) for x in c.chromsizes.values],
                    "dtype": c.chromsizes.dtype.name,
                }
            )
        ),
        "binsize": capture(
            lambda: value_result(None if c.binsize is None else int(c.binsize))
        ),
        "storage_mode": capture(lambda: value_result(c.storage_mode)),
        "shape": capture(lambda: value_result([int(x) for x in c.shape])),
        "root": value_result(c.root),
    }
    return {"kind": "multi", "items": items}


def op_extent(spec, output):
    import cooler

    c = cooler.Cooler(spec["uri"])
    items = []
    for region in spec["regions"]:
        if spec["op"] == "offset":
            items.append(capture(lambda: value_result(int(c.offset(region_from(region))))))
        else:
            items.append(
                capture(lambda: value_result([int(x) for x in c.extent(region_from(region))]))
            )
    return {"kind": "list", "items": items}


def op_table(spec, output):
    import cooler

    c = cooler.Cooler(spec["uri"])
    kwargs = spec.get("selector_kwargs") or {}
    which = spec["selector"]
    if which == "chroms":
        sel = c.chroms(**({"convert_enum": kwargs["convert_enum"]} if "convert_enum" in kwargs else {}))
    elif which == "bins":
        sel = c.bins(**({"convert_enum": kwargs["convert_enum"]} if "convert_enum" in kwargs else {}))
    else:
        extra = {"convert_enum": kwargs["convert_enum"]} if "convert_enum" in kwargs else {}
        sel = c.pixels(join=kwargs.get("join", False), **extra)
    fields = spec.get("fields")
    if fields is not None:
        sel = sel[fields]
    access = spec["access"]
    if "slice" in access:
        return output.table(sel[slice_from(access["slice"])])
    if "row" in access:
        return output.table(sel[access["row"]])
    if "fetch" in access:
        return output.table(sel.fetch(region_from(access["fetch"])))
    if "columns" in access:
        return value_result([str(x) for x in sel.columns])
    if "dtypes" in access:
        dtypes = sel.dtypes
        if isinstance(dtypes, pd.Series):
            pairs = [[str(k), "category" if isinstance(v, pd.CategoricalDtype) else str(v)] for k, v in dtypes.items()]
        else:
            name = fields if isinstance(fields, str) else ""
            pairs = [[name, "category" if isinstance(dtypes, pd.CategoricalDtype) else str(dtypes)]]
        return value_result(pairs)
    if "len" in access:
        return value_result(len(sel))
    raise ValueError("unknown table access")


def op_matrix(spec, output):
    import cooler

    c = cooler.Cooler(spec["uri"])
    options = dict(spec.get("options") or {})
    sel = c.matrix(**options)
    access = spec["access"]
    if "slice" in access:
        axes = access["slice"]
        if len(axes) == 1:
            return output.matrix(sel[axis_from(axes[0])])
        return output.matrix(sel[axis_from(axes[0]), axis_from(axes[1])])
    if "fetch" in access:
        regions = access["fetch"]
        if len(regions) > 1 and regions[1] is not None:
            return output.matrix(sel.fetch(region_from(regions[0]), region_from(regions[1])))
        return output.matrix(sel.fetch(region_from(regions[0])))
    raise ValueError("unknown matrix access")


def op_parse_region(spec, output):
    import cooler
    from cooler.util import parse_region

    sizes = None
    if spec.get("chromsizes_from"):
        sizes = cooler.Cooler(spec["chromsizes_from"]).chromsizes
    items = []
    for region in spec["regions"]:
        def body(region=region):
            chrom, start, end = parse_region(region_from(region), sizes)
            return value_result([str(chrom), int(start), int(end)])

        items.append(capture(body))
    return {"kind": "list", "items": items}


def op_dump(spec, output):
    import cooler

    c = cooler.Cooler(spec["uri"])
    items = {
        "info": capture(lambda: value_result(tagged(c.info))),
        "chroms": capture(lambda: output.table(c.chroms()[:])),
        "bins": capture(lambda: output.table(c.bins()[:])),
        "pixels": capture(lambda: output.table(c.pixels()[:])),
    }
    return {"kind": "multi", "items": items}


def op_create(spec, output):
    import cooler

    for step in spec["steps"]:
        bins = read_table(step["bins"])
        options = dict(step.get("options") or {})
        if "h5opts" in options and options["h5opts"] is not None:
            options["h5opts"] = {
                k: (tuple(v) if isinstance(v, list) else v) for k, v in options["h5opts"].items()
            }
        if "pixels" in step:
            pixels = read_table(step["pixels"])
        else:
            pixels = (read_table(d) for d in step["pixel_chunks"])
        cooler.create_cooler(step["uri"], bins, pixels, **options)
    return {"kind": "created"}


OPS = {
    "info": op_info,
    "extent": op_extent,
    "offset": op_extent,
    "table": op_table,
    "matrix": op_matrix,
    "parse_region": op_parse_region,
    "dump": op_dump,
    "create": op_create,
}


def run(case_path, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    start = time.monotonic()
    try:
        with open(case_path) as handle:
            spec = json.load(handle)
        output = Output(out_dir)
        if spec["op"] == "list_coolers":
            import cooler.fileops

            result = capture(lambda: value_result(cooler.fileops.list_coolers(spec["path"])))
        else:
            result = capture(lambda: OPS[spec["op"]](spec, output))
    except Exception as exc:  # noqa: BLE001
        result = {"kind": "error", "type": HARNESS_ERROR, "message": repr(exc)}
    document = {"result": result, "op_seconds": time.monotonic() - start}
    with open(os.path.join(out_dir, "result.json"), "w") as handle:
        json.dump(document, handle)


# ---------------------------------------------------------------------------
# input preparation for create cases


def prep(spec_path, out_dir):
    import cooler

    with open(spec_path) as handle:
        spec = json.load(handle)
    src = cooler.Cooler(spec["source"])
    bins = src.bins()[:]
    if spec.get("bins_columns") is not None:
        bins = bins[spec["bins_columns"]]
    if spec.get("chrom") == "string":
        bins["chrom"] = bins["chrom"].astype(str).astype(object)
    pixels = src.pixels()[:]
    if spec.get("pixel_columns") is not None:
        pixels = pixels[spec["pixel_columns"]]
    n_bins = len(bins)
    for transform in spec.get("transforms", []):
        op = transform["op"]
        if op == "scale_count":
            pixels["count"] = (pixels["count"] * transform["factor"]).astype(transform.get("dtype", "float64"))
        elif op == "add_column":
            source = pixels[transform.get("from", "count")]
            pixels[transform["name"]] = (source * transform.get("factor", 1.0)).astype(transform["dtype"])
        elif op == "astype":
            pixels[transform["column"]] = pixels[transform["column"]].astype(transform["dtype"])
        elif op == "shuffle":
            order = np.random.default_rng(transform.get("seed", 0)).permutation(len(pixels))
            pixels = pixels.iloc[order].reset_index(drop=True)
        elif op == "head":
            pixels = pixels.iloc[: transform["n"]].reset_index(drop=True)
        elif op == "duplicate_row":
            row = pixels.iloc[[transform["row"]]]
            pixels = pd.concat([pixels, row], ignore_index=True)
        elif op == "swap_ids_row":
            r = transform["row"]
            b1, b2 = pixels.at[r, "bin1_id"], pixels.at[r, "bin2_id"]
            pixels.at[r, "bin1_id"], pixels.at[r, "bin2_id"] = b2, b1
        elif op == "set_value":
            pixels.at[transform["row"], transform["column"]] = (
                n_bins if transform["value"] == "n_bins" else transform["value"]
            )
        elif op == "full_square":
            off = pixels[pixels["bin1_id"] != pixels["bin2_id"]].copy()
            off["bin1_id"], off["bin2_id"] = off["bin2_id"].to_numpy(), off["bin1_id"].to_numpy()
            pixels = (
                pd.concat([pixels, off], ignore_index=True)
                .sort_values(["bin1_id", "bin2_id"], kind="stable")
                .reset_index(drop=True)
            )
        else:
            raise ValueError(f"unknown transform {op}")
    write_input_table(bins.reset_index(drop=True), os.path.join(out_dir, "bins"))
    chunks = spec.get("chunks")
    if chunks:
        parts = np.array_split(np.arange(len(pixels)), chunks)
        if spec.get("reverse_chunks"):
            parts = parts[::-1]
        for i, positions in enumerate(parts):
            part = pixels.iloc[positions].reset_index(drop=True)
            write_input_table(part, os.path.join(out_dir, f"pixels_{i}"))
    else:
        write_input_table(pixels.reset_index(drop=True), os.path.join(out_dir, "pixels"))


if __name__ == "__main__":
    if len(sys.argv) != 4 or sys.argv[1] not in ("run", "prep"):
        sys.stderr.write(__doc__)
        sys.exit(2)
    if sys.argv[1] == "run":
        run(sys.argv[2], sys.argv[3])
    else:
        prep(sys.argv[2], sys.argv[3])
