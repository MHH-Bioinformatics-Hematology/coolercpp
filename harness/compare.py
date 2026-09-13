"""Comparators of the coolercpp equivalence harness.

compare_results(py, cpp, py_dir, cpp_dir, options) compares two result
documents produced by oracle.py and the C++ driver. describe_file(path)
renders an HDF5 file as a canonical structure (object tree, datatypes,
shapes, chunking, filter pipelines, fill values, decoded data digests and
attributes) and compare_descriptions diffs two of them.
"""

import hashlib
import math
import os

import h5py
import numpy as np

# Attributes whose value changes from run to run and is therefore not
# compared. Their presence and datatype still are.
NORMALIZED_ATTRIBUTES = ["creation-date"]

# Float agreement classes, strictest first. The owner's gate is three
# significant digits per item.
FLOAT_CLASSES = [
    ("bit-exact", 0.0),
    ("rel<=1e-15", 1e-15),
    ("rel<=1e-12", 1e-12),
    ("rel<=1e-9", 1e-9),
    ("rel<=1e-6", 1e-6),
    ("3-significant-digits", 1e-3),
]
CLASS_RANK = {name: i for i, (name, _) in enumerate(FLOAT_CLASSES)}

PYTHON_TYPE_ALIASES = {"FileNotFoundError": "OSError", "JSONDecodeError": "ValueError"}


class Comparison:
    def __init__(self):
        self.problems = []
        self.float_class = None
        self.message_class = None

    def fail(self, where, text):
        self.problems.append(f"{where}: {text}")

    def floats(self, klass):
        if self.float_class is None or CLASS_RANK[klass] > CLASS_RANK[self.float_class]:
            self.float_class = klass

    @property
    def ok(self):
        return not self.problems


def float_agreement(a, b):
    """The strictest class two float arrays of equal shape reach, or None."""
    a = np.asarray(a, dtype=np.float64).ravel()
    b = np.asarray(b, dtype=np.float64).ravel()
    nan_a, nan_b = np.isnan(a), np.isnan(b)
    if not np.array_equal(nan_a, nan_b):
        return None
    keep = ~nan_a
    a, b = a[keep], b[keep]
    if a.tobytes() == b.tobytes():
        return "bit-exact"
    inf = np.isinf(a) | np.isinf(b)
    if not np.array_equal(a[inf], b[inf]):
        return None
    a, b = a[~inf], b[~inf]
    scale = np.maximum(np.abs(a), np.abs(b))
    diff = np.abs(a - b)
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.where(diff == 0, 0.0, diff / np.where(scale == 0, 1.0, scale))
    worst = float(rel.max()) if rel.size else 0.0
    for name, bound in FLOAT_CLASSES[1:]:
        if worst <= bound:
            return name
    return None


def compare_arrays(cmp, where, a, b):
    if a.shape != b.shape:
        cmp.fail(where, f"shape {a.shape} != {b.shape}")
        return
    if a.dtype.kind in "fc" or b.dtype.kind in "fc":
        klass = float_agreement(a, b)
        if klass is None:
            cmp.fail(where, "float values differ beyond three significant digits")
        else:
            cmp.floats(klass)
        return
    if not np.array_equal(a, b):
        mismatch = int(np.count_nonzero(a != b)) if a.shape == b.shape else -1
        cmp.fail(where, f"integer values differ ({mismatch} items)")


def load(directory, name):
    return np.load(os.path.join(directory, name))


def compare_tables(cmp, where, py, cpp, py_dir, cpp_dir):
    names_py = [c["name"] for c in py["columns"]]
    names_cpp = [c["name"] for c in cpp["columns"]]
    if names_py != names_cpp:
        cmp.fail(where, f"columns {names_py} != {names_cpp}")
        return
    for cp, cc in zip(py["columns"], cpp["columns"]):
        at = f"{where}.{cp['name']}"
        if cp["dtype"] != cc["dtype"]:
            cmp.fail(at, f"dtype {cp['dtype']} != {cc['dtype']}")
            continue
        if cp["dtype"] == "object":
            if cp["values"] != cc["values"]:
                cmp.fail(at, "string values differ")
        elif cp["dtype"] == "category":
            if cp["categories"] != cc["categories"]:
                cmp.fail(at, "categories differ")
            if cp.get("ordered") != cc.get("ordered"):
                cmp.fail(at, "ordered flag differs")
            compare_arrays(cmp, at + ".codes", load(py_dir, cp["codes"]), load(cpp_dir, cc["codes"]))
        else:
            a, b = load(py_dir, cp["file"]), load(cpp_dir, cc["file"])
            if a.dtype != b.dtype:
                cmp.fail(at, f"array dtype {a.dtype} != {b.dtype}")
            compare_arrays(cmp, at, a, b)
    index_py = expand_index(py["index"], py_dir)
    index_cpp = expand_index(cpp["index"], cpp_dir)
    if not np.array_equal(index_py, index_cpp):
        cmp.fail(where + ".index", "index differs")


def expand_index(entry, directory):
    if entry["kind"] == "range":
        return np.arange(entry["start"], entry["start"] + entry["size"], dtype=np.int64)
    return load(directory, entry["file"]).astype(np.int64)


def compare_values(cmp, where, a, b):
    if isinstance(a, float) or isinstance(b, float):
        if not (isinstance(a, (int, float)) and isinstance(b, (int, float))):
            cmp.fail(where, f"{a!r} != {b!r}")
        elif not (a == b or (math.isnan(a) and math.isnan(b))):
            klass = float_agreement(np.array([a]), np.array([b]))
            if klass is None:
                cmp.fail(where, f"{a!r} != {b!r}")
            else:
                cmp.floats(klass)
        else:
            cmp.floats("bit-exact")
        return
    if type(a) is not type(b) and not (isinstance(a, int) and isinstance(b, int)):
        cmp.fail(where, f"{a!r} != {b!r}")
        return
    if isinstance(a, dict):
        if list(a.keys()) != list(b.keys()):
            cmp.fail(where, f"keys {list(a)} != {list(b)}")
            return
        for key in a:
            compare_values(cmp, f"{where}.{key}", a[key], b[key])
        return
    if isinstance(a, list):
        if len(a) != len(b):
            cmp.fail(where, f"length {len(a)} != {len(b)}")
            return
        for i, (x, y) in enumerate(zip(a, b)):
            compare_values(cmp, f"{where}[{i}]", x, y)
        return
    if a != b:
        cmp.fail(where, f"{a!r} != {b!r}")


def compare_results(py, cpp, py_dir, cpp_dir, messages="exact", where="result", cmp=None):
    cmp = cmp or Comparison()
    if py.get("type") == "HarnessError" or cpp.get("type") == "HarnessError":
        cmp.fail(where, f"harness error: py={py.get('message')} cpp={cpp.get('message')}")
        return cmp
    if py["kind"] != cpp["kind"]:
        detail = ""
        for side, doc in (("py", py), ("cpp", cpp)):
            if doc["kind"] == "error":
                detail += f" {side} raised {doc['type']}: {doc['message']}"
        cmp.fail(where, f"kind {py['kind']} != {cpp['kind']}.{detail}")
        return cmp
    kind = py["kind"]
    if kind == "error":
        type_py = PYTHON_TYPE_ALIASES.get(py["type"], py["type"])
        if type_py != cpp["type"]:
            cmp.fail(where, f"error type {py['type']} != {cpp['type']}")
        if py["message"] == cpp["message"]:
            if cmp.message_class is None:
                cmp.message_class = "exact message"
        elif messages == "type-only":
            cmp.message_class = "type only"
        else:
            cmp.fail(where, f"error message {py['message']!r} != {cpp['message']!r}")
    elif kind == "value":
        compare_values(cmp, where, py["value"], cpp["value"])
    elif kind == "list":
        if len(py["items"]) != len(cpp["items"]):
            cmp.fail(where, "item count differs")
        for i, (a, b) in enumerate(zip(py["items"], cpp["items"])):
            compare_results(a, b, py_dir, cpp_dir, messages, f"{where}[{i}]", cmp)
    elif kind == "multi":
        if list(py["items"]) != list(cpp["items"]):
            cmp.fail(where, "item names differ")
        for key in py["items"]:
            if key in cpp["items"]:
                compare_results(py["items"][key], cpp["items"][key], py_dir, cpp_dir, messages,
                                f"{where}.{key}", cmp)
    elif kind == "table":
        compare_tables(cmp, where, py, cpp, py_dir, cpp_dir)
    elif kind == "dense":
        if py["shape"] != cpp["shape"] or py["dtype"] != cpp["dtype"]:
            cmp.fail(where, f"dense {py['shape']} {py['dtype']} != {cpp['shape']} {cpp['dtype']}")
        else:
            compare_arrays(cmp, where, load(py_dir, py["file"]), load(cpp_dir, cpp["file"]))
    elif kind == "sparse":
        if py["shape"] != cpp["shape"]:
            cmp.fail(where, f"shape {py['shape']} != {cpp['shape']}")
        for part in ("row", "col", "data"):
            if py[part]["dtype"] != cpp[part]["dtype"]:
                cmp.fail(f"{where}.{part}", f"dtype {py[part]['dtype']} != {cpp[part]['dtype']}")
            else:
                compare_arrays(cmp, f"{where}.{part}", load(py_dir, py[part]["file"]),
                               load(cpp_dir, cpp[part]["file"]))
    elif kind == "created":
        pass
    else:
        cmp.fail(where, f"unknown result kind {kind}")
    return cmp


# ---------------------------------------------------------------------------
# HDF5 structure


def type_description(tid):
    cls = tid.get_class()
    desc = {"class": int(cls), "size": int(tid.get_size())}
    if cls == h5py.h5t.INTEGER:
        desc["order"] = int(tid.get_order())
        desc["sign"] = int(tid.get_sign())
    elif cls == h5py.h5t.FLOAT:
        desc["order"] = int(tid.get_order())
    elif cls == h5py.h5t.STRING:
        desc["variable"] = bool(tid.is_variable_str())
        desc["cset"] = int(tid.get_cset())
        desc["strpad"] = int(tid.get_strpad())
    elif cls == h5py.h5t.ENUM:
        desc["base"] = type_description(tid.get_super())
        desc["members"] = [
            [tid.get_member_name(i).decode("latin-1"), int(tid.get_member_value(i))]
            for i in range(tid.get_nmembers())
        ]
    return desc


def attribute_value(value):
    if isinstance(value, np.ndarray):
        return [attribute_value(v) for v in value.tolist()]
    if isinstance(value, (np.generic,)):
        value = value.item()
    if isinstance(value, bytes):
        return {"bytes": value.decode("latin-1")}
    if isinstance(value, float) and math.isnan(value):
        return "NaN"
    return value


def attributes_description(obj):
    out = {}
    for name in sorted(obj.attrs.keys()):
        aid = h5py.h5a.open(obj.id, name.encode())
        entry = {
            "type": type_description(aid.get_type()),
            "shape": list(aid.get_space().shape) if aid.get_space().get_simple_extent_type() == h5py.h5s.SIMPLE else None,
        }
        if name in NORMALIZED_ATTRIBUTES:
            entry["value"] = "<normalized>"
        else:
            entry["value"] = attribute_value(obj.attrs[name])
        out[name] = entry
    return out


def data_digest(dataset):
    digest = hashlib.sha256()
    n = dataset.shape[0] if dataset.shape else 0
    step = 1 << 22
    if dataset.shape == ():
        digest.update(np.asarray(dataset[()]).tobytes())
        return digest.hexdigest()
    for lo in range(0, n, step):
        block = dataset.id
        array = np.empty((min(n, lo + step) - lo,), dtype=dataset.id.dtype)
        selection = h5py.h5s.create_simple(dataset.shape)
        selection.select_hyperslab((lo,), (array.shape[0],))
        memory = h5py.h5s.create_simple(array.shape)
        block.read(memory, selection, array)
        digest.update(array.tobytes())
    return digest.hexdigest()


def dataset_description(dataset):
    dcpl = dataset.id.get_create_plist()
    filters = []
    for i in range(dcpl.get_nfilters()):
        code, flags, values, name = dcpl.get_filter(i)
        filters.append([int(code), int(flags), [int(v) for v in values], name.decode("latin-1")])
    return {
        "kind": "dataset",
        "type": type_description(dataset.id.get_type()),
        "shape": list(dataset.shape),
        "maxshape": [None if m is None else int(m) for m in dataset.maxshape],
        "chunks": None if dataset.chunks is None else list(dataset.chunks),
        "layout": int(dcpl.get_layout()),
        "filters": filters,
        "fill_value_defined": int(dcpl.fill_value_defined()),
        "fill_time": int(dcpl.get_fill_time()),
        "alloc_time": int(dcpl.get_alloc_time()),
        "fillvalue": attribute_value(dataset.fillvalue),
        "track_times": bool(dcpl.get_obj_track_times()),
        "data_sha256": data_digest(dataset),
        "attrs": attributes_description(dataset),
    }


def describe_file(path):
    tree = {}
    with h5py.File(path, "r") as f:
        tree["/"] = {"kind": "group", "attrs": attributes_description(f)}

        def visit(name, obj):
            key = "/" + name
            if isinstance(obj, h5py.Dataset):
                tree[key] = dataset_description(obj)
            else:
                tree[key] = {"kind": "group", "attrs": attributes_description(obj)}

        f.visititems(visit)
    return tree


# Differences in these dataset properties are recorded but do not fail a case:
# h5py leaves object time tracking of datasets off, coolercpp does too; the
# property is listed so a change becomes visible.
def compare_descriptions(py, cpp):
    problems = []
    for key in sorted(set(py) | set(cpp)):
        if key not in py:
            problems.append(f"{key}: only in the coolercpp file")
            continue
        if key not in cpp:
            problems.append(f"{key}: only in the Python file")
            continue
        a, b = py[key], cpp[key]
        for field in sorted(set(a) | set(b)):
            if a.get(field) != b.get(field):
                if field == "attrs":
                    for name in sorted(set(a["attrs"]) | set(b["attrs"])):
                        if a["attrs"].get(name) != b["attrs"].get(name):
                            problems.append(
                                f"{key} attribute {name}: {a['attrs'].get(name)} != {b['attrs'].get(name)}"
                            )
                else:
                    problems.append(f"{key} {field}: {a.get(field)} != {b.get(field)}")
    return problems
