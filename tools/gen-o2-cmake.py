#!/usr/bin/env python3
"""Generate the CMake build of the imported O2 subset.

Input:  extern/o2 (the imported O2 files, with the upstream CMakeLists.txt
        and LinkDef files of each library directory).
Output: cmake/o2/CMakeLists.txt and cmake/o2/linkdef/<Library>LinkDef.h.

The generated build keeps the upstream library split: one shared library
and one ROOT dictionary for each upstream library. It compiles only the
imported sources. A dictionary contains only the classes whose headers
are imported.
"""
import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parent.parent
O2 = ROOT / "extern" / "o2"
OUT = ROOT / "cmake" / "o2"

# Upstream library directory -> list of (library, sources subset or None).
# None means: all imported sources of the directory.
LIBS = {
    "Common/Utils": [("CommonUtils", None)],
    "Common/MathUtils": [("MathUtils", None)],
    "DataFormats/common": [("CommonDataFormat", None)],
    "DataFormats/Detectors/Common": [("DetectorsCommonDataFormats", None)],
    "Detectors/3rdparty/FairRootBase": [("FairRootBase", None)],
    "Common/Field": [("Field", None)],
    "DataFormats/Parameters": [("DataFormatsParameters", None)],
    "Framework/Core": [("Framework", None)],
    "DataFormats/Calibration": [("DataFormatsCalibration", None)],
    "DataFormats/simulation": [("SimulationDataFormat", None)],
    "DataFormats/Detectors/FIT/FV0": [("DataFormatsFV0", None)],
    "Common/SimConfig": [("SimConfig", None)],
    "Detectors/Base": [("DetectorsBase", None)],
    "Detectors/FastSim": [("FastSim", None)],
    "Generators": [("Generators", None)],
    "Steer": [("Steer", None)],
    "Detectors/Passive": [("DetectorsPassive", None)],
    "Detectors/FIT/FV0/base": [("FV0Base", None)],
    "Detectors/FIT/FV0/simulation": [("FV0Simulation", None)],
    "Detectors/gconfig": [
        ("SimSetup", ["src/SimSetup.cxx", "src/SetCuts.cxx"]),
        ("G4Setup", ["src/G4Config.cxx", "src/G4RunConfiguration.cxx",
                     "src/G4LocalFieldConstruction.cxx", "src/VecGeomNavigation.cxx",
                     "src/O2MonopolePhysics.cxx"]),
    ],
}

# Upstream gets these dependencies through libraries that are not part of
# the subset (for example O2::DataFormatsCTP). Add them directly.
EXTRA_DEPENDS = {
    "SimulationDataFormat": ["CommonDataFormat"],
}

PRAGMA = re.compile(r"#pragma\s+link\s+C\+\+\s+(class|struct|function|enum|namespace|typedef|global)\s+(.+?)\s*;")


def upstream_block(cmake_text, lib):
    """Return the text of o2_add_library(<lib> ...) and its dictionary call."""
    flat = re.sub(r"#.*", "", cmake_text)
    add = re.search(r"o2_add_library\(\s*" + lib + r"\b(.*?)\)", flat, re.S)
    dic = re.search(r"o2_target_root_dictionary\(\s*" + lib + r"\b(.*?)\)", flat, re.S)
    return (add.group(1) if add else ""), (dic.group(1) if dic else None)


STD_TEMPLATES = {"vector", "map", "unordered_map", "pair", "array", "set", "unique_ptr"}


def declared(name, texts):
    base = re.sub(r"<.*", "", name).split("::")[-1].strip()
    if not base:
        return False
    if base in STD_TEMPLATES and "<" in name:
        args = re.findall(r"((?:o2::)?[A-Za-z_][\w:]*)", name.split("<", 1)[1])
        user = [a for a in args if a.startswith("o2::") or a.startswith("Fair")]
        return bool(user) and all(declared(a, texts) for a in user)
    pat = re.compile(r"\b(class|struct|enum|namespace|using|typedef)\s+(class\s+)?" + re.escape(base) + r"\b")
    if any(pat.search(t) for t in texts):
        return True
    return re.search(r"\b" + re.escape(base) + r"\b", " ".join(texts)) is not None and "(" in name


def template_args_ok(name, all_texts):
    for arg in re.findall(r"(o2::[\w:]+)", name.split("<", 1)[1] if "<" in name else ""):
        if not declared(arg, all_texts):
            return False
    return True


def reduce_linkdef(src, dict_headers, all_texts):
    texts = [h.read_text(errors="replace") for h in dict_headers]
    out = []
    kept = 0
    for line in src.read_text().splitlines():
        inc = re.match(r'\s*#include\s+[<"]([^>"]+)[>"]', line)
        if inc and not any(h.as_posix().endswith("/" + inc.group(1)) for h in ALL_HEADERS):
            continue
        m = PRAGMA.search(line)
        if not m:
            out.append(line)
            continue
        kind, name = m.groups()
        name = re.sub(r"[+\-!]+$", "", name).strip()
        is_instance = "<" in name and not re.sub(r"<.*", "", name).split("::")[-1].strip() in STD_TEMPLATES
        visible = declared(name, texts) or (is_instance and declared(name, all_texts))
        if kind == "global" or (visible and template_args_ok(name, all_texts)):
            out.append(line)
            kept += 1
    return "\n".join(out) + "\n", kept


ALL_HEADERS = sorted(O2.rglob("*.h"))


def main():
    all_headers = ALL_HEADERS
    all_texts = [h.read_text(errors="replace") for h in all_headers]
    include_dirs = sorted({str(p.relative_to(O2)) for p in O2.rglob("include") if p.is_dir()})
    include_dirs += ["GPU/Common", "GPU/Utils"]

    (OUT / "linkdef").mkdir(parents=True, exist_ok=True)
    lines = []
    w = lines.append
    w("# Generated by tools/gen-o2-cmake.py. Do not edit.")
    w("# Build of the imported O2 subset with the upstream library split.")
    w("")
    w("set(O2_INCLUDE_DIRS")
    for d in include_dirs:
        w(f"  ${{O2_SOURCE_DIR}}/{d}")
    w(")")
    w("")
    built = []
    for d, libs in LIBS.items():
        cm = (O2 / d / "CMakeLists.txt").read_text()
        for lib, subset in libs:
            add, dic = upstream_block(cm, lib)
            srcs = subset or sorted(str(p.relative_to(O2 / d)) for p in (O2 / d / "src").glob("*.cxx"))
            srcs = [s for s in srcs if (O2 / d / s).exists()]
            deps = [x for x in re.findall(r"O2::(\w+)", add) if x in built]
            deps += [x for x in EXTRA_DEPENDS.get(lib, []) if x not in deps]
            w(f"o2_library({lib}")
            w(f"  DIR {d}")
            w("  SOURCES " + " ".join(srcs))
            if deps:
                w("  DEPENDS " + " ".join(deps))
            if dic is not None:
                hdrs = [h for h in re.findall(r"(include/[\w/.]+\.h)", dic) if (O2 / d / h).exists()]
                if not hdrs and "GLOB" in cm:
                    hdrs = sorted(str(h.relative_to(O2 / d)) for h in (O2 / d / "include").glob("*.h"))
                if lib == "Generators":
                    hdrs = sorted(h for h in re.findall(r"(include/Generators/\w+\.h)", cm)
                                  if (O2 / d / h).exists())
                ld = re.search(r"LINKDEF\s+(\S+)", dic)
                ld = O2 / d / (ld.group(1) if ld else f"src/{lib}LinkDef.h")
                if not ld.exists():
                    cands = list((O2 / d / "src").glob("*LinkDef*.h"))
                    ld = cands[0] if len(cands) == 1 else None
                if hdrs and ld:
                    text, kept = reduce_linkdef(ld, [O2 / d / h for h in hdrs], all_texts)
                    out = OUT / "linkdef" / f"{lib}LinkDef.h"
                    out.write_text(f"// Reduced from {d}/{ld.relative_to(O2 / d)} by tools/gen-o2-cmake.py.\n" + text)
                    w("  DICTIONARY " + " ".join(hdrs))
            w(")")
            w("")
            built.append(lib)
    (OUT / "o2libs.cmake").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
