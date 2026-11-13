#!/usr/bin/env python3
"""Split a callgrind profile into the layers of the simulation stack.

For each layer, the script reports the self cost (instructions executed
in the code of that layer) and the number of calls that cross from one
layer into another. The layer of a function is the layer of the shared
object that holds it.

Usage: layers.py <callgrind.out> [--minus <callgrind.out>] [--per N] [--json out.json]

With --minus and --per, the numbers are the cost per event of the event loop:
profile a run with N1 events and one with N2 events, and use --per N1-N2.
"""
import argparse
import collections
import json
import os
import re
import sys

LAYERS = [
    ("VecGeom", re.compile(r"lib(vecgeom|tgeo2vecgeom)\w*\.so")),
    ("VGM", re.compile(r"lib(BaseVGM|ClhepVGM|Geant4GM|RootGM|XmlVGM)\w*\.so")),
    ("Geant4", re.compile(r"libG4\w+\.so")),
    ("geant4_vmc", re.compile(r"lib(geant4vmc|g4root|mtroot|g4vmc)\w*\.so")),
    ("VMC", re.compile(r"libVMC\w*\.so")),
    ("FairRoot", re.compile(r"libO2FairRootBase\.so")),
    ("O2", re.compile(r"(libO2\w+\.so|/fv0sim$)")),
    ("Pythia8", re.compile(r"libpythia8\.so")),
    ("ROOT", re.compile(r"lib(Core|Cling|RIO|Tree|Geom|Hist|MathCore|Matrix|Physics|EG|GenVector|Thread|Imt|Net|MultiProc|Graf|Gpad|TreePlayer|Rint|Gdml|XMLIO)\.so")),
    ("support", re.compile(r"lib(tbb|boost_|fmt|FairLogger)\w*\.so")),
    ("system", re.compile(r"(libc\.so|libm\.so|libstdc\+\+\.so|libgcc_s|ld-linux|libpthread|libdl|vgpreload|\?\?\?)")),
]


def is_plt(key):
    """Callgrind names a PLT stub by its address, in the unknown object '???'."""
    ob, fn = key
    return (ob or "???") == "???" or (fn or "").startswith("0x")


# The layer of a function by its source file. This also works when several
# layers are linked into one binary (T2b and later).
FILE_LAYERS = [
    ("VecGeom", re.compile(r"lib(vecgeom|tgeo2vecgeom)\w*\.so")),
    ("VGM", re.compile(r"lib(BaseVGM|ClhepVGM|Geant4GM|RootGM|XmlVGM)\w*\.so")),
    ("VecGeom", re.compile(r"/extern/(vecgeom|veccore|tgeo2vecgeom)/")),
    ("VGM", re.compile(r"/extern/vgm/")),
    ("geant4_vmc", re.compile(r"/extern/geant4_vmc/")),
    ("Geant4", re.compile(r"/extern/geant4/")),
    ("VMC", re.compile(r"/extern/vmc/")),
    ("FairRoot", re.compile(r"/extern/o2/Detectors/3rdparty/FairRootBase/")),
    ("O2", re.compile(r"/extern/o2/|/app/")),
    ("ROOT", re.compile(r"/extern/root/")),
    ("Pythia8", re.compile(r"/extern/pythia8/")),
    ("support", re.compile(r"/extern/(tbb|boost|fmt|fairlogger)/")),
    # Inline code of the C++ standard library, compiled into the caller's object.
    ("system", re.compile(r"^/usr/(include|lib/gcc)/")),
]


def layer_of_key(key, files):
    """Layer of a function key (object, function): by namespace for templates of
    VecGeom (instantiated in the caller's file), then by source file, else by object."""
    if re.match(r"(?:\S+ )?vecgeom::", key[1] or ""):
        return "VecGeom"
    f = files.get(key)
    if f:
        for name, pat in FILE_LAYERS:
            if pat.search(f):
                return name
    return layer_of(key[0])


def layer_of(obj):
    base = obj or "???"
    if base == "???":
        return "PLT stubs"
    for name, pat in LAYERS:
        if pat.search(base):
            return name
    return "other:" + os.path.basename(base)


def parse(path):
    """Self cost per (object, function) and call counts per (caller, callee) function."""
    names = {}
    self_cost = collections.Counter()
    calls = collections.Counter()
    files = {}
    ob = fn = cob = cfn = fl = None
    pending_call = None
    events = None

    def resolve(kind, value):
        m = re.match(r"\((\d+)\)(?:\s+(.*))?", value)
        if not m:
            return value
        if m.group(2):
            names[(kind, m.group(1))] = m.group(2)
        return names.get((kind, m.group(1)), value)

    with open(path, errors="replace") as f:
        for line in f:
            if line.startswith("events:"):
                events = line.split()[1:]
            elif line.startswith("ob="):
                ob = resolve("ob", line[3:].strip())
            elif line.startswith("cob="):
                cob = resolve("ob", line[4:].strip())
            elif line.startswith("fl="):
                fl = resolve("fl", line[3:].strip())
            elif line.startswith(("fi=", "fe=")):
                resolve("fl", line[3:].strip())
            elif line.startswith("cfi=") or line.startswith("cfl="):
                resolve("fl", line[4:].strip())
            elif line.startswith("fn="):
                fn = resolve("fn", line[3:].strip())
                files.setdefault((ob, fn), fl)
                cob = cfn = None
            elif line.startswith("cfn="):
                cfn = resolve("fn", line[4:].strip())
            elif line.startswith("calls="):
                pending_call = int(line.split("=")[1].split()[0])
            elif line and (line[0].isdigit() or line[0] in "+-*"):
                if pending_call is not None:
                    calls[((ob, fn), (cob or ob, cfn))] += pending_call
                    pending_call = None
                    cob = None
                else:
                    parts = line.split()
                    self_cost[(ob, fn)] += int(parts[1]) if len(parts) > 1 else 0
    return self_cost, calls, files, (events or ["Ir"])[0]


def short(fn, width=70):
    fn = re.sub(r"\(.*", "", fn or "?")
    return fn if len(fn) <= width else fn[: width - 3] + "..."


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("profile")
    ap.add_argument("--minus", help="subtract this profile, for example a run with fewer events")
    ap.add_argument("--per", type=float, default=1.0, help="divide all numbers by this, for example the event difference")
    ap.add_argument("--top", type=int, default=15)
    ap.add_argument("--json")
    args = ap.parse_args()

    self_cost, calls, files, event = parse(args.profile)
    label = os.path.basename(args.profile)
    if args.minus:
        s2, c2, f2, _ = parse(args.minus)
        for k, v in f2.items():
            files.setdefault(k, v)
        self_cost.subtract(s2)
        calls.subtract(c2)
        label += " minus " + os.path.basename(args.minus)
    per = args.per

    layer_cost = collections.Counter()
    for (ob, fn), c in self_cost.items():
        layer_cost[layer_of_key((ob, fn), files)] += c
    # Look through PLT stubs: a call to a stub counts as calls to the targets of the stub.
    plt_targets = collections.defaultdict(collections.Counter)
    for (caller, callee), n in calls.items():
        if is_plt(caller) and n > 0:
            plt_targets[caller][callee] += n
    resolved = collections.Counter()
    plt_calls = collections.Counter()
    for (caller, callee), n in calls.items():
        if is_plt(caller):
            continue
        if is_plt(callee) and plt_targets[callee]:
            out = plt_targets[callee]
            total_out = sum(out.values())
            for target, m in out.items():
                resolved[(caller, target)] += n * m / total_out
                plt_calls[(layer_of_key(caller, files), layer_of_key(target, files))] += n * m / total_out
        else:
            resolved[(caller, callee)] += n

    crossings = collections.Counter()
    pair_calls = collections.Counter()
    for ((ob, fn), (cob, cfn)), n in resolved.items():
        a, b = layer_of_key((ob, fn), files), layer_of_key((cob, cfn), files)
        if a != b:
            crossings[(a, b)] += n
            pair_calls[(a, b, short(fn, 50), short(cfn, 50))] += n
    total = sum(layer_cost.values())

    print(f"{label}: self cost per layer ({event}), divided by {per:g}")
    for name, c in layer_cost.most_common():
        if abs(c) / per >= 1:
            print(f"  {name:22s} {c / per:16,.0f}  {100 * c / total:6.2f} %")
    print(f"  {'total':22s} {total / per:16,.0f}")
    print(f"top {args.top} functions by self cost")
    for (ob, fn), c in self_cost.most_common(args.top):
        print(f"  {100 * c / total:6.2f} %  {layer_of_key((ob, fn), files):10s} {short(fn)}")
    print("calls across layer boundaries (caller -> callee)")
    for (a, b), n in crossings.most_common(args.top):
        if n / per >= 1:
            print(f"  {a:12s} -> {b:12s} {n / per:14,.0f}")
    print("of these, calls through a PLT stub (calls into another shared object)")
    for (a, b), n in plt_calls.most_common(args.top):
        if n / per >= 1:
            print(f"  {a:12s} -> {b:12s} {n / per:14,.0f}")
    print(f"top {args.top} boundary call pairs")
    for (a, b, f1, f2), n in pair_calls.most_common(args.top):
        print(f"  {n / per:12,.0f}  {a}:{f1}  ->  {b}:{f2}")
    if args.json:
        with open(args.json, "w") as f:
            json.dump({"per": per, "total": total / per,
                       "self": {k: v / per for k, v in layer_cost.items()},
                       "crossings": {f"{a}->{b}": n / per for (a, b), n in crossings.items()},
                       "plt_calls": {f"{a}->{b}": n / per for (a, b), n in plt_calls.items()},
                       "top_pairs": [[a, b, f1, f2, n / per] for (a, b, f1, f2), n in pair_calls.most_common(50)],
                       "top_functions": [[layer_of_key((ob, fn), files), short(fn, 200), c / per] for (ob, fn), c in self_cost.most_common(50)]},
                      f, indent=1)


if __name__ == "__main__":
    sys.exit(main())
