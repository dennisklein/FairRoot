#!/usr/bin/env python3
"""Run the FV0 benchmark for one build variant.

For each seed, the script runs fv0sim, measures time and memory, and
writes the physics observables with fv0dump. Results go to
<out>/<label>/:
  runs.csv              one line per run (times, memory, output hash)
  summary.json          fixed cost and time per event
  physics/seed<N>.*.csv observables for compare.py
  callgrind.<N>.out     with --callgrind, for layers.py

Usage:
  run.py --prefix <install prefix> --label <name> --seeds 1-5 --events 1000
         [--callgrind] [--out results]
"""
import argparse
import csv
import hashlib
import json
import os
import re
import resource
import shutil
import statistics
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))

# The benchmark: charged pions into the FV0 acceptance, nominal L3 field.
SIM_ARGS = ["-m", "PIPE", "FV0", "-e", "TGeant4", "-g", "boxgen"]
FIELD_KG = "5"
KEY_VALUES = ";".join([
    "BoxGun.pdg=211",
    "BoxGun.number=10",
    "BoxGun.eta[0]=2.2", "BoxGun.eta[1]=5.1",
    "BoxGun.prange[0]=0.5", "BoxGun.prange[1]=10",
    "BoxGun.phirange[0]=0", "BoxGun.phirange[1]=360",
    "align-geom.mDetectors=none",
])


def seeds(spec):
    out = []
    for part in spec.split(","):
        a, _, b = part.partition("-")
        out += list(range(int(a), int(b or a) + 1))
    return out


def env_for(prefix):
    """Environment of an installed stack, as thisroot.sh and geant4.sh set it."""
    script = (f"source {prefix}/bin/thisroot.sh >/dev/null && "
              f"source {prefix}/bin/geant4.sh >/dev/null && env -0")
    out = subprocess.run(["bash", "-c", script], capture_output=True, check=True).stdout
    env = dict(kv.split("=", 1) for kv in out.decode().split("\0") if "=" in kv)
    env["LD_LIBRARY_PATH"] = f"{prefix}/lib:" + env.get("LD_LIBRARY_PATH", "")
    env["O2_ROOT"] = prefix
    env["VMCWORKDIR"] = f"{prefix}/share"
    env["ROOT_INCLUDE_PATH"] = f"{prefix}/include"
    env["PYTHIA8DATA"] = f"{prefix}/share/Pythia8/xmldoc"
    # Static TLS reserve for the dlopened Geant4 (as SimSetup/StaticTLSReserve.h sets it). Without it,
    # fv0sim restarts itself, and valgrind would profile only the first process.
    env["GLIBC_TUNABLES"] = "glibc.rtld.optional_static_tls=32768"
    return env


def run_one(prefix, env, workdir, seed, events, callgrind, field, extra_kv):
    cmd = [f"{prefix}/bin/fv0sim", *SIM_ARGS, "--field", field, "-n", str(events), "--seed", str(seed),
           "-o", f"seed{seed}", "--configKeyValues", ";".join(filter(None, [KEY_VALUES, extra_kv]))]
    if callgrind:
        cmd = ["valgrind", "--tool=callgrind", f"--callgrind-out-file=callgrind.{seed}.out",
               "--dump-instr=no", "--collect-jumps=no", *cmd]
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    t0 = time.perf_counter()
    with open(os.path.join(workdir, f"seed{seed}.log"), "w") as log:
        rc = subprocess.run(cmd, cwd=workdir, env=env, stdout=log, stderr=subprocess.STDOUT).returncode
    wall = time.perf_counter() - t0
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    text = open(os.path.join(workdir, f"seed{seed}.log"), errors="replace").read()
    if rc != 0:
        sys.exit(f"fv0sim failed for seed {seed}, see {workdir}/seed{seed}.log")
    init = re.search(r"Init: Real time ([\d.]+) s", text)
    run = re.search(r"Run: Real time ([\d.]+) s, CPU time ([\d.]+)", text)
    return {
        "seed": seed, "events": events, "wall_s": round(wall, 3),
        "cpu_s": round((after.ru_utime - before.ru_utime) + (after.ru_stime - before.ru_stime), 3),
        "init_s": float(init.group(1)) if init else "",
        "run_s": float(run.group(1)) if run else "",
        "run_cpu_s": float(run.group(2)) if run else "",
        "maxrss_mb": round(after.ru_maxrss / 1024, 1),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--seeds", default="1-5")
    ap.add_argument("--events", type=int, default=1000)
    ap.add_argument("--fixed-runs", type=int, default=3,
                    help="1-event runs that measure the fixed cost (start-up, physics tables)")
    ap.add_argument("--callgrind", action="store_true")
    ap.add_argument("--field", default=FIELD_KG,
                    help="field in kG; change it only to test the power of the physics gate")
    ap.add_argument("--key-values", default="",
                    help="more --configKeyValues; use it only to test the power of the physics gate")
    ap.add_argument("--out", default=os.path.join(HERE, "..", "results"))
    args = ap.parse_args()

    prefix = os.path.abspath(args.prefix)
    out = os.path.abspath(os.path.join(args.out, args.label))
    work = os.path.join(out, "work")
    phys = os.path.join(out, "physics")
    os.makedirs(work, exist_ok=True)
    os.makedirs(phys, exist_ok=True)
    env = env_for(prefix)

    rows = []
    for s in seeds(args.seeds):
        r = run_one(prefix, env, work, s, args.events, args.callgrind, args.field, args.key_values)
        dump = os.path.join(phys, f"seed{s}")
        subprocess.run([f"{prefix}/bin/fv0dump", os.path.join(work, f"seed{s}.root"), dump],
                       env=env, check=True)
        r["hits_sha1"] = hashlib.sha1(open(dump + ".hits.csv", "rb").read()).hexdigest()[:12]
        if args.callgrind:
            shutil.move(os.path.join(work, f"callgrind.{s}.out"), os.path.join(out, f"callgrind.{s}.out"))
        rows.append(r)
        print(f"{args.label} seed {s}: {r}", flush=True)

    with open(os.path.join(out, "runs.csv" if not args.callgrind else "runs-callgrind.csv"), "w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    if args.callgrind or args.fixed_runs < 1:
        return

    # Fixed cost: runs with 1 event. Per-event time = (wall - fixed) / (events - 1).
    fixed = []
    for i in range(args.fixed_runs):
        fixed.append(run_one(prefix, env, work, seeds(args.seeds)[0], 1, False, args.field, args.key_values)["wall_s"])
    fixed_s = statistics.median(fixed)
    per_event = [(r["wall_s"] - fixed_s) / (r["events"] - 1) * 1e3 for r in rows]
    summary = {
        "label": args.label, "events": args.events, "seeds": len(rows),
        "fixed_s": fixed_s, "fixed_runs_s": fixed,
        "ms_per_event": statistics.mean(per_event),
        "ms_per_event_err": statistics.stdev(per_event) / len(per_event) ** 0.5 if len(per_event) > 1 else 0.0,
        "maxrss_mb": max(r["maxrss_mb"] for r in rows),
    }
    with open(os.path.join(out, "summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
    print(f"{args.label}: fixed {fixed_s:.2f} s, {summary['ms_per_event']:.2f} "
          f"+- {summary['ms_per_event_err']:.2f} ms per event", flush=True)


if __name__ == "__main__":
    main()
