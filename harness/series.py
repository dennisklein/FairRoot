#!/usr/bin/env python3
"""Measure several builds in one interleaved series.

The machine drifts by several percent between sessions. That is more than
the effect of one transformation. This script thus runs the builds round
robin: seed 1 on every build, then seed 2 on every build, and so on, plus
the 1-event runs for the fixed cost. Every build sees the same drift.

Usage:
  series.py --out results/series-<name> --seeds 101-105 --events 1000 \\
            baseline=/home/user/build/stage t1a=/home/user/build-isa/stage ...

Writes <out>/<label>/runs.csv and physics/, and <out>/summary.csv with the
CPU time per event (CPU time minus the fixed cost, divided by events - 1).
"""
import argparse
import csv
import hashlib
import os
import statistics
import subprocess

import run as harness


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("builds", nargs="+", help="label=prefix")
    ap.add_argument("--out", required=True)
    ap.add_argument("--seeds", default="101-105")
    ap.add_argument("--events", type=int, default=1000)
    ap.add_argument("--fixed-runs", type=int, default=3)
    args = ap.parse_args()

    builds = [b.split("=", 1) for b in args.builds]
    envs = {label: harness.env_for(os.path.abspath(prefix)) for label, prefix in builds}
    rows = {label: [] for label, _ in builds}
    fixed = {label: [] for label, _ in builds}
    seed_list = harness.seeds(args.seeds)

    for label, _ in builds:
        for sub in ("work", "physics"):
            os.makedirs(os.path.join(args.out, label, sub), exist_ok=True)

    def one(label, prefix, seed, events):
        work = os.path.join(args.out, label, "work")
        return harness.run_one(os.path.abspath(prefix), envs[label], work, seed, events,
                               False, harness.FIELD_KG, "")

    for i in range(max(len(seed_list), args.fixed_runs)):
        for label, prefix in builds:
            if i < args.fixed_runs:
                fixed[label].append(one(label, prefix, seed_list[0], 1))
            if i < len(seed_list):
                s = seed_list[i]
                r = one(label, prefix, s, args.events)
                dump = os.path.join(args.out, label, "physics", f"seed{s}")
                subprocess.run([f"{prefix}/bin/fv0dump", os.path.join(args.out, label, "work", f"seed{s}.root"), dump],
                               env=envs[label], check=True)
                r["hits_sha1"] = hashlib.sha1(open(dump + ".hits.csv", "rb").read()).hexdigest()[:12]
                rows[label].append(r)
                print(f"{label} seed {s}: cpu {r['cpu_s']} s", flush=True)

    summary = []
    for label, _ in builds:
        with open(os.path.join(args.out, label, "runs.csv"), "w") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[label][0]))
            w.writeheader()
            w.writerows(rows[label])
        fixed_cpu = statistics.median(r["cpu_s"] for r in fixed[label])
        per = [(r["cpu_s"] - fixed_cpu) / (r["events"] - 1) * 1e3 for r in rows[label]]
        summary.append({"label": label, "cpu_ms_per_event": round(statistics.mean(per), 3),
                        "err": round(statistics.stdev(per) / len(per) ** 0.5, 3) if len(per) > 1 else 0,
                        "fixed_cpu_s": round(fixed_cpu, 3),
                        "maxrss_mb": max(r["maxrss_mb"] for r in rows[label])})
    base = summary[0]["cpu_ms_per_event"]
    with open(os.path.join(args.out, "summary.csv"), "w") as f:
        w = csv.DictWriter(f, fieldnames=list(summary[0]) + ["vs_first_pct"])
        w.writeheader()
        for s in summary:
            s["vs_first_pct"] = round(100 * (s["cpu_ms_per_event"] / base - 1), 1)
            w.writerow(s)
            print(f"{s['label']:12s} {s['cpu_ms_per_event']:7.2f} +- {s['err']:.2f} ms/event (cpu), "
                  f"{s['vs_first_pct']:+.1f} %, fixed {s['fixed_cpu_s']} s, {s['maxrss_mb']} MB")


if __name__ == "__main__":
    main()
