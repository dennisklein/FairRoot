#!/usr/bin/env python3
"""Physics-equality gate of the flatsim harness.

Compares the FV0 observables of a candidate build with a reference build.
Each input is a directory with the fv0dump CSV files of several runs
(seed<N>.hits.csv, seed<N>.events.csv). The two samples must come from
different seeds, so that they are statistically independent.

Event observables are independent from event to event, so the script
uses the asymptotic two-sample KS test for them. Hit observables are not
independent: one particle gives several hits in neighbouring cells, and a
shower gives many. A test that treats each hit as independent gives
p-values that are much too small. For hit observables, the script thus
uses a permutation test over whole events: it mixes the events of the two
samples, splits them again at random, and counts how often the statistic
(KS distance or chi-square) of a random split is at least as large as the
observed one.

The gate passes when no observable differs significantly: every p-value
must be above alpha / (number of tests) (Bonferroni correction).

Usage: compare.py <reference-dir> <candidate-dir> [--alpha 0.01] [--json out.json]
"""
import argparse
import csv
import glob
import json
import math
import os
import re
import sys

import numpy as np


def read(directory, kind):
    """Rows of all seeds. The directory is a result of run.py or its physics/ subdirectory."""
    if os.path.isdir(os.path.join(directory, "physics")):
        directory = os.path.join(directory, "physics")
    paths = sorted(glob.glob(os.path.join(directory, f"seed*.{kind}.csv")))
    if not paths:
        sys.exit(f"no seed*.{kind}.csv in {directory}")
    rows = []
    for path in paths:
        seed = int(re.search(r"seed(\d+)\.", os.path.basename(path)).group(1))
        with open(path) as f:
            for r in csv.DictReader(f):
                row = {k: float(v) for k, v in r.items()}
                row["key"] = (seed, int(row["event"]))
                rows.append(row)
    return rows


def ks_2samp(a, b):
    """Two-sample Kolmogorov-Smirnov test, asymptotic p-value."""
    a, b = sorted(a), sorted(b)
    n, m = len(a), len(b)
    if n == 0 or m == 0:
        sys.exit("KS test with an empty sample")
    i = j = 0
    d = 0.0
    while i < n and j < m:
        x = min(a[i], b[j])
        while i < n and a[i] <= x:
            i += 1
        while j < m and b[j] <= x:
            j += 1
        d = max(d, abs(i / n - j / m))
    ne = n * m / (n + m)
    lam = (math.sqrt(ne) + 0.12 + 0.11 / math.sqrt(ne)) * d
    p = 0.0
    for k in range(1, 101):
        term = 2 * (-1) ** (k - 1) * math.exp(-2 * k * k * lam * lam)
        p += term
        if abs(term) < 1e-12:
            break
    return d, min(max(p, 0.0), 1.0)


def gammaq(a, x):
    """Regularized upper incomplete gamma function Q(a, x)."""
    if x <= 0:
        return 1.0
    gln = math.lgamma(a)
    if x < a + 1:  # series for P
        s = term = 1.0 / a
        ap = a
        for _ in range(1000):
            ap += 1
            term *= x / ap
            s += term
            if abs(term) < abs(s) * 1e-15:
                break
        return 1.0 - s * math.exp(-x + a * math.log(x) - gln)
    b = x + 1 - a  # continued fraction for Q
    c = 1e300
    d = 1 / b
    h = d
    for i in range(1, 1000):
        an = -i * (i - a)
        b += 2
        d = an * d + b
        d = 1e-300 if abs(d) < 1e-300 else d
        c = b + an / c
        c = 1e-300 if abs(c) < 1e-300 else c
        d = 1 / d
        delta = d * c
        h *= delta
        if abs(delta - 1) < 1e-15:
            break
    return math.exp(-x + a * math.log(x) - gln) * h


def chi2_homogeneity(ca, cb, min_expected=5.0):
    """Chi-square test that two histograms come from the same distribution."""
    keys = sorted(set(ca) | set(cb))
    na, nb = sum(ca.values()), sum(cb.values())
    chi2, ndf = 0.0, -1
    for k in keys:
        a, b = ca.get(k, 0), cb.get(k, 0)
        tot = a + b
        ea, eb = tot * na / (na + nb), tot * nb / (na + nb)
        if ea < min_expected or eb < min_expected:
            continue
        chi2 += (a - ea) ** 2 / ea + (b - eb) ** 2 / eb
        ndf += 1
    if ndf < 1:
        return float("nan"), 0, 1.0
    return chi2, ndf, gammaq(ndf / 2, chi2 / 2)


def event_index(ref_e, can_e, ref_h, can_h):
    """Pooled event numbering, label of each event (True = reference), and event of each hit."""
    keys = [("r",) + r["key"] for r in ref_e] + [("c",) + r["key"] for r in can_e]
    index = {k: i for i, k in enumerate(keys)}
    is_ref = np.array([k[0] == "r" for k in keys])
    hit_event = np.array([index[("r",) + h["key"]] for h in ref_h] + [index[("c",) + h["key"]] for h in can_h])
    return is_ref, hit_event


def permutations(is_ref, n, rng, batch):
    """Random relabelings of the events that keep the number of reference events."""
    for start in range(0, n, batch):
        b = min(batch, n - start)
        yield np.argsort(rng.random((b, is_ref.size)), axis=1) < is_ref.sum()


def perm_ks(values, hit_event, is_ref, n, rng, batch=64):
    """KS distance of hit values between the samples, with an event permutation p-value."""
    order = np.argsort(values, kind="stable")
    v = values[order]
    ev = hit_event[order]
    last = np.append(v[1:] != v[:-1], True)  # compare the CDFs only after the last of equal values

    def dist(labels):  # labels: (b, events) bool
        lab = labels[:, ev]
        ca = np.cumsum(lab, axis=1)
        cb = np.arange(1, lab.shape[1] + 1) - ca
        na, nb = ca[:, -1:], cb[:, -1:]
        return np.abs(ca / na - cb / nb)[:, last].max(axis=1)

    d0 = dist(is_ref[None, :])[0]
    exceed = sum(int((dist(p) >= d0 - 1e-12).sum()) for p in permutations(is_ref, n, rng, batch))
    return d0, (exceed + 1) / (n + 1)


def perm_chi2(categories, hit_event, is_ref, n, rng, min_expected=5.0, batch=256):
    """Chi-square homogeneity of hit categories, with an event permutation p-value."""
    cats, cat_index = np.unique(categories, return_inverse=True)
    per_event = np.zeros((is_ref.size, cats.size))
    np.add.at(per_event, (hit_event, cat_index), 1)
    total = per_event.sum(axis=0)

    def chi2(labels):
        a = labels.astype(float) @ per_event
        b = total - a
        na, nb = a.sum(axis=1, keepdims=True), b.sum(axis=1, keepdims=True)
        ea, eb = total * na / (na + nb), total * nb / (na + nb)
        keep = (total * min(na.min(), nb.min()) / (na + nb).max()) >= min_expected
        return (((a - ea) ** 2 / ea + (b - eb) ** 2 / eb)[:, keep]).sum(axis=1), int(keep.sum()) - 1

    c0, ndf = chi2(is_ref[None, :])
    exceed = sum(int((chi2(p)[0] >= c0[0] - 1e-9).sum()) for p in permutations(is_ref, n, rng, batch))
    return c0[0], ndf, (exceed + 1) / (n + 1)


def counts(rows, key):
    c = {}
    for r in rows:
        c[int(r[key])] = c.get(int(r[key]), 0) + 1
    return c


def mean_err(v):
    n = len(v)
    if n < 2:
        return float("nan"), float("nan")
    m = sum(v) / n
    s = math.sqrt(sum((x - m) ** 2 for x in v) / (n - 1))
    return m, s / math.sqrt(n)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("reference")
    ap.add_argument("candidate")
    ap.add_argument("--alpha", type=float, default=0.01)
    ap.add_argument("--permutations", type=int, default=4000)
    ap.add_argument("--json")
    args = ap.parse_args()
    rng = np.random.default_rng(12345)

    ref_e, can_e = read(args.reference, "events"), read(args.candidate, "events")
    ref_h, can_h = read(args.reference, "hits"), read(args.candidate, "hits")
    tests = []

    for key in ["nhits", "edep", "ntracks", "nrefs"]:
        a, b = [r[key] for r in ref_e], [r[key] for r in can_e]
        d, p = ks_2samp(a, b)
        tests.append({"observable": f"event.{key}", "test": "KS", "stat": d, "p": p,
                      "ref_mean": mean_err(a), "cand_mean": mean_err(b)})
    for r in ref_h + can_h:
        r["r"] = math.hypot(r["x"], r["y"])
    is_ref, hit_event = event_index(ref_e, can_e, ref_h, can_h)
    for r in ref_h + can_h:
        r["dphi_prim"] = r["dphi"] if r["primary"] == 1 else None
        r["dr_prim"] = r["dr"] if r["primary"] == 1 else None  # primaries: the bending is clean
    for key in ["edep", "t", "r", "z", "e", "dphi_prim", "dr_prim"]:
        sel = np.array([r[key] is not None for r in ref_h + can_h])
        a = [r[key] for r in ref_h if r[key] is not None]
        b = [r[key] for r in can_h if r[key] is not None]
        d, p = perm_ks(np.array(a + b), hit_event[sel], is_ref, args.permutations, rng)
        tests.append({"observable": f"hit.{key}", "test": "KS perm", "stat": d, "p": p,
                      "ref_mean": mean_err(a), "cand_mean": mean_err(b)})
    for key in ["cell", "pdg", "primary"]:
        cats = np.array([int(r[key]) for r in ref_h + can_h])
        chi2, ndf, p = perm_chi2(cats, hit_event, is_ref, args.permutations, rng)
        tests.append({"observable": f"hit.{key}", "test": f"chi2 perm/ndf={ndf}", "stat": chi2, "p": p})

    threshold = args.alpha / len(tests)
    passed = all(t["p"] > threshold for t in tests)
    print(f"reference: {len(ref_e)} events, {len(ref_h)} hits; "
          f"candidate: {len(can_e)} events, {len(can_h)} hits")
    print(f"{'observable':16s} {'test':14s} {'stat':>10s} {'p':>9s}  ref mean / cand mean")
    for t in tests:
        flag = "" if t["p"] > threshold else "  <-- FAIL"
        means = ""
        if "ref_mean" in t:
            means = f"  {t['ref_mean'][0]:.5g}±{t['ref_mean'][1]:.2g} / {t['cand_mean'][0]:.5g}±{t['cand_mean'][1]:.2g}"
        print(f"{t['observable']:16s} {t['test']:14s} {t['stat']:10.4g} {t['p']:9.3g}{means}{flag}")
    print(f"gate: {'PASS' if passed else 'FAIL'} (pass needs every p > {threshold:.2g})")
    if args.json:
        with open(args.json, "w") as f:
            json.dump({"passed": passed, "threshold": threshold, "tests": tests}, f, indent=1)
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
