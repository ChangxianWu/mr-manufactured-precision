# -*- coding: utf-8 -*-
"""Generate manuscript tables from the frozen pair-level results."""
import csv
import math
import statistics
import sys
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = Path("D:/Horzontal_project/20B3/work/crossdomain")
OUT = ROOT / "output"
PAIRS = OUT / "pairs.csv"
FIRST_ENDPOINTS = {"SEPSIS", "PNEUMONIA", "GASTROENT"}


def read_rows():
    with open(PAIRS, encoding="utf-8", newline="") as fh:
        return [r for r in csv.DictReader(fh, delimiter="\t") if int(r["k"]) >= 2]


def rank(v):
    o = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(o):
        j = i
        while j + 1 < len(o) and v[o[j + 1]] == v[o[i]]:
            j += 1
        a = (i + j) / 2.0 + 1
        for k in range(i, j + 1):
            r[o[k]] = a
        i = j + 1
    return r


def spearman(x, y):
    rx, ry = rank(x), rank(y)
    n = len(x)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    return num / (dx * dy) if dx and dy else float("nan")


def criteria(rows):
    f = [float(r["f"]) for r in rows]
    share = [float(r["max_window_share"]) for r in rows]
    keff = [float(r["k_eff"]) / float(r["k"]) for r in rows]
    kk = [float(r["k_kept"]) for r in rows]
    p1 = spearman(share, f)
    p3 = spearman([math.log10(float(r["k"])) for r in rows],
                  [math.log10(max(x, 1e-9)) for x in keff])
    p4 = spearman(kk, f)
    hi = [v for v, s in zip(f, share) if s >= 0.8]
    lo = [v for v, s in zip(f, share) if s < 0.8]
    p5 = [(float(r["f"]), float(r["f_block"])) for r in rows if r["f_block"] not in ("", None)]
    wins = sum(1 for a, b in p5 if b > a)
    return {
        "n": len(rows), "p1": p1, "p3": p3, "p4": p4,
        "hi_med": statistics.median(hi) if hi else float("nan"), "hi_n": len(hi),
        "lo_med": statistics.median(lo) if lo else float("nan"), "lo_n": len(lo),
        "p5_wins": wins, "p5_n": len(p5),
        "p5_gain": statistics.median([b - a for a, b in p5]) if p5 else float("nan"),
    }


def main():
    rows = read_rows()
    first = [r for r in rows if r["endpoint"] in FIRST_ENDPOINTS]
    extra = [r for r in rows if r["endpoint"] not in FIRST_ENDPOINTS]

    by = {}
    for r in rows:
        by.setdefault(r["trait"], []).append(r)
    order = sorted(by, key=lambda t: statistics.median([float(x["f"]) for x in by[t]]))
    t1 = ["| Trait | n pairs | k (median) | blocks after clumping | k_eff/k | f (median) | f (range) | largest 10 Mb window | median F |",
          "|---|---|---|---|---|---|---|---|---|"]
    for t in order:
        rs = by[t]
        ks = [float(r["k"]) for r in rs]
        kk = [float(r["k_kept"]) for r in rs]
        ke = [float(r["k_eff"]) / float(r["k"]) for r in rs]
        fs = [float(r["f"]) for r in rs]
        sh = max(float(r["max_window_share"]) for r in rs)
        mf = statistics.median([float(r["median_F"]) for r in rs])
        t1.append("| %s | %d | %.0f | %.0f | %.3f | %.3f | %.3f–%.3f | %.2f | %.0f |"
                  % (t, len(rs), statistics.median(ks), statistics.median(kk),
                     statistics.median(ke), statistics.median(fs), min(fs), max(fs), sh, mf))

    allc, fc, ec = criteria(rows), criteria(first), criteria(extra)
    t2 = ["| Criterion | All (n=%d) | Original 3 endpoints (n=%d) | Added 8 endpoints (n=%d) |"
          % (allc["n"], fc["n"], ec["n"]),
          "|---|---|---|---|",
          "| P1 concentration ↔ f (Spearman) | %.3f | %.3f | %.3f |"
          % (allc["p1"], fc["p1"], ec["p1"]),
          "| P2 f, single-locus vs multi-locus | %.3f (n=%d) vs %.3f (n=%d) | %.3f (n=%d) vs %.3f (n=%d) | %.3f (n=%d) vs %.3f (n=%d) |"
          % (allc["hi_med"], allc["hi_n"], allc["lo_med"], allc["lo_n"],
             fc["hi_med"], fc["hi_n"], fc["lo_med"], fc["lo_n"],
             ec["hi_med"], ec["hi_n"], ec["lo_med"], ec["lo_n"]),
          "| P3 log k ↔ log(k_eff/k) | %.3f | %.3f | %.3f |"
          % (allc["p3"], fc["p3"], ec["p3"]),
          "| P4 blocks ↔ f | %.3f | %.3f | %.3f |"
          % (allc["p4"], fc["p4"], ec["p4"]),
          "| P5 single-block f > whole-trait f | %d/%d (+%.3f) | %d/%d (+%.3f) | %d/%d (+%.3f) |"
          % (allc["p5_wins"], allc["p5_n"], allc["p5_gain"],
             fc["p5_wins"], fc["p5_n"], fc["p5_gain"],
             ec["p5_wins"], ec["p5_n"], ec["p5_gain"]),
          ]

    txt = "\n".join([
        "# companion     11_companion_tables.py        ",
        "",
        "## Table 1. Metabolite traits, ordered by median f",
        "",
    ] + t1 + ["", "## Table 2. Pre-declared criteria, by stratum", ""] + t2 + [""])
    (OUT / "tables.md").write_text(txt, encoding="utf-8", newline="")
    print(txt)


if __name__ == "__main__":
    main()
