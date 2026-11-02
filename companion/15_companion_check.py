# -*- coding: utf-8 -*-
"""Recompute every reported number from the result files and compare with the manuscript text and docx."""
import csv
import gzip
import json
import math
import re
import statistics
import sys
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = Path("D:/Horzontal_project/20B3/work/crossdomain")
OUT = ROOT / "output"
MD = Path("D:/Horzontal_project/20B3/   /20B3_companion_manuscript_v0.3.md")
DOCX = OUT / "20B3_companion_manuscript_v0.3.docx"

FAIL = []


def check(label, ok, detail=""):
    print(("  [OK] " if ok else "  [FAIL] ") + label + (("  " + detail) if detail else ""))
    if not ok:
        FAIL.append(label)


def read_pairs(path):
    with open(path, encoding="utf-8", newline="") as fh:
        return [r for r in csv.DictReader(fh, delimiter="\t") if int(r["k"]) >= 2]


def norm_p(b, se):
    return math.erfc(abs(b / se) / math.sqrt(2.0)) if se > 0 else 1.0


def bh_count(rows, col_b, col_se, alpha=0.05):
    ps = [norm_p(float(r[col_b]), float(r[col_se])) for r in rows]
    n = len(ps)
    order = sorted(range(n), key=lambda i: ps[i])
    th = 0
    for rank, i in enumerate(order, 1):
        if ps[i] <= alpha * rank / n:
            th = rank
    return th, min(ps)


def spearman(x, y):
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
    rx, ry = rank(x), rank(y)
    n = len(x)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = (sum((a - mx) ** 2 for a in rx)) ** 0.5
    dy = (sum((b - my) ** 2 for b in ry)) ** 0.5
    return num / (dx * dy) if dx and dy else float("nan")


def main():
    md = MD.read_text(encoding="utf-8")
    pairs = read_pairs(OUT / "pairs.csv")
    og = read_pairs(OUT / "pairs_opengwas.csv")
    jin = json.loads((OUT / "companion.json").read_text(encoding="utf-8"))

    f = [float(r["f"]) for r in pairs]
    share = [float(r["max_window_share"]) for r in pairs]
    kk = [float(r["k_kept"]) for r in pairs]
    p1 = spearman(share, f)
    p4 = spearman(kk, f)
    thr_un, min_un = bh_count(pairs, "b_unclumped", "se_unclumped")
    thr_cl, min_cl = bh_count(pairs, "b_clumped", "se_clumped")
    p5 = [(float(r["f"]), float(r["f_block"])) for r in pairs if r["f_block"] not in ("", None)]
    wins = sum(1 for a, b in p5 if b > a)
    fb = sorted(float(b) for a, b in p5)
    med_fb = statistics.median(fb)

    print("[15]     FinnGen %d   " % len(pairs))
    check("    110", len(pairs) == 110, str(len(pairs)))
    check("  clump    BH 90", thr_un == 90, str(thr_un))
    check("clumped    BH 0", thr_cl == 0, str(thr_cl))
    check("P1 rho=0.884", abs(p1 - 0.884) < 0.001, "%.3f" % p1)
    check("P4 rho=-0.969", abs(p4 + 0.969) < 0.001, "%.3f" % p4)
    check("P5 94/109", (wins, len(p5)) == (94, 109), "%d/%d" % (wins, len(p5)))
    check("f_block    0.758", abs(med_fb - 0.758) < 0.001, "%.3f" % med_fb)
    check("f    0.421", abs(statistics.median(f) - 0.421) < 0.001,
          "%.3f" % statistics.median(f))

    print("[15]    UKB %d   " % len(og))
    of = [float(r["f"]) for r in og]
    o1 = spearman([float(r["max_window_share"]) for r in og], of)
    o4 = spearman([float(r["k_kept"]) for r in og], of)
    othr_un, _ = bh_count(og, "b_unclumped", "se_unclumped")
    othr_cl, _ = bh_count(og, "b_clumped", "se_clumped")
    check("   20  ", len(og) == 20, str(len(og)))
    check("   f    0.440", abs(statistics.median(of) - 0.440) < 0.001,
          "%.3f" % statistics.median(of))
    check("   P1 0.950", abs(o1 - 0.950) < 0.001, "%.3f" % o1)
    check("   P4 -0.969", abs(o4 + 0.969) < 0.001, "%.3f" % o4)
    check("   14→1", (othr_un, othr_cl) == (14, 1), "%d→%d" % (othr_un, othr_cl))

    print("[15]     ")
    must = ["90 of 110", "81.8%", "0 of 110", "0.884", "0.969", "94 of 109",
            "0.758", "0.440", "0.950", "14", "Apolipoprotein A1"]
    for m in must:
        check("    %r" % m, m in md)
    stale = ["0.49–0.87", "ρ=0.922", "−0.957", "18/24", "0.583   "]
    for s in stale:
        check("     %r" % s, s not in md)

    print("[15] docx   md   ")
    try:
        from docx import Document
        d = Document(DOCX)
        text = "\n".join(p.text for p in d.paragraphs)
        cells = "\n".join(c.text for t in d.tables for row in t.rows for c in row.cells)
        allt = text + "\n" + cells
        for m in ["90 of 110", "0.884", "0.758", "0.440", "94 of 109"]:
            check("docx   %r" % m, m in allt)
        check("docx   markdown   ", "**" not in allt and "|" not in text)
        check("docx    = 3", len(d.tables) == 3, str(len(d.tables)))
        check("docx    = 1", len(d.inline_shapes) == 1, str(len(d.inline_shapes)))
    except Exception as exc:  # noqa: BLE001
        check("docx   ", False, str(exc)[:80])

    print()
    if FAIL:
        print("[15] FAILED %d  " % len(FAIL))
        sys.exit(1)
    print("[15] COMPANION_CHECK_OK")


if __name__ == "__main__":
    main()
