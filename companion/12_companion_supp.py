# -*- coding: utf-8 -*-
"""Build the supplementary workbook (per-pair, criteria, sensitivity, methods)."""
import csv
import statistics
import sys
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

from openpyxl import Workbook
from openpyxl.styles import Font, Alignment

ROOT = Path("D:/Horzontal_project/20B3/work/crossdomain")
OUT = ROOT / "output"
XLSX = OUT / "Supplementary_Tables_companion.xlsx"
FIRST = {"SEPSIS", "PNEUMONIA", "GASTROENT"}

PAIR_COLS = ["trait", "gcst", "endpoint", "k", "k_kept", "se_unclumped", "se_clumped",
             "k_eff", "f", "max_window_share", "b_unclumped", "b_clumped",
             "median_F", "n_F_lt10", "k_block", "f_block"]


def rows_of(path):
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def sheet(wb, title, header, data, widths=None):
    ws = wb.create_sheet(title)
    ws.append(header)
    for c in range(1, len(header) + 1):
        ws.cell(row=1, column=c).font = Font(bold=True)
        ws.cell(row=1, column=c).alignment = Alignment(horizontal="center")
    for r in data:
        ws.append(r)
    ws.freeze_panes = "A2"
    for i in range(1, len(header) + 1):
        ws.column_dimensions[ws.cell(row=1, column=i).column_letter].width = (
            widths[i - 1] if widths and i <= len(widths) else 14)
    return ws


def main():
    pairs = [r for r in rows_of(OUT / "pairs.csv") if int(r["k"]) >= 2]
    r01 = {r["endpoint"] + "|" + r["gcst"]: r
           for r in rows_of(OUT / "pairs_r01.csv") if int(r["k"]) >= 2}

    wb = Workbook()
    wb.remove(wb.active)

    # S1
    sheet(wb, "S1_pairs", PAIR_COLS,
          [[r.get(c, "") for c in PAIR_COLS] for r in pairs],
          widths=[24, 14, 13] + [11] * 13)

    # S2 criteria by stratum
    def crit(rs):
        def rank(v):
            o = sorted(range(len(v)), key=lambda i: v[i])
            out = [0.0] * len(v)
            i = 0
            while i < len(o):
                j = i
                while j + 1 < len(o) and v[o[j + 1]] == v[o[i]]:
                    j += 1
                a = (i + j) / 2.0 + 1
                for k in range(i, j + 1):
                    out[o[k]] = a
                i = j + 1
            return out

        def sp(x, y):
            rx, ry = rank(x), rank(y)
            n = len(x)
            mx, my = sum(rx) / n, sum(ry) / n
            num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
            dx = (sum((a - mx) ** 2 for a in rx)) ** 0.5
            dy = (sum((b - my) ** 2 for b in ry)) ** 0.5
            return num / (dx * dy) if dx and dy else float("nan")

        f = [float(r["f"]) for r in rs]
        sh = [float(r["max_window_share"]) for r in rs]
        p5 = [(float(r["f"]), float(r["f_block"])) for r in rs if r["f_block"]]
        return {
            "n": len(rs),
            "P1": sp(sh, f),
            "P2_single": statistics.median([v for v, s in zip(f, sh) if s >= 0.8]) if any(s >= 0.8 for s in sh) else "",
            "P2_multi": statistics.median([v for v, s in zip(f, sh) if s < 0.8]),
            "P3": sp([__import__("math").log10(float(r["k"])) for r in rs],
                     [__import__("math").log10(max(float(r["k_eff"]) / float(r["k"]), 1e-9)) for r in rs]),
            "P4": sp([float(r["k_kept"]) for r in rs], f),
            "P5_wins": sum(1 for a, b in p5 if b > a),
            "P5_n": len(p5),
            "P5_gain": statistics.median([b - a for a, b in p5]) if p5 else "",
        }

    first = [r for r in pairs if r["endpoint"] in FIRST]
    extra = [r for r in pairs if r["endpoint"] not in FIRST]
    s2 = []
    for label, rs in [("All", pairs), ("Original 3 endpoints", first),
                      ("Added 8 endpoints", extra)]:
        c = crit(rs)
        s2.append([label, c["n"], round(c["P1"], 3),
                   round(c["P2_single"], 3) if c["P2_single"] != "" else "",
                   round(c["P2_multi"], 3), round(c["P3"], 3), round(c["P4"], 3),
                   c["P5_wins"], c["P5_n"], round(c["P5_gain"], 3)])
    sheet(wb, "S2_criteria", ["stratum", "n_pairs", "P1_rho", "P2_single_locus_f",
                              "P2_multi_locus_f", "P3_rho", "P4_rho", "P5_wins",
                              "P5_evaluable", "P5_median_gain"], s2)

    # S3 threshold sensitivity
    s3 = []
    for r in pairs:
        k = r["endpoint"] + "|" + r["gcst"]
        if k not in r01:
            continue
        a = r01[k]
        s3.append([r["trait"], r["endpoint"], r["k"], r["k_kept"], round(float(r["f"]), 4),
                   a["k_kept"], round(float(a["f"]), 4),
                   round(float(a["f"]) - float(r["f"]), 4)])
    sheet(wb, "S3_r2_sensitivity",
          ["trait", "endpoint", "k_unclumped", "k_kept_r2_0.001", "f_r2_0.001",
           "k_kept_r2_0.01", "f_r2_0.01", "delta_f"], s3)

    # S4 per-trait summary
    by = {}
    for r in pairs:
        by.setdefault(r["trait"], []).append(r)
    s4 = []
    for t in sorted(by, key=lambda x: statistics.median([float(r["f"]) for r in by[x]])):
        rs = by[t]
        fs = [float(r["f"]) for r in rs]
        s4.append([t, len(rs),
                   round(statistics.median([float(r["k"]) for r in rs]), 1),
                   round(statistics.median([float(r["k_kept"]) for r in rs]), 1),
                   round(statistics.median([float(r["k_eff"]) / float(r["k"]) for r in rs]), 3),
                   round(statistics.median(fs), 3), round(min(fs), 3), round(max(fs), 3),
                   round(max(float(r["max_window_share"]) for r in rs), 2),
                   round(statistics.median([float(r["median_F"]) for r in rs]), 1)])
    sheet(wb, "S4_traits", ["trait", "n_pairs", "k_median", "k_kept_median",
                            "keff_over_k_median", "f_median", "f_min", "f_max",
                            "max_window_share", "median_F"], s4)

    # S5b external outcome resource (OpenGWAS / UKB)
    og = OUT / "pairs_opengwas.csv"
    if og.exists():
        ogr = [r for r in rows_of(og) if int(r["k"]) >= 2]
        sheet(wb, "S6_opengwas",
              ["trait", "endpoint", "k", "k_kept", "se_unclumped", "se_clumped",
               "k_eff", "f", "max_window_share", "b_unclumped", "b_clumped",
               "median_F", "n_F_lt10"],
              [[r[c] for c in ["trait", "endpoint", "k", "k_kept", "se_unclumped",
                               "se_clumped", "k_eff", "f", "max_window_share",
                               "b_unclumped", "b_clumped", "median_F", "n_F_lt10"]]
               for r in ogr], widths=[24] + [12] * 12)

    # S5 methods facts
    facts = [
        ("exposure", "Kettunen 2016 NMR metabolomics (GWAS Catalog), GRCh37, rsIDs present"),
        ("traits", "10 (ApoB, total cholesterol, GlycA, ApoA1, phenylalanine, isoleucine, glucose, pyruvate, lactate, 3-hydroxybutyrate)"),
        ("instrument_threshold", "P<5e-8"),
        ("clumping", "greedy; r2<=0.001 within 10 Mb; 1000 Genomes phase 3 EUR (503 samples)"),
        ("clumping_sensitivity", "r2<=0.01, same criteria"),
        ("outcomes", "FinnGen R12: 11 infection endpoints"),
        ("alignment", "ok / flip / palindromic-ambiguous / mismatch; ambiguous palindromes with |EAF-0.5|<0.05 dropped"),
        ("estimators", "Wald ratio (single variant); fixed-effect IVW (multi-variant)"),
        ("diagnostics", "k_eff=(SE_clumped/SE_unclumped)^2; f=(SE_clumped/SE_unclumped)/sqrt(k)"),
        ("criteria", "P1 concentration-f; P2 single vs multi locus; P3 k-k_eff/k; P4 blocks-f; P5 single-block f"),
        ("multiplicity", "Benjamini-Hochberg across 110 pairs within each instrument definition"),
        ("permutation", "20,000 permutations shuffling concentration between traits"),
    ]
    sheet(wb, "S5_methods", ["item", "detail"], [list(x) for x in facts], widths=[22, 110])

    order = ["S1_pairs", "S2_criteria", "S3_r2_sensitivity", "S4_traits",
             "S5_methods", "S6_opengwas"]
    wb._sheets = [wb[t] for t in order if t in wb.sheetnames]
    wb.save(XLSX)
    print("[supp] -> %s S1 %d   / S3 %d   / S4 %d   "
          % (XLSX, len(pairs), len(s3), len(s4)))


if __name__ == "__main__":
    main()
