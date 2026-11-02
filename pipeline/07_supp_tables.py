# -*- coding: utf-8 -*-
"""Build the supplementary tables workbook from the frozen outputs."""
import csv
import sys
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

import openpyxl
from openpyxl.styles import Font, Alignment

ROOT = Path(r"D:/Horzontal_project/20B3")
OUT = ROOT / "work/output"
TOOLS = ROOT / "work/tools"
DEST = ROOT / "work/attachments/Supplementary_Tables_20B3.xlsx"
DEST.parent.mkdir(parents=True, exist_ok=True)


def read_tsv(p):
    """    tab    tab          R data.table::fwrite       """
    with open(p, encoding="utf-8-sig", newline="") as fh:
        head = fh.readline()
    delim = "\t" if "\t" in head else ","
    with open(p, encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter=delim))


def sheet(wb, name, header, rows, first=False, widths=None):
    ws = wb.active if first else wb.create_sheet()
    ws.title = name
    ws.append(header)
    for c in ws[1]:
        c.font = Font(bold=True, size=10)
        c.alignment = Alignment(vertical="center", wrap_text=True)
    for r in rows:
        ws.append(list(r))
    if widths:
        for i, w in enumerate(widths, start=1):
            ws.column_dimensions[openpyxl.utils.get_column_letter(i)].width = w
    ws.freeze_panes = "A2"
    return ws


wb = openpyxl.Workbook()

comp = read_tsv(TOOLS / "b3_stage0_competitors_2026-10-11.tsv")
_seen = set()
s1_rows = []
for i, r in enumerate(comp, 1):
    ab = (r.get("abstract_head") or "").replace("\t", " ").strip()
    pmid = (r.get("pmid") or "").strip()
    title = (r.get("title") or "").strip()
    key = title.lower()
    reason = ""
    if not pmid:
        reason = "Excluded: no PubMed identifier (preprint or unindexed)"
    elif key in _seen:
        reason = "Excluded: duplicate title"
    elif "skin microbiota" in key or "skin microbiome" in key:
        reason = "Excluded: body site outside gut microbiome scope"
    else:
        _seen.add(key)
    s1_rows.append([
        i, pmid, r.get("first_pub", ""), r.get("journal", ""),
        title, (r.get("doi") or ""),
        ("Included" if not reason else reason), ab[:280],
    ])
sheet(wb, "S1_competitors",
      ["#", "PMID", "First publication", "Journal", "Title", "DOI",
       "Inclusion in primary count (33 of 38 raw records)",
       "Abstract (truncated to 280 characters)"],
      s1_rows, first=True, widths=[5, 12, 16, 30, 70, 30, 34, 90])

cs = read_tsv(OUT / "clump_summary.csv")
s2_rows = []
for r in cs:
    s2_rows.append([
        float(r["r2_max"]), int(r["n_taxa"]), int(r["n_snp_before"]),
        int(r["n_snp_after"]), float(r["inflation_factor"]),
    ])
sheet(wb, "S2_clump_sensitivity",
      ["r2 threshold", "Taxa", "SNPs before clumping", "SNPs after clumping",
       "Compression factor"],
      s2_rows, widths=[14, 8, 22, 22, 20])

col = read_tsv(OUT / "coloc_focal.csv")
s3_rows = []
for r in col:
    s3_rows.append([
        r["locus"], r["method"], r["p12"], r["endpoint"], r["endpoint_label"],
        int(r["n_snp_shared"]),
        (round(float(r["PP_H0"]), 4) if r["PP_H0"] else "NA"),
        (round(float(r["PP_H1"]), 4) if r["PP_H1"] else "NA"),
        (round(float(r["PP_H2"]), 4) if r["PP_H2"] else "NA"),
        (round(float(r["PP_H3"]), 4) if r["PP_H3"] else "NA"),
        (round(float(r["PP_H4"]), 4) if r["PP_H4"] else "NA"),
        (round(float(r["ratio_H4_H3"]), 3) if r["ratio_H4_H3"] else "NA"),
        r["verdict"],
    ])
sheet(wb, "S3_colocalization",
      ["Locus", "Method", "p12 prior", "Endpoint key", "Endpoint", "Shared SNPs",
       "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4", "PP.H4 / PP.H3", "Verdict"],
      s3_rows, widths=[8, 13, 10, 14, 46, 12, 9, 9, 9, 9, 9, 14, 20])

stg = read_tsv(OUT / "sens_steiger_retained.csv")
sst = read_tsv(OUT / "sens_steiger_summary.csv")
s4_rows = [[k["metric"], k["value"]] for k in sst]
s4_rows.append(["", ""])
s4_rows.append(["Retained pair", "Value"])
for r in stg:
    s4_rows.append([
        "%s x %s" % (r["name"], r["endpoint"]),
        "R2_exp=%.5f; R2_out_liab=%.5f; one-sided P=%s; direction=%s"
        % (float(r["R2_exp"]), float(r["R2_out_liab"]),
           r["p_steiger_one_sided"], r["direction"])])
s4_rows.append(["", ""])
s4_rows.append(["Effect-estimate displacement", "Value"])
si = read_tsv(OUT / "sens_inflation.csv")
import statistics as _st
_all = [abs(float(r["b_raw"]) - float(r["b_clumped"])) for r in si]
_b100 = [abs(float(r["b_raw"]) - float(r["b_clumped"]))
         for r in si if r["n_snp_bin"] == ">=100"]
s4_rows.append(["Median |change in log-OR|, all matched pairs", "%.4f" % _st.median(_all)])
s4_rows.append(["Maximum |change in log-OR|, all matched pairs", "%.2f" % max(_all)])
s4_rows.append(["Median |change in log-OR|, >=100 instruments", "%.3f" % _st.median(_b100)])
s4_rows.append(["Maximum |change in log-OR|, >=100 instruments", "%.2f" % max(_b100)])
sheet(wb, "S4_direction_and_shift_summary",
      ["Quantity", "Value"], s4_rows, widths=[62, 62])

rb = read_tsv(OUT / "robust_preclump.csv")
s5_rows = []
for r in sorted(rb, key=lambda x: float(x["p_ivw"])):
    s5_rows.append([
        r["endpoint"], r["name"], r["level"], int(r["n_snp"]),
        ("lactase axis" if r["in_lct"] == "1" else "non-lactase"),
        round(float(r["p_ivw"]), 6), round(float(r["p_egger"]), 6),
        round(float(r["p_intercept"]), 6), round(float(r["p_wmed"]), 6),
        round(float(r["I2"]), 2)])
sheet(wb, "S5_robust_preclump",
      ["Endpoint", "Taxon", "Level", "Instruments", "Axis", "IVW P", "MR-Egger P",
       "Egger intercept P", "Weighted median P", "I-squared (%)"],
      s5_rows, widths=[12, 30, 10, 12, 14, 13, 13, 18, 18, 14])

al = read_tsv(OUT / "align_audit_unclumped.csv")
cl = read_tsv(OUT / "align_audit_clumped.csv")
s6_rows = []
for r in al:
    s6_rows.append(["unclumped", r["endpoint"], int(r["n_rows"]), int(r["n_ok"]),
                    int(r["n_flip"]), int(r["n_palin_ambig"]), int(r["n_palin_kept"]),
                    int(r["n_mismatch"])])
for r in cl:
    s6_rows.append(["clumped", r["endpoint"], int(r["n_rows"]), "",
                    "", int(r["n_palin_ambig"]), "", int(r["n_mismatch"])])
sheet(wb, "S6_align_audit",
      ["Definition", "Endpoint", "Rows", "Aligned same strand", "Aligned flipped",
       "Dropped palindromic near MAF 0.5", "Palindromic kept", "Dropped mismatch"],
      s6_rows, widths=[12, 12, 9, 20, 17, 28, 18, 17])

_lc = read_tsv(OUT / "loose_clump_summary.csv")
s7_rows = [[r["r2_max"], r["n_taxa_with_instruments"], r["n_slots"], r["n_tests"],
            r["bh_survivors"], r["min_p"], r["max_n_snp"], r["survivor_names"]]
           for r in _lc]
sheet(wb, "S7_loose_clump",
      ["r2 threshold", "Taxa with instruments", "Instrument slots", "Tests",
       "BH survivors", "Minimum P", "Max variants per taxon", "Surviving taxa"],
      s7_rows, widths=[12, 20, 16, 9, 13, 13, 20, 70])

_lp = read_tsv(OUT / "inflation_law.csv")
s8_rows = []
for r in _lp:
    s8_rows.append([
        r["endpoint"], r["name"], "", r["axis"], int(float(r["k"])),
        round(float(r["se_ratio"]), 4), round(float(r["pred"]), 3),
        round(float(r["se_ratio"]) / float(r["pred"]), 3),
        round(float(r["k_eff"]), 1),
        round(float(r["k_eff_ratio"]), 3),
    ])
sheet(wb, "S8_law_pairs_mibiogen",
      ["Endpoint", "Taxon", "Level", "Axis", "k (instruments)", "SE ratio",
       "Perfect-redundancy bound sqrt(k)", "Observed / bound",
       "Effective count k_eff", "k_eff / k"],
      s8_rows, widths=[12, 30, 10, 24, 14, 10, 28, 16, 18, 10])

_ev = read_tsv(OUT / "law_external_validation.csv")
s9_rows = []
for r in _ev:
    s9_rows.append([
        r["endpoint"], r["name"], r["level"],
        ("lactase axis" if r["in_lct"] == "1" else "non-lactase"),
        int(float(r["k_finnish_matched"])), int(float(r["k_clumped"])),
        round(float(r["b_unclumped"]), 4), round(float(r["se_unclumped"]), 5),
        round(float(r["b_clumped"]), 4), round(float(r["se_clumped"]), 5),
        round(float(r["se_ratio_observed"]), 3), round(float(r["pred_sqrt_k"]), 3),
        round(float(r["k_eff"]), 1), round(float(r["k_eff_over_k"]), 3),
    ])
sheet(wb, "S9_law_pairs_external",
      ["Endpoint", "Taxon", "Level", "Axis", "k matched", "k clumped",
       "b unclumped", "SE unclumped", "b clumped", "SE clumped",
       "SE ratio", "Bound sqrt(k)", "k_eff", "k_eff / k"],
      s9_rows, widths=[12, 26, 9, 15, 11, 11, 12, 14, 12, 14, 10, 14, 10, 10])

_se = read_tsv(OUT / "se_abs_by_bin.csv")
s10_rows = [[r["n_snp_bin"], r["n_pairs"], r["median_se_unclumped"],
             r["median_se_clumped"], r["min_se_unclumped"], r["n_p_underflow"]]
            for r in _se]
sheet(wb, "S10_abs_se",
      ["Instrument bin", "Pairs", "Median SE unclumped", "Median SE clumped",
       "Minimum SE unclumped", "Tests with P underflow"],
      s10_rows, widths=[15, 9, 22, 22, 22, 22])

_cal_rows = []
if (OUT / "keff_calibration.csv").exists():
    for r in read_tsv(OUT / "keff_calibration.csv"):
        _cal_rows.append([float(r["rho"]), int(r["k"]), int(r["n_sim"]),
                          float(r["IF_mean"]), float(r["IF_theory_sqrt_1_plus_k1rho"]),
                          float(r["keff_over_k_mean"])])
    sheet(wb, "S11_keff_calibration",
          ["Instrument correlation rho", "k", "Simulations",
           "Inflation factor (simulated)", "Theoretical sqrt(1+(k-1)rho)",
           "k_eff / k (simulated)"],
          _cal_rows, widths=[24, 8, 13, 26, 30, 22])

_mc = read_tsv(OUT / "mr_all_endpoints_clumped.csv")
_keep = {(r["taxon"], r["endpoint"]) for r in _mc if r.get("sig_bh") == "TRUE"}
_lt = read_tsv(OUT / "loose_clump_tests.csv")
s12_rows = []
for _th in ("0.001", "0.5", "0.8"):
    if _th == "0.001":
        _l = sum(1 for r in _mc if r.get("sig_bh") == "TRUE" and r["in_lct"] == "1")
        _n = sum(1 for r in _mc if r.get("sig_bh") == "TRUE" and r["in_lct"] != "1")
        s12_rows.append([_th, _l + _n, _l, _n, 0, "main analysis"])
        continue
    _sub = [r for r in _lt if r["r2_max"] == _th]
    _ps = sorted(float(r["p"]) for r in _sub)
    _m = len(_ps); _bh = 0
    for _i, _p in enumerate(_ps, 1):
        if _p <= 0.05 * _i / _m:
            _bh = _i
    _thr = _ps[_bh - 1] if _bh else -1.0
    _surv = [r for r in _sub if float(r["p"]) <= _thr]
    _l = sum(1 for r in _surv if r["in_lct"] == "1")
    _n = len(_surv) - _l
    _new = sum(1 for r in _surv if (r["taxon"], r["endpoint"]) not in _keep)
    s12_rows.append([_th, len(_surv), _l, _n, _new, "naive IVW at relaxed threshold"])
sheet(wb, "S12_survivors_by_axis",
      ["Clumping threshold", "BH survivors", "Lactase-axis projections",
       "Non-lactase taxa", "Not significant in the main analysis", "Estimator"],
      s12_rows, widths=[18, 14, 24, 18, 30, 34])

s13_rows = [
    ["Instrument-availability funnel (nominal taxa, taxa with any significant variant, taxa after "
     "removing hierarchical projections)", "Pre-submission", "Author holds"],
    ["Clumped and unclumped estimates side by side, with instrument counts and distinct-variant count",
     "Pre-submission", "Author holds"],
    ["k_eff and f (information-equivalence and attained fraction)", "Reader-side",
     "Two standard errors (clumped, unclumped)"],
    ["Direction-of-effect test", "Reader-side",
     "Per-variant exposure effect sizes and allele frequencies"],
    ["Colocalisation at a retained locus", "Reader-side",
     "Regional summary statistics on both sides"],
    ["Cross-cohort replication of the inflation contrast", "Reader-side",
     "An independent exposure cohort matched by taxon name"],
]
sheet(wb, "S13_reporting_checks",
      ["Check", "When", "Data required"], s13_rows, widths=[70, 16, 46])

_lc_sum = OUT / "loose_coloc_summary.csv"
_lc_raw = OUT / "loose_coloc.csv"
s14_rows = []
_meta = {}
if _lc_sum.exists():
    for r in read_tsv(_lc_sum):
        _meta[r.get("locus", "")] = r
if _lc_raw.exists():
    for r in read_tsv(_lc_raw):
        m = _meta.get(r.get("locus", ""), {})
        s14_rows.append([
            r.get("locus", ""), r.get("endpoint", ""), r.get("taxon", ""),
            m.get("n_exp_region", ""), m.get("n_out_region", ""),
            r.get("n_snp_shared", ""), r.get("method", ""), r.get("p12", ""),
            r.get("PP_H4", ""), r.get("PP_H3", ""), r.get("verdict", ""), ""])
for loc, m in _meta.items():
    if not any(x[0] == loc for x in s14_rows):
        _st = m.get("status", "")
        if _st == "exposure_region_unavailable":
            _st = ("not assessed (exposure table is ordered by P value and has no "
                   "tabix index; a regional subset needs a complete scan)")
        s14_rows.append([
            loc, m.get("endpoint", ""), m.get("taxon", ""),
            m.get("n_exp_region", ""), m.get("n_out_region", ""),
            m.get("n_snp_shared", ""), "not run", "", "", "", "",
            "%s — %s" % (_st, m.get("note", ""))])
sheet(wb, "S14_loose_coloc",
      ["Locus", "Endpoint", "Taxon", "Exposure region n", "Outcome region n",
       "Shared SNPs", "Method", "p12", "PP.H4", "PP.H3", "Verdict", "Status / note"],
      s14_rows, widths=[8, 14, 44, 14, 14, 12, 12, 8, 10, 10, 16, 60])

wb.save(DEST)
print("[07]     -> %s" % DEST)
print("     S1      %d   S2    %d   S3     %d   S4 %d   S5 %d   S6 %d  "
      % (len(s1_rows), len(s2_rows), len(s3_rows), len(s4_rows),
         len(s5_rows), len(s6_rows)))
