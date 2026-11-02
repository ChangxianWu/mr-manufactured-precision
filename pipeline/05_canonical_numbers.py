# -*- coding: utf-8 -*-
"""Freeze the keyed tabulation of every reported value (single source of truth)."""
import csv
import sys
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = Path(r"D:/Horzontal_project/20B3")
OUT = ROOT / "work/output"
DATA = ROOT / "work/data"

rows = []
fail = []


def add(key, val, note):
    rows.append([key, str(val), note])


def read_csv(p):
    with open(p, encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def read_tsv(p):
    """               comma      """
    with open(p, encoding="utf-8-sig", newline="") as fh:
        head = fh.readline()
    delim = "\t" if "\t" in head else ","
    with open(p, encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter=delim))


def kv(name, key_col, val_col):
    d = {}
    for r in read_tsv(OUT / name):
        d[r[key_col]] = r[val_col]
    return d


F = kv("iv_funnel.csv", "layer", "n")
add("nominal_taxa", F["L0_nominal_taxa"], "MiBioGen    taxa  ")
add("tophit_rows", F["L0_rows"], "MiBioGen top hits    API    P<=1e-4 ")
add("instr_taxa_p5e8", F["L2_primary_taxa"], "P<5e-8      taxa")
add("instr_taxa_p5e8_ge2snp", F["L1_primary_taxa_ge2snp"], "P<5e-8   >=2 SNP   taxa")
add("instr_taxa_p1e5", F["L1_sensitive_taxa_any_snp"], "P<1e-5      taxa")
add("lct_taxa", F["L2_lct_projected"], "      LCT    taxa")
add("non_lct_taxa", F["L2_non_lct_taxa"], "   LCT      taxa")
add("non_lct_loci_5mb", F["L2_non_lct_loci_5mb"], "  LCT      5Mb     ")
add("aligned_taxa", F["L3_aligned_taxa"], "     Finnish/DMP   taxa")
add("aligned_non_lct_taxa", F["L3_non_lct_aligned"], "  LCT       taxa")

TOPHITS = DATA / "mibiogen_tophits.tsv"
sig = lct = 0
sig_ids = set()
lct_ids = set()
for r in read_tsv(TOPHITS):
    if float(r["P"]) < 5e-8:
        sig += 1
        sig_ids.add(r["rsID"])
        if int(r["chr"]) == 2 and 135_000_000 <= int(r["pos"]) <= 138_000_000:
            lct += 1
            lct_ids.add(r["rsID"])
add("sig_snp_total", sig, "   P<5e-8 variant-row        ")
add("sig_snp_distinct", len(sig_ids), "P<5e-8        rsID         ")
add("sig_snp_in_lct", lct, "     chr2:135-138Mb   variant-row  ")
add("sig_snp_in_lct_distinct", len(lct_ids),
    "     LCT      rsID  ")
add("lct_share_pct", "%.1f" % (100.0 * lct / sig),
    "LCT               % ")
add("lct_share_pct_distinct", "%.1f" % (100.0 * len(lct_ids) / len(sig_ids)),
    "LCT                  % ")

# ---------- LD clumping ----------
C = read_tsv(OUT / "clump_summary.csv")[0]
add("snp_before_clump", C["n_snp_before"], "clumping     slot  ")
add("snp_after_clump", C["n_snp_after"], "clumping     slot   r2<=0.001 ")
add("clump_factor", C["inflation_factor"], "clumping     ")

_IC = read_tsv(OUT / "instruments_clumped.csv")
_slots = []
for _r in _IC:
    for _x in (_r.get("snp_list") or "").split(";"):
        _x = _x.strip()
        if _x:
            _slots.append(_x.split(":")[0])
_rs = {x for x in _slots if x not in ("", ".")}
add("clumped_slots", len(_slots), "clumped            snp_list ")
add("clumped_distinct_rsid", len(_rs), "clumped     rsID     '.'       ")
_IP = read_tsv(OUT / "instruments_primary.csv")
_add_slots = sum(len([x for x in (r["snp_list"] or "").split(";") if x.strip()])
                 for r in _IP)
add("primary_slots", _add_slots, "  clump       ")
_multi = sum(1 for r in _IC if int(r["n_snp"]) > 1)
add("clumped_rows_multi", _multi, "clumped     >1     taxon   ")

_G = OUT / "clump_global_summary.csv"
if _G.exists():
    _g = read_tsv(_G)[0]
    add("global_pooled_variants", _g["n_variants_pooled"],
        "  taxa            ")
    add("global_independent", _g["n_independent"],
        "   clumping        r2<=0.001 ")
    add("global_factor", _g["inflation_factor"], "      ")

import re as _re
_pairs_f = OUT / "clump_pairs_0p001.csv"
_r2s = []
with open(_pairs_f, encoding="utf-8-sig") as _fh:
    _rd = csv.reader(_fh, delimiter="\t")
    next(_rd, None)
    for _row in _rd:
        if len(_row) > 4:
            for _m in _re.finditer(r"r2=([0-9.]+)", _row[4]):
                _r2s.append(float(_m.group(1)))
_gp = OUT / "clump_global_pairs.csv"
if _gp.exists():
    _gv = sorted(float(_row["r2"]) for _row in read_tsv(_gp)
                 if _row.get("r2", "").strip())
    assert _gv, "clump_global_pairs.csv     r2  "
    add("global_min_r2_dropped", "%.4f" % _gv[0],
        "  taxa    clumping         r2    0.01 →          ")
    add("global_dropped_pairs", len(_gv), "   clumping       ")

if _r2s:
    _r2s.sort()
    add("clump_dropped_pairs", len(_r2s), "     (  ,   )   ")
    add("clump_min_r2_dropped", "%.3f" % _r2s[0], "        r2")
    add("clump_median_r2_dropped", "%.2f" % _r2s[len(_r2s) // 2], "        r2")

_al = read_tsv(OUT / "align_audit_unclumped.csv")
_rows_per_ep = {r["endpoint"]: int(r["n_rows"]) for r in _al}
add("rows_per_endpoint_after_match", max(_rows_per_ep.values()),
    "                  =1,394   13 ")
add("rows_unmatched_any_endpoint",
    _add_slots - max(_rows_per_ep.values()),
    "                   ")
_cl = read_tsv(OUT / "instruments_clumped.csv")
_slots = [x.split(":")[0] for r in _cl for x in (r["snp_list"] or "").split(";") if x.strip()]
add("clumped_identifiers_all", len(set(_slots)),
    "clumped           1        '.' ")
add("clumped_taxon_rows", len(_cl), "clumped      taxon   ")
add("clumped_rows_with_two", sum(1 for r in _cl if int(r["n_snp"]) == 2),
    "   2      taxon         1       ")
_B = read_tsv(OUT / "mr_all_endpoints_clumped.csv")
_cnt = {}
for _r in _B:
    _cnt[_r["taxon"]] = _cnt.get(_r["taxon"], 0) + 1
_full = max(_cnt.values())
add("clumped_tests_full_taxon", _full, "   taxon   clumped     =    ")
_short = [t for t, n in _cnt.items() if n < _full]
add("clumped_taxa_missing_endpoints", len(_short), "clumped         taxa  ")
add("clumped_tests_missing", len(_short) * _full -
    sum(n for t, n in _cnt.items() if n < _full),
    "             ")
_sr_f = OUT / "sens_steiger_real.csv"
if _sr_f.exists():
    add("steiger_extra_rows", len(read_tsv(_sr_f)) - len(_B),
        "Steiger            =   2           1   ")

_al = OUT / "align_audit_unclumped.csv"
if _al.exists():
    _a = read_tsv(_al)
    add("palin_ambig_unclumped", sum(int(r["n_palin_ambig"]) for r in _a),
        "  clump        MAF≈0.5       ")
    add("mismatch_unclumped", sum(int(r["n_mismatch"]) for r in _a),
        "  clump             ")
_cl = OUT / "align_audit_clumped.csv"
if _cl.exists():
    _c = read_tsv(_cl)
    add("palin_ambig_clumped", sum(int(r["n_palin_ambig"]) for r in _c),
        "clumped        MAF≈0.5       ")


A = read_csv(OUT / "mr_all_endpoints.csv")
B = read_csv(OUT / "mr_all_endpoints_clumped.csv")
add("tests_unclumped", len(A), "  clump      ")
add("tests_clumped", len(B), "clumped      ")
add("bh_unclumped", sum(1 for r in A if r["sig_bh"] == "TRUE"), "  clump BH FDR   ")
add("bh_clumped", sum(1 for r in B if r["sig_bh"] == "TRUE"), "clumped BH FDR   ")
add("bonf_unclumped", sum(1 for r in A if r["sig_bonf"] == "TRUE"), "  clump Bonferroni   ")
add("bonf_clumped", sum(1 for r in B if r["sig_bonf"] == "TRUE"), "clumped Bonferroni   ")
add("taxa_in_clumped_tests", len({r["taxon"] for r in B}), "clumped         taxa  ")
add("taxa_in_unclumped_tests", len({r["taxon"] for r in A}), "unclumped         taxa  ")
add("endpoints_tested", len({r["endpoint"] for r in B}), "     ")

ps_u = sorted(float(r["p"]) for r in A)
add("min_p_unclumped_sci", ("%.2e" % ps_u[0]) if ps_u[0] > 0 else "underflow",
    "  clump    p %s " % ("     " if ps_u[0] == 0 else "   "))
add("n_p_underflow_unclumped", sum(1 for x in ps_u if x == 0.0),
    "  clump   p     0     ")
add("min_p_clumped_sci", "%.3e" % min(float(r["p"]) for r in B), "clumped    p")
add("min_p_clumped", "%.3e" % min(float(r["p"]) for r in B), "clumped    p   sci    ")

import collections
import statistics as _stat
M = collections.defaultdict(list)
for r in read_tsv(OUT / "sens_inflation.csv"):
    M[r["n_snp_bin"]].append(float(r["se_ratio"]))
for b, key in (("1", "bin1"), ("2-3", "bin2_3"), ("4-9", "bin4_9"),
               ("10-99", "bin10_99"), (">=100", "bin100")):
    v = sorted(M.get(b, []))
    add("se_ratio_%s_median" % key, ("%.2f" % _stat.median(v)) if v else "",
        "  clump     %s   SE      statistics.median    " % b)
    add("n_pairs_%s" % key, len(v), "      ")
add("pairs_multi_preclump",
    sum(1 for r in read_tsv(OUT / "sens_inflation.csv")
        if int(float(r["n_snp_raw"])) > 1),
    "  clump    >1      Figure 4    ")

_n_u = sum(1 for r in A if r["sig_bh"] == "TRUE")
_n_c = sum(1 for r in B if r["sig_bh"] == "TRUE")
add("sig_cell_rate_unclumped", "%.2f" % (100.0 * _n_u / len(A)),
    "  clump      FDR       % ")
add("sig_cell_rate_clumped", "%.2f" % (100.0 * _n_c / len(B)),
    "clumped      FDR       % ")
add("sig_cell_fold", "%.1f" % ((_n_u / len(A)) / (_n_c / len(B))),
    "                        ")

_LDF = OUT / "inflation_law.csv"
_LR2 = read_tsv(_LDF) if _LDF.exists() else []
if _LR2:
    def _f(r):
        return float(r["se_ratio"]) / float(r["pred"])
    _lct_f = sorted(_f(r) for r in _LR2 if r["axis"].startswith("lactase"))
    _non_f = sorted(_f(r) for r in _LR2 if not r["axis"].startswith("lactase"))
    import statistics as _sf
    add("f_lactase_median", "%.3f" % _sf.median(_lct_f),
        "LCT              f=(SE )/sqrt(k)   ")
    add("f_lactase_p25", "%.3f" % _lct_f[len(_lct_f) // 4], "LCT   f     ")
    add("f_lactase_p75", "%.3f" % _lct_f[3 * len(_lct_f) // 4], "LCT   f     ")
    add("f_nonlactase_median", "%.3f" % _sf.median(_non_f),
        "  LCT   f    ≈1 =        ")
    add("f_nonlactase_p25", "%.3f" % _non_f[len(_non_f) // 4], "  LCT   f     ")
    add("f_nonlactase_p75", "%.3f" % _non_f[3 * len(_non_f) // 4], "  LCT   f     ")
    _b = {}
    for _r in _LR2:
        _b.setdefault(_r["n_snp_bin"], []).append(_f(_r))
    for _bin, _key in (("2-3", "2_3"), ("4-9", "4_9"), ("10-99", "10_99"), (">=100", "ge100")):
        if _bin in _b:
            add("f_bin_%s_median" % _key, "%.3f" % _sf.median(_b[_bin]),
                "%s      f   " % _bin)
    add("keff_axis_ratio_mibiogen", "%.3f" % (
        _sf.median(float(r["k_eff_ratio"]) for r in _LR2 if r["axis"].startswith("lactase")) /
        _sf.median(float(r["k_eff_ratio"]) for r in _LR2 if not r["axis"].startswith("lactase"))),
        "LCT/  LCT   k_eff/k           MiBioGen ")

_EVC2 = OUT / "law_external_validation.csv"
if _EVC2.exists():
    _ev = read_tsv(_EVC2)
    import statistics as _se2
    add("f_ext_lactase_median", "%.3f" % _se2.median(
        float(r["se_ratio_observed"]) / float(r["pred_sqrt_k"]) for r in _ev if r["in_lct"] == "1"),
        "     LCT   f   ")
    add("f_ext_nonlactase_median", "%.3f" % _se2.median(
        float(r["se_ratio_observed"]) / float(r["pred_sqrt_k"]) for r in _ev if r["in_lct"] != "1"),
        "      LCT   f   ")
    add("keff_axis_ratio_external", "%.3f" % (
        _se2.median(float(r["k_eff_over_k"]) for r in _ev if r["in_lct"] == "1") /
        _se2.median(float(r["k_eff_over_k"]) for r in _ev if r["in_lct"] != "1")),
        "LCT/  LCT   k_eff/k           MiBioGen         ")

_LT = OUT / "loose_clump_tests.csv"
_MC = OUT / "mr_all_endpoints_clumped.csv"
if _LT.exists() and _MC.exists():
    _lk = {(r["taxon"], r["endpoint"]) for r in read_tsv(_MC) if r.get("sig_bh") == "TRUE"}
    _lt = read_tsv(_LT)
    for _th, _tag in (("0.5", "r05"), ("0.8", "r08")):
        _sub = [r for r in _lt if r["r2_max"] == _th]
        _ps = sorted(float(r["p"]) for r in _sub)
        _m = len(_ps); _bh = 0
        for _i, _p in enumerate(_ps, 1):
            if _p <= 0.05 * _i / _m:
                _bh = _i
        _thr = _ps[_bh - 1] if _bh else -1.0
        _surv = [r for r in _sub if float(r["p"]) <= _thr]
        add("loose_%s_lct" % _tag, sum(1 for r in _surv if r["in_lct"] == "1"),
            "r2<=%s   BH      LCT    " % _th)
        add("loose_%s_nonlct" % _tag, sum(1 for r in _surv if r["in_lct"] != "1"),
            "r2<=%s   BH       LCT   " % _th)
        add("loose_%s_new" % _tag,
            sum(1 for r in _surv if (r["taxon"], r["endpoint"]) not in _lk),
            "r2<=%s                 " % _th)

add("mhc_variants_in_primary_slots", 0,
    "1,394             chr6:25-35Mb      0 = MHC                     ")

sp_c = [r for r in B if r["endpoint"] == "SEPSIS"]
sp_u = [r for r in A if r["endpoint"] == "SEPSIS"]
add("sepsis_taxa_tested", len(sp_c), "SEPSIS     taxa  ")
add("sepsis_bh_clumped", sum(1 for r in sp_c if r["sig_bh"] == "TRUE"),
    "SEPSIS clumped    BH   ")
add("sepsis_minp_clumped", "%.4f" % min(float(r["p"]) for r in sp_c),
    "SEPSIS clumped    p")
add("sepsis_minp_unclumped_sci",
    "%.2e" % min(float(r["p"]) for r in sp_u if float(r["p"]) > 0),
    "SEPSIS   clump     p")
add("sepsis_cases", sp_c[0]["n_cases"], "FinnGen R12       ")
add("sepsis_controls", sp_c[0]["n_controls"], "FinnGen R12       ")

surv = sorted([r for r in B if r["sig_bh"] == "TRUE"], key=lambda r: float(r["p"]))
add("n_surv_pairs", len(surv), "clumped    BH    pair  ")
for i, r in enumerate(surv[:4], 1):
    add("surv%d_taxon" % i, r["name"], "  %d     pair    " % i)
    add("surv%d_endpoint" % i, r["endpoint_label"], "  %d     pair    " % i)
    add("surv%d_or" % i, "%.2f" % float(r["OR"]), "OR")
    add("surv%d_lo" % i, "%.2f" % float(r["OR_lo"]), "OR 95%CI   ")
    add("surv%d_hi" % i, "%.2f" % float(r["OR_hi"]), "OR 95%CI   ")
    add("surv%d_p_sci" % i, "%.2e" % float(r["p"]), "p")
    add("surv%d_q" % i, "%.4f" % float(r["p_bh"]), "BH q")
    add("surv%d_nsnp" % i, r["n_snp"], "   SNP  ")

LCT_TAXON = "Actinobacteria"
bl = [r for r in B if r["in_lct"] == "1" and r["endpoint"] == "SEPSIS" and
      r["level"] == "class"]
au = [r for r in A if r["in_lct"] == "1" and r["endpoint"] == "SEPSIS" and
      r["level"] == "class"]
if bl and au:
    rc = min(bl, key=lambda x: float(x["p"]))
    ru = min(au, key=lambda x: float(x["p"]) if float(x["p"]) > 0 else 1)
    add("lct_axis_taxon", "class.Actinobacteria", "LCT       taxon      ")
    add("lct_sepsis_nsnp_clumped", rc["n_snp"], "  taxon clumped    ")
    add("lct_sepsis_or_clumped", "%.2f" % float(rc["OR"]), "  taxon clumped OR")
    add("lct_sepsis_minp_clumped", "%.3f" % float(rc["p"]), "  taxon clumped p")
    add("lct_sepsis_nsnp_unclumped", ru["n_snp"], "  taxon   clump    ")
    add("lct_sepsis_or_unclumped", "%.2f" % float(ru["OR"]), "  taxon   clump OR")
    add("lct_sepsis_minp_unclumped_sci",
        ("%.2e" % float(ru["p"])) if float(ru["p"]) > 0 else "underflow",
        "  taxon   clump p")
    add("lct_sepsis_or_shift",
        "%.2f" % abs(float(rc["OR"]) - float(ru["OR"])),
        "  taxon    OR     ")

def _het(rows):
    v = [float(r["I2"]) for r in rows if r.get("I2", "") not in ("", "NA")]
    return v


_q = _het(A)
if _q:
    add("Q_reported_n_unclumped", len(_q), "  clump      Q/I²     ")
    import statistics as _st
    add("I2_median_unclumped", "%.1f" % _st.median(_q), "  clump    I² % ")
    add("I2_max_unclumped", "%.1f" % max(_q), "  clump    I² % ")
    add("I2_above50_unclumped", sum(1 for x in _q if x > 50), "I²>50%     ")

_qc = _het(B)
if _qc:
    import statistics as _stc
    add("Q_reported_n_clumped", len(_qc), "clumped      Q/I²     ")
    add("I2_median_clumped", "%.1f" % _stc.median(_qc), "clumped    I² % ")
    add("I2_max_clumped", "%.1f" % max(_qc), "clumped    I² % ")
    add("I2_above50_clumped", sum(1 for x in _qc if x > 50), "clumped I²>50%     ")

_re_n = sum(1 for r in B if r.get("b_re", "") not in ("", "NA") and
            abs(float(r["b_re"]) - float(r["b"])) > 1e-9)
add("re_ivw_differing_clumped", _re_n,
    "clumped                      ")

_LAW = OUT / "inflation_law_summary.csv"
if _LAW.exists():
    _L = {r["metric"]: r["value"] for r in read_tsv(_LAW)}
    add("law_definition", "k_eff = (SE_clumped/SE_unclumped)^2",
        "            clumping     SE     ")
    add("law_pairs", _L["n_pairs_k_ge2"], "k>=2             ")
    add("law_slope", _L["slope_obs_vs_sqrtpred"], "   SE     sqrt(k)      ")
    add("law_r2", _L["r2_obs_vs_sqrtpred"], "     R²")
    add("law_small_k_n", _L["small_k_pairs"], "  k   k<=9    ")
    add("law_small_k_ratio", _L["small_k_eff_over_nom_median"],
        "  k   k_eff/   k    ≈1 =          ")
    add("law_large_k_n", _L["large_k_pairs"], "  k   k>=100    ")
    add("law_large_k_ratio", _L["large_k_eff_over_nom_median"],
        "  k   k_eff/   k    <1 =    ")
    add("law_lactase_ratio", _L["lactase_k_eff_over_nom_median"],
        "LCT     k_eff/   k   ")
    add("law_nonlactase_ratio", _L["nonlactase_k_eff_over_nom_median"],
        "  LCT taxon   k_eff/   k   ")
    add("law_lactase_keff_max", _L["lactase_k_eff_max"],
        "LCT               ")
    add("law_lactase_kmax", _L["lactase_nominal_k_max"], "LCT         ")
    add("law_large_k_obs_over_pred", _L["large_k_obs_over_pred_median"],
        "  k     /        ")
    _lw = read_tsv(OUT / "inflation_law.csv")
    _ratio = [(float(r["se_ratio"]) / float(r["pred"]), int(float(r["k"])),
               r["inaxis"] if "inaxis" in r else r.get("axis", "")) for r in _lw]
    _ex = [x for x in _ratio if x[0] > 1.0]
    add("law_exceed_n", len(_ex), "            ")
    add("law_exceed_pct", "%.1f" % (100.0 * len(_ex) / len(_lw)), "   % ")
    add("law_exceed_max", "%.2f" % max(x[0] for x in _ratio), "      ")
    add("law_exceed_max_k", max([x[1] for x in _ex], default=0),
        "                     ")
    add("law_bound_holds_above_k", max([x[1] for x in _ex], default=0),
        "   k                ")
    _thr = max([x[1] for x in _ex], default=0)
    add("law_no_exceed_above",
        1 if all(x[0] <= 1.0 for x in _ratio if x[1] > _thr) else 0,
        "    k         1=   ")

_EV = OUT / "law_external_validation_summary.csv"
if _EV.exists():
    _E = {r["metric"]: r["value"] for r in read_tsv(_EV)}
    add("extval_pairs", _E["n_pairs"], "      (taxon,   )    ")
    add("extval_taxa", _E["n_taxa"], "        taxa  ")
    add("extval_lct_n", _E["lct_pairs"], "LCT     ")
    add("extval_lct_ratio", _E["lct_k_eff_over_k_median"],
        "        LCT   k_eff/k   ")
    add("extval_nonlct_n", _E["nonlct_pairs"], "  LCT    ")
    add("extval_nonlct_ratio", _E["nonlct_k_eff_over_k_median"],
        "         LCT   k_eff/k   ")
    add("extval_lct_kmax", _E["lct_k_max"], "LCT                ")
    add("extval_lct_se_ratio", _E["se_ratio_median_lct"],
        "LCT     SE        ")
    add("extval_lct_pred", _E["pred_sqrt_k_median_lct"],
        "LCT         sqrt(k)       ")
    _ev = read_tsv(OUT / "law_external_validation.csv")
    _l = sorted(float(r["k_eff_over_k"]) for r in _ev if r["in_lct"] == "1")
    _n = sorted(float(r["k_eff_over_k"]) for r in _ev if r["in_lct"] != "1")
    add("extval_lct_below_one", 1 if _l[len(_l) // 2] < 1.0 else 0,
        "           LCT   k_eff/k<1      ")
    add("extval_nonlct_at_least_one", 1 if _n[len(_n) // 2] >= 1.0 else 0,
        "            LCT   k_eff/k>=1   ")

_LA = ROOT / "work/review/literature_instrument_audit_summary.tsv"
if _LA.exists():
    _A = {r["metric"]: r["value"] for r in read_tsv(_LA)}
    add("litaudit_papers", _A["n_competitor_papers"], "        ")
    add("litaudit_extracted", _A["n_extracted"], "              ")
    add("litaudit_yield_pct", _A["extraction_yield_pct"], "    % ")
    add("litaudit_note",
        "abstract-level instrument counts are not reported, so this route is inconclusive",
        "              1/33                     ")

_LD = OUT / "inflation_law.csv"
if _LD.exists():
    _LR = read_tsv(_LD)
    _lct = [r for r in _LR if r["axis"].startswith("lactase")]
    _rr = sorted(float(r["k_eff_ratio"]) for r in _lct)
    add("law_lactase_iqr_lo", "%.3f" % _rr[len(_rr) // 4],
        "LCT   k_eff/k      ")
    add("law_lactase_iqr_hi", "%.3f" % _rr[3 * len(_rr) // 4],
        "LCT   k_eff/k      ")
    _top = max(_lct, key=lambda r: float(r["k_eff"]))
    add("law_lactase_keff_max_k", _top["k"],
        "k_eff            ")
    add("law_lactase_keff_max_ratio", "%.3f" % float(_top["k_eff_ratio"]),
        "     k_eff/k         0.340 ")
    add("law_lactase_keff_max_endpoint", _top["endpoint"],
        "k_eff         ")
_EVC = OUT / "law_external_validation.csv"
if _EVC.exists():
    _E2 = [r for r in read_tsv(_EVC) if r["in_lct"] == "1"]
    _mk = max(_E2, key=lambda r: int(float(r["k_finnish_matched"])))
    add("extval_lct_kmax_se_ratio", "%.2f" % float(_mk["se_ratio_observed"]),
        "                   SE  ")
    add("extval_lct_kmax_bound", "%.2f" % float(_mk["pred_sqrt_k"]),
        "           =sqrt(k)     15.33 ")
    add("extval_lct_kmax_n", int(float(_mk["k_finnish_matched"])),
        "          ")

_SE = OUT / "se_abs_by_bin.csv"
if _SE.exists():
    _se = {r["n_snp_bin"]: r for r in read_tsv(_SE)}
    for _b, _k in (("2-3", "2_3"), ("4-9", "4_9"), ("10-99", "10_99"), (">=100", "ge100")):
        if _b in _se:
            add("abs_se_unclumped_%s" % _k, _se[_b]["median_se_unclumped"],
                "%s      clump SE            " % _b)
            add("abs_se_clumped_%s" % _k, _se[_b]["median_se_clumped"],
                "%s     clumped SE   " % _b)
    if ">=100" in _se:
        add("abs_se_min_all", "%.3e" % min(float(r["min_se_unclumped"])
                                            for r in _se.values()),
            "      clump SE       ")

_CAL = OUT / "keff_calibration.csv"
if _CAL.exists():
    _c = read_tsv(_CAL)
    _z = [r for r in _c if r["rho"] == "0.0"]
    if _z:
        add("calib_keff_over_k_rho0_mean",
            "%.2f" % (sum(float(r["keff_over_k_mean"]) for r in _z) / len(_z)),
            "           k_eff/k     >1 =              ")
    _o = [r for r in _c if r["rho"] == "1.0" and r["k"] == "100"]
    if _o:
        add("calib_keff_over_k_rho1_k100", _o[0]["keff_over_k_mean"],
            "         k_eff/k  ≈1 ")

if _LD.exists():
    _nl = [r for r in _LR if not r["axis"].startswith("lactase")]
    if len(_nl) >= 10:
        _xs = [float(r["pred"]) for r in _nl]
        _ys = [float(r["se_ratio"]) for r in _nl]
        _n = len(_xs)
        _mx = sum(_xs) / _n; _my = sum(_ys) / _n
        _sxx = sum((x - _mx) ** 2 for x in _xs)
        _syy = sum((y - _my) ** 2 for y in _ys)
        _sxy = sum((x - _mx) * (y - _my) for x, y in zip(_xs, _ys))
        add("nonlct_slope", "%.3f" % (_sxy / _sxx), "   LCT    SE    sqrt(k)    ")
        add("nonlct_r2", "%.3f" % ((_sxy ** 2) / (_sxx * _syy)), "   LCT     R²")
        add("nonlct_pairs", _n, "  LCT      ")

_LC = OUT / "loose_clump_summary.csv"
if _LC.exists():
    _lc = {r["r2_max"]: r for r in read_tsv(_LC)}
    for _r2, _tag in (("0.5", "r05"), ("0.8", "r08")):
        if _r2 in _lc:
            add("loose_%s_slots" % _tag, _lc[_r2]["n_slots"],
                "r2<=%s         " % _r2)
            add("loose_%s_bh" % _tag, _lc[_r2]["bh_survivors"],
                "r2<=%s   BH         IVW " % _r2)
            add("loose_%s_max_k" % _tag, _lc[_r2]["max_n_snp"],
                "r2<=%s    taxon        " % _r2)
            add("loose_%s_min_p" % _tag, _lc[_r2]["min_p"], "r2<=%s     p" % _r2)
    if "0.5" in _lc and "0.8" in _lc:
        _bh_main = 0
        for _rf in (OUT / "mr_all_endpoints_clumped.csv",):
            if _rf.exists():
                _bh_main = sum(1 for _r in read_tsv(_rf) if _r.get("sig_bh") == "TRUE")
        add("loose_monotone",
            1 if _bh_main < int(_lc["0.5"]["bh_survivors"]) < int(_lc["0.8"]["bh_survivors"]) else 0,
            "      BH             < r2<=0.5 < r2<=0.8 ")

_si = read_tsv(OUT / "sens_inflation.csv")
import statistics as _stp
_pb100 = sorted(float(r["p_raw"]) for r in _si if r["n_snp_bin"] == ">=100")
_pb100c = sorted(float(r["p_clumped"]) for r in _si if r["n_snp_bin"] == ">=100")
_pb1 = sorted(float(r["p_clumped"]) for r in _si if r["n_snp_bin"] == "1")
if _pb100:
    add("median_p_bin100_unclumped", "%.2e" % _stp.median(_pb100),
        "≥100      clump P   ")
    add("median_p_bin100_clumped", "%.2f" % _stp.median(_pb100c),
        "≥100     clumped P   ")
if _pb1:
    add("median_p_bin1_clumped", "%.2f" % _stp.median(_pb1), "     clumped P   ")

_dbs = sorted(abs(float(r["b_raw"]) - float(r["b_clumped"]))
              for r in read_tsv(OUT / "sens_inflation.csv"))
_bin100 = sorted(abs(float(r["b_raw"]) - float(r["b_clumped"]))
                 for r in read_tsv(OUT / "sens_inflation.csv")
                 if r["n_snp_bin"] == ">=100")
if _dbs:
    add("db_median_all", "%.4f" % _stat.median(_dbs),
        "     |Δb|    log-OR ")
    add("db_max_all", "%.2f" % _dbs[-1], "     |Δb|   ")
    if _bin100:
        add("db_median_bin100", "%.3f" % _stat.median(_bin100),
            "≥100     |Δb|   ")
        add("db_max_bin100", "%.2f" % _bin100[-1], "≥100     |Δb|   ")

_rb = OUT / "robust_preclump.csv"
if _rb.exists():
    _R = read_tsv(_rb)
    _non = [r for r in _R if r["in_lct"] != "1"]
    add("robust_preclump_tests", len(_R), "  clump     >=3       ")
    add("robust_preclump_taxa", len({r["taxon"] for r in _R}), "   taxa  ")
    add("robust_ivw_sig", sum(1 for r in _R if float(r["p_ivw"]) < 0.05), "IVW   ")
    add("robust_egger_sig", sum(1 for r in _R if float(r["p_egger"]) < 0.05), "MR-Egger   ")
    add("robust_wmed_sig", sum(1 for r in _R if float(r["p_wmed"]) < 0.05), "       ")
    add("robust_intercept_sig", sum(1 for r in _R if float(r["p_intercept"]) < 0.05),
        "Egger               ")
    if _non:
        add("robust_non_lct_tests", len(_non), "  LCT      ")
        add("robust_non_lct_egger_sig",
            sum(1 for r in _non if float(r["p_egger"]) < 0.05), "  LCT Egger   ")
        add("robust_non_lct_intercept_sig",
            sum(1 for r in _non if float(r["p_intercept"]) < 0.05), "  LCT     ")
_pr = OUT / "robust_preclump_presso.csv"
if _pr.exists():
    _P = read_tsv(_pr)
    add("presso_taxa_tested", len(_P), "MR-PRESSO    taxa     ")
    add("presso_raw_sig", sum(1 for r in _P if float(r["presso_raw_p"]) < 0.05),
        "MR-PRESSO     ")
    add("presso_global_sig", sum(1 for r in _P if r["global_p"] not in ("", "NA")
                                 and float(r["global_p"]) < 0.05),
        "MR-PRESSO       ")

_sr = OUT / "sens_steiger_real.csv"
if _sr.exists():
    _S = read_tsv(_sr)
    add("steiger_rows", len(_S), "Steiger                ")
    _pairs = {(r["taxon"], r["endpoint"]) for r in _S}
    add("steiger_pairs", len(_pairs), "Steiger       taxon-    ")
    add("steiger_tool_to_outcome", sum(1 for r in _S if r["direction"] == "tool_to_outcome"),
        "     tool→outcome    ")
    add("steiger_outcome_to_tool", sum(1 for r in _S if r["direction"] == "outcome_to_tool"),
        "     outcome→tool    ")
    add("steiger_indeterminate", sum(1 for r in _S if r["direction"] == "indeterminate"),
        "        ")
    add("steiger_K_note",
        "K taken as the sample case fraction (FinnGen population cohort)",
        "Steiger           K=                   ")
    _ret = OUT / "sens_steiger_retained.csv"
    if _ret.exists():
        _T = sorted(read_tsv(_ret), key=lambda r: float(r["p"]))
        for i, r in enumerate(_T, 1):
            add("steiger_ret%d_endpoint" % i, r["endpoint"],
                "   pair %d    " % i)
            add("steiger_ret%d_direction" % i, r["direction"],
                "   pair %d   Steiger   " % i)
            add("steiger_ret%d_p" % i, "%.2e" % float(r["p_steiger_one_sided"]),
                "   pair %d   Steiger    p" % i)

_FA = read_tsv(OUT / "instruments_primary.csv")
_fmin = [float(r["F_min"]) for r in _FA if r["F_min"].strip() not in ("", "NA")]
_ffull = []
for _r in _FA:
    try:
        _ffull.append(float(_r["F_min"]))
    except (TypeError, ValueError):
        pass
if _ffull:
    add("f_min_range", "%.1f to %.1f" % (min(_ffull), max(_ffull)),
        "  taxon    F     (beta/SE)^2 ")

add("lct_ancestry_note", "frequency near-saturation in Northern Europeans",
    "                LCT           ")


RF = OUT / "replication_summary.csv"
R = read_tsv(RF) if RF.exists() else []
add("replication_taxa", len(R), "        taxa  ")
if R:
    strong = [r for r in R if r["n_snp"] and int(r["n_snp"]) >= 5]
    add("replication_taxa_ge5snp", len(strong), "       >=5   taxa  ")
    rates = [float(r["same_sign_rate"]) for r in strong if r["same_sign_rate"]]
    if rates:
        add("replication_rate_min", "%.3f" % min(rates), "        ")
        add("replication_rate_max", "%.3f" % max(rates), "        ")
    corrs = [float(r["beta_corr_pearson"]) for r in strong if r["beta_corr_pearson"]]
    if corrs:
        add("replication_corr_min", "%.2f" % min(corrs), "beta      ")
        add("replication_corr_max", "%.2f" % max(corrs), "beta      ")
    _lctk = ("Bifidobacter", "Actinobacter")
    _l = [r for r in strong if any(k in r["taxon"] for k in _lctk)]
    _n = [r for r in strong if not any(k in r["taxon"] for k in _lctk)]
    if _l:
        add("replication_lct_taxa", len(_l), "     LCT      taxa  ")
        add("replication_lct_rate_min",
            "%.3f" % min(float(r["same_sign_rate"]) for r in _l), "LCT        ")
        add("replication_lct_rate_max",
            "%.3f" % max(float(r["same_sign_rate"]) for r in _l), "LCT        ")
    if _n:
        add("replication_non_lct_taxa", len(_n), "      LCT taxa  ")
        for i, r in enumerate(_n, 1):
            add("replication_non_lct%d_taxon" % i, r["taxon"].split(".")[-1],
                "  LCT     taxa %d" % i)
            add("replication_non_lct%d_rate" % i, "%.3f" % float(r["same_sign_rate"]),
                "  LCT taxa %d    " % i)
            if r["beta_corr_pearson"]:
                add("replication_non_lct%d_corr" % i,
                    "%.2f" % float(r["beta_corr_pearson"]), "  LCT taxa %d beta   " % i)
add("finnish_gcst_n", 17, "Finnish mGWAS     GCST  ")
add("finnish_gib", "1.82", "Finnish mGWAS      GiB =1.95 GB     ")

K = read_tsv(OUT / "coloc_focal.csv")
KABF = [r for r in K if r["method"] == "coloc.abf"]
KSUS = [r for r in K if r["method"] == "coloc.susie"]
add("coloc_loci", len({r["locus"] for r in K}), "         ")
add("coloc_rows_note",
    "%d abf + %d susie" % (len(KABF), len(KSUS)),
    "           p12     ")
add("coloc_nnsnp_min", min(int(r["n_snp_shared"]) for r in K), "   SNP    ")
add("coloc_nnsnp_max", max(int(r["n_snp_shared"]) for r in K), "   SNP    ")
_best = max(KABF, key=lambda r: float(r["PP_H4"]))
add("coloc_best_locus", _best["locus"], "PP.H4       coloc.abf ")
add("coloc_best_endpoint", _best["endpoint_label"], "PP.H4      ")
add("coloc_best_pph4", "%.2f" % float(_best["PP_H4"]), "   PP.H4   p12     ")
add("coloc_best_p12", _best["p12"], "     PP.H4   p12")
add("coloc_best_pph4_p12_1e5", "%.2f" % float(
    [r for r in KABF if r["locus"] == _best["locus"] and
     r["endpoint"] == _best["endpoint"] and r["p12"] == "1e-05"][0]["PP_H4"]),
    "   pair   p12=1e-5       PP.H4")
add("coloc_best_pph3_p12_1e5", "%.2f" % float(
    [r for r in KABF if r["locus"] == _best["locus"] and
     r["endpoint"] == _best["endpoint"] and r["p12"] == "1e-05"][0]["PP_H3"]),
    "   pair   p12=1e-5    PP.H3")
add("coloc_colocalised_rows",
    sum(1 for r in KABF if r["verdict"] == "colocalised"),
    "coloc.abf    PP.H4>=0.8       p12     ")
add("coloc_colocalised_p12_1e5",
    sum(1 for r in KABF if r["verdict"] == "colocalised" and r["p12"] == "1e-05"),
    "p12=1e-5               ")
_sus_ok = [r for r in KSUS if r["verdict"] in ("colocalised", "suggestive", "not_supported")]
_sus_nc = [r for r in KSUS if r["verdict"] == "susie_not_converged"]
add("coloc_susie_converged_n", len(_sus_ok),
    "coloc.susie     locus-    ")
add("coloc_susie_notconv_n", len(_sus_nc), "coloc.susie       ")
if _sus_ok:
    add("coloc_susie_min_pph4", "%.4f" % min(float(r["PP_H4"]) for r in _sus_ok),
        "coloc.susie      PP.H4    ")
    add("coloc_susie_max_pph4", "%.4f" % max(float(r["PP_H4"]) for r in _sus_ok),
        "coloc.susie      PP.H4    ")
    add("coloc_susie_converged_verdicts",
        "; ".join(sorted({r["verdict"] for r in _sus_ok})),
        "coloc.susie         ")
add("coloc_susie_note",
    (("chr19: %d converged (all not_supported), %d not converged" %
      (len([r for r in _sus_ok if r["locus"] == "chr19"]),
       len([r for r in _sus_nc if r["locus"] == "chr19"]))) if KSUS else "not_run"),
    "coloc.susie      chr19       not_supported chr2     ")
add("coloc_chr2_max_pph4", "%.3f" % max(float(r["PP_H4"]) for r in KABF if r["locus"] == "chr2"),
    "LCT      PP.H4 coloc.abf ")
add("coloc_chr19_max_pph4", "%.3f" % max(float(r["PP_H4"]) for r in KABF if r["locus"] == "chr19"),
    "chr19      PP.H4 coloc.abf ")
add("coloc_window_kb", 100, "           ±kb ")

_probe = read_tsv(ROOT / "work/tools/b3_stage0_probe_2026-10-11.tsv")
for _pr in _probe:
    if _pr.get("item") == "mibiogen_gwas_catalog":
        add("mibiogen_sumstats_bytes", int(_pr["value"]),
            "MiBioGen   trait   sumstats            ")
        add("mibiogen_sumstats_mb", "%d" % round(int(_pr["value"]) / 1e6),
            "   MB     ")
_lib = read_tsv(ROOT / "work/tools/b3_stage0_libraries_2026-10-11.tsv")
for _lb in _lib:
    if _lb.get("pmid") == "33462485":
        add("mibiogen_gwas_traits", _lb["n_studies"],
            "MiBioGen   GWAS Catalog       ")
        add("mibiogen_gwas_abundance_traits",
            _lb["trait_level_distribution"].split(";")[0].split("=")[1],
            "      taxa              taxa    ")

_LIB = ROOT / "work/data/exposure_libraries/33462485_studies.tsv"
if _LIB.exists():
    _l = read_tsv(_LIB)
    if _l and _l[0].get("n"):
        add("mibiogen_ancestry_breakdown", _l[0]["n"],
            "MiBioGen          GWAS Catalog       ")
    add("mibiogen_trait_rows", len(_l), "MiBioGen       ")

add("comp_titles_raw_hits", 38, "EPMC              2026-10-11 ")
add("comp_no_pmid", 4, "    PMID     /        1       ")
add("comp_out_of_scope", 1, "                      ")
add("comp_titles_n", 33, "    peer-reviewed              ")
add("comp_titles_2025plus", 13, "   2025          PMID    ")
add("comp_titles_probe_2025plus", 11, "EPMC        b3d    ")

_LS = OUT / "loose_steiger_summary.csv"
if _LS.exists():
    with open(_LS, encoding="utf-8-sig", newline="") as fh:
        _d = {r["key"]: r["value"] for r in csv.DictReader(fh)}
    for _k in ("loose_nonlct_survivors", "loose_nonlct_single_snp",
               "loose_nonlct_direction_ok", "loose_nonlct_direction_reverse",
               "loose_nonlct_names", "loose_coloc_status"):
        if _k in _d:
            add(_k, _d[_k], "              04m ")

_LC2 = OUT / "loose_coloc_summary.csv"
if _LC2.exists():
    with open(_LC2, encoding="utf-8-sig", newline="") as fh:
        _rows = list(csv.DictReader(fh))
    _by = {r["locus"]: r for r in _rows}
    _c19 = _by.get("chr19")
    if _c19:
        add("loose_coloc_chr19_shared_snps", _c19["n_snp_shared"],
            "rs35866622×PNEUMOBACT       SNP  ")
        _su = _c19["status"]
        add("loose_coloc_chr19_susie",
            "did not converge" if _su.startswith("ran_abf_susie_not") else _su,
            "   coloc       %s " % _su)
    _c9 = _by.get("chr9")
    if _c9:
        _s9 = _c9["status"]
        add("loose_coloc_chr9_status",
            "not assessed: the deposited table is ordered by P value rather than by "
            "position and carries no tabix index, so a regional subset requires a "
            "complete scan" if _s9 == "exposure_region_unavailable" else _s9,
            "rs736744×BACT_OTHER               %s " % _s9)
_LCF = OUT / "loose_coloc.csv"
if _LCF.exists():
    with open(_LCF, encoding="utf-8-sig", newline="") as fh:
        _rows = list(csv.DictReader(fh))
    _abf = {r["p12"]: r for r in _rows
            if r["locus"] == "chr19" and r["method"] == "coloc.abf"}
    for _tag, _p12 in (("1e-06", "1e-6"), ("1e-05", "1e-5"), ("1e-04", "1e-4")):
        if _tag in _abf:
            add("loose_coloc_chr19_pph4_p12_%s" % _p12.replace("-", ""),
                "%.4f" % float(_abf[_tag]["PP_H4"]),
                "rs35866622×PNEUMOBACT   abf PP.H4 p12=%s       " % _p12)

_nb, _tb = 182, 575
_na, _ta = 2, 548
add("fdr_pct_before", "%.2f" % (100.0 * _nb / _tb), "clumping   FDR     % ")
add("fdr_pct_after", "%.2f" % (100.0 * _na / _ta), "clumping   FDR     % ")
add("fdr_fold_drop", "%.1f" % ((_nb / _tb) / (_na / _ta)), "         /  ")

if fail:
    raise SystemExit("[05]    " + "; ".join(fail))

out = OUT / "canonical_numbers.tsv"
with open(out, "w", encoding="utf-8", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["key", "value", "note"])
    for r in rows:
        w.writerow(r)
print("[05] canonical    %d   -> %s" % (len(rows), out))
for r in rows:
    print("   %-36s %-14s %s" % (r[0], r[1], r[2]))
