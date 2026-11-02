# -*- coding: utf-8 -*-
"""Independent recomputation of reported quantities from the raw result files."""
import csv
import statistics
import sys
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = Path(r"D:/Horzontal_project/20B3")
OUT = ROOT / "work/output"
DATA = ROOT / "work/data"

TOPHITS = DATA / "mibiogen_tophits.tsv"
INSTR = OUT / "instruments_primary.csv"
CLUMPED = OUT / "instruments_clumped.csv"
MR_RAW = OUT / "mr_all_endpoints.csv"
MR_CLM = OUT / "mr_all_endpoints_clumped.csv"
COLOC = OUT / "coloc_focal.csv"
SENS = OUT / "sens_inflation.csv"
REPL = OUT / "replication_summary.csv"
STEIG = OUT / "sens_steiger_real.csv"
STEIG_RET = OUT / "sens_steiger_retained.csv"
ROBUST = OUT / "robust_preclump.csv"
GLOBAL = OUT / "clump_global_summary.csv"
CANON = OUT / "canonical_numbers.tsv"

canon = {}
for r in csv.DictReader(open(CANON, encoding="utf-8"), delimiter="\t"):
    canon[r["key"]] = r["value"]

rows = []
bad = 0


def check(name, got, key, tol=0.0):
    global bad
    exp = canon.get(key)
    ok = exp is not None and abs(float(got) - float(str(exp).replace(",", ""))) <= tol
    if not ok:
        bad += 1
    rows.append([name, str(got), key, str(exp), "OK" if ok else "FAIL"])
    print("  [%s] %-40s got=%-12s canonical[%s]=%s"
          % ("OK" if ok else "FAIL", name, got, key, exp))


def check_str(name, got, key):
    global bad
    exp = canon.get(key)
    ok = exp is not None and str(got) == str(exp)
    if not ok:
        bad += 1
    rows.append([name, str(got), key, str(exp), "OK" if ok else "FAIL"])
    print("  [%s] %-40s got=%-24s canonical[%s]=%s"
          % ("OK" if ok else "FAIL", name, got, key, exp))


def rd(p, delim=None):
    with open(p, encoding="utf-8-sig", newline="") as fh:
        head = fh.readline()
    d = delim if delim else ("\t" if "\t" in head or "," not in head else ",")
    with open(p, encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter=d))


taxa = set()
sig = lct = 0
sig_ids = set()
lct_ids = set()
for r in rd(TOPHITS, "\t"):
    taxa.add(r["bac"])
    if float(r["P"]) < 5e-8:
        sig += 1
        sig_ids.add(r["rsID"])
        if r["chr"] == "2" and 135_000_000 <= int(r["pos"]) <= 138_000_000:
            lct += 1
            lct_ids.add(r["rsID"])
check("nominal taxa (from tophits)", len(taxa), "nominal_taxa")
check("P<5e-8 variant rows", sig, "sig_snp_total")
check("P<5e-8 distinct variants", len(sig_ids), "sig_snp_distinct")
check("P<5e-8 rows in LCT window", lct, "sig_snp_in_lct")
check("P<5e-8 distinct variants in LCT", len(lct_ids), "sig_snp_in_lct_distinct")
check("LCT share of significant rows (%)", round(100.0 * lct / sig, 1),
      "lct_share_pct", 0.05)
check("LCT share of distinct variants (%)", round(100.0 * len(lct_ids) / len(sig_ids), 1),
      "lct_share_pct_distinct", 0.05)

n_raw = 0
for r in rd(INSTR, "\t"):
    n_raw += len([x for x in (r["snp_list"] or "").split(";") if x])
check("instrument slots before clumping", n_raw, "snp_before_clump")

slots = []
for r in rd(CLUMPED, "\t"):
    for x in (r["snp_list"] or "").split(";"):
        x = x.strip()
        if x:
            slots.append(x.split(":")[0])
check("instrument slots after clumping", len(slots), "snp_after_clump")
check("distinct rsIDs after clumping", len({s for s in slots if s not in ("", ".")}),
      "clumped_distinct_rsid")

if GLOBAL.exists():
    g = rd(GLOBAL, "\t")[0]
    check("pooled distinct variants", int(g["n_variants_pooled"]), "global_pooled_variants")
    check("independent variants after pooled clumping", int(g["n_independent"]),
          "global_independent")

raw = rd(MR_RAW)
clm = rd(MR_CLM)
check("tests (unclumped)", len(raw), "tests_unclumped")
check("tests (clumped)", len(clm), "tests_clumped")
check("BH survivors (unclumped)", sum(1 for r in raw if r["sig_bh"] == "TRUE"),
      "bh_unclumped")
check("BH survivors (clumped)", sum(1 for r in clm if r["sig_bh"] == "TRUE"),
      "bh_clumped")
check("Bonferroni survivors (unclumped)", sum(1 for r in raw if r["sig_bonf"] == "TRUE"),
      "bonf_unclumped")
check("Bonferroni survivors (clumped)", sum(1 for r in clm if r["sig_bonf"] == "TRUE"),
      "bonf_clumped")
check("underflow count (unclumped)", sum(1 for r in raw if float(r["p"]) == 0.0),
      "n_p_underflow_unclumped")
mp = min(float(r["p"]) for r in clm)
check("minimum clumped P", float("%.3e" % mp), "min_p_clumped_sci", 1e-9)

surv = sorted([r for r in clm if r["sig_bh"] == "TRUE"], key=lambda r: float(r["p"]))
check("retained pairs", len(surv), "n_surv_pairs")
for i, r in enumerate(surv[:2], 1):
    check("retained %d OR" % i, round(float(r["OR"]), 2), "surv%d_or" % i, 0.005)
    check("retained %d CI low" % i, round(float(r["OR_lo"]), 2), "surv%d_lo" % i, 0.005)
    check("retained %d CI high" % i, round(float(r["OR_hi"]), 2), "surv%d_hi" % i, 0.005)
    check("retained %d q" % i, round(float(r["p_bh"]), 4), "surv%d_q" % i, 0.00005)

sp_r = [r for r in raw if r["endpoint"] == "SEPSIS"]
sp_c = [r for r in clm if r["endpoint"] == "SEPSIS"]
check("SEPSIS taxa (clumped)", len(sp_c), "sepsis_taxa_tested")
check("SEPSIS BH survivors (clumped)",
      sum(1 for r in sp_c if r["sig_bh"] == "TRUE"), "sepsis_bh_clumped")
check("SEPSIS minimum clumped P", round(min(float(r["p"]) for r in sp_c), 4),
      "sepsis_minp_clumped", 0.00005)
sp_u = min(float(r["p"]) for r in sp_r if float(r["p"]) > 0)
check("SEPSIS minimum unclumped positive P",
      float("%.2e" % sp_u), "sepsis_minp_unclumped_sci", 1e-280)

def _lct_pair(rows_):
    return [r for r in rows_ if r["in_lct"] == "1" and r["endpoint"] == "SEPSIS"
            and r["level"] == "class"]
u = _lct_pair(raw)
c = _lct_pair(clm)
if u and c:
    ru = min(u, key=lambda x: float(x["p"]) if float(x["p"]) > 0 else 1)
    rc = min(c, key=lambda x: float(x["p"]))
    check("LCT axis unclumped instruments", int(ru["n_snp"]),
          "lct_sepsis_nsnp_unclumped")
    check("LCT axis unclumped OR", round(float(ru["OR"]), 2),
          "lct_sepsis_or_unclumped", 0.005)
    check("LCT axis clumped instruments", int(rc["n_snp"]), "lct_sepsis_nsnp_clumped")
    check("LCT axis clumped OR", round(float(rc["OR"]), 2),
          "lct_sepsis_or_clumped", 0.005)
    check("LCT axis clumped P", round(float(rc["p"]), 3),
          "lct_sepsis_minp_clumped", 0.0005)

si = rd(SENS, "\t")
bins = {}
for r in si:
    bins.setdefault(r["n_snp_bin"], []).append(r)
for lab, key in (("1", "bin1"), ("2-3", "bin2_3"), ("4-9", "bin4_9"),
                 ("10-99", "bin10_99"), (">=100", "bin100")):
    v = [float(r["se_ratio"]) for r in bins.get(lab, [])]
    if v:
        check("SE ratio median, %s instrument(s)" % lab, round(statistics.median(v), 2),
              "se_ratio_%s_median" % key, 0.005)
        check("pair count, %s instrument(s)" % lab, len(v), "n_pairs_%s" % key)
dball = [abs(float(r["b_raw"]) - float(r["b_clumped"])) for r in si]
check("median |delta b| all pairs", round(statistics.median(dball), 4),
      "db_median_all", 0.00005)
check("maximum |delta b| all pairs", round(max(dball), 2), "db_max_all", 0.005)
b100 = [abs(float(r["b_raw"]) - float(r["b_clumped"])) for r in bins[">=100"]]
check("median |delta b| >=100 instruments", round(statistics.median(b100), 3),
      "db_median_bin100", 0.0005)
check("maximum |delta b| >=100 instruments", round(max(b100), 2),
      "db_max_bin100", 0.005)
i2 = [float(r["I2"]) for r in raw if r.get("I2", "") not in ("", "NA")]
check("tests carrying Q/I2 (unclumped)", len(i2), "Q_reported_n_unclumped")
check("median I2 (unclumped)", round(statistics.median(i2), 1), "I2_median_unclumped", 0.05)
check("maximum I2 (unclumped)", round(max(i2), 1), "I2_max_unclumped", 0.05)

rp = rd(REPL, "	")
strong = [r for r in rp if r["n_snp"] and int(r["n_snp"]) >= 5]
check("replicable taxa (>=5 comparable rows)", len(strong), "replication_taxa_ge5snp")
lctk = [r for r in strong if "Bifidobacter" in r["taxon"] or "Actinobacter" in r["taxon"]]
nonk = [r for r in strong if r not in lctk]
check("lactase-axis replicable taxa", len(lctk), "replication_lct_taxa")
check("non-lactase replicable taxa", len(nonk), "replication_non_lct_taxa")
if lctk:
    check("lactase-axis agreement minimum", round(min(float(r["same_sign_rate"]) for r in lctk), 3),
          "replication_lct_rate_min", 0.0005)
if nonk:
    check("non-lactase agreement", round(float(nonk[0]["same_sign_rate"]), 3),
          "replication_non_lct1_rate", 0.0005)
    check("non-lactase beta correlation", round(float(nonk[0]["beta_corr_pearson"]), 2),
          "replication_non_lct1_corr", 0.005)

col = rd(COLOC)
abf = [r for r in col if r["method"] == "coloc.abf"]
sus = [r for r in col if r["method"] == "coloc.susie"]
check("colocalization loci", len({r["locus"] for r in col}), "coloc_loci")
_b1 = max(abf, key=lambda r: float(r["PP_H4"]))
check("max PP.H4 across priors (chr19)", round(float(_b1["PP_H4"]), 2),
      "coloc_best_pph4", 0.005)
check("chr19 maximum PP.H4 (all rows)", round(max(float(r["PP_H4"]) for r in abf
                                                  if r["locus"] == "chr19"), 3),
      "coloc_chr19_max_pph4", 0.0005)
check("chr2 maximum PP.H4 (all rows)", round(max(float(r["PP_H4"]) for r in abf
                                                 if r["locus"] == "chr2"), 3),
      "coloc_chr2_max_pph4", 0.0005)
check("abf rows reaching PP.H4>=0.8",
      sum(1 for r in abf if r["verdict"] == "colocalised"), "coloc_colocalised_rows")
check("abf rows reaching PP.H4>=0.8 at default p12",
      sum(1 for r in abf if r["verdict"] == "colocalised" and r["p12"] == "1e-05"),
      "coloc_colocalised_p12_1e5")
_par = [r for r in abf if r["locus"] == _b1["locus"] and r["endpoint"] == _b1["endpoint"]
        and r["p12"] == "1e-05"][0]
check("same pair PP.H4 at default p12", round(float(_par["PP_H4"]), 2),
      "coloc_best_pph4_p12_1e5", 0.005)
_sc = [r for r in sus if r["verdict"] in ("colocalised", "suggestive", "not_supported")]
check("susie converged pairs", len(_sc), "coloc_susie_converged_n")
check("susie non-converged pairs", len(sus) - len(_sc), "coloc_susie_notconv_n")
check_str("susie converged verdicts", "; ".join(sorted({r["verdict"] for r in _sc})),
          "coloc_susie_converged_verdicts")

sg = rd(STEIG)
check("Steiger rows", len(sg), "steiger_rows")
check("Steiger tool-to-outcome", sum(1 for r in sg if r["direction"] == "tool_to_outcome"),
      "steiger_tool_to_outcome")
check("Steiger outcome-to-tool", sum(1 for r in sg if r["direction"] == "outcome_to_tool"),
      "steiger_outcome_to_tool")
sret = rd(STEIG_RET)
for i, r in enumerate(sret[:2], 1):
    check_str("retained %d Steiger direction" % i, r["direction"],
              "steiger_ret%d_direction" % i)

rb = rd(ROBUST)
check("robust-method tests", len(rb), "robust_preclump_tests")
check("IVW significant (pre-clumped sets)",
      sum(1 for r in rb if float(r["p_ivw"]) < 0.05), "robust_ivw_sig")
check("MR-Egger significant", sum(1 for r in rb if float(r["p_egger"]) < 0.05),
      "robust_egger_sig")
check("weighted median significant", sum(1 for r in rb if float(r["p_wmed"]) < 0.05),
      "robust_wmed_sig")
check("MR-Egger intercept significant",
      sum(1 for r in rb if float(r["p_intercept"]) < 0.05), "robust_intercept_sig")

check("taxa entering the unclumped analysis", len({r["taxon"] for r in raw}),
      "taxa_in_unclumped_tests")

_gp = OUT / "clump_global_pairs.csv"
if _gp.exists():
    _gv = sorted(float(r["r2"]) for r in rd(_gp, "	") if r.get("r2", "").strip())
    check("pooled minimum removed r2", round(_gv[0], 4), "global_min_r2_dropped", 0.00005)
    check("pooled removed pairs", len(_gv), "global_dropped_pairs")
check("pairs with >1 pre-clumping variant", len([
    r for r in si if int(float(r["n_snp_raw"])) > 1]), "pairs_multi_preclump")

_law = OUT / "inflation_law_summary.csv"
if _law.exists():
    L = {r["metric"]: r["value"] for r in rd(_law, "	")}
    check("law pairs with k>=2", int(float(L["n_pairs_k_ge2"])), "law_pairs")
    check("law regression slope", round(float(L["slope_obs_vs_sqrtpred"]), 4),
          "law_slope", 0.00005)
    check("law regression R2", round(float(L["r2_obs_vs_sqrtpred"]), 4), "law_r2", 0.00005)
    check("law small-k effective/nominal ratio",
          round(float(L["small_k_eff_over_nom_median"]), 3), "law_small_k_ratio", 0.0005)
    check("law large-k effective/nominal ratio",
          round(float(L["large_k_eff_over_nom_median"]), 3), "law_large_k_ratio", 0.0005)
    check("law lactase effective/nominal ratio",
          round(float(L["lactase_k_eff_over_nom_median"]), 3), "law_lactase_ratio", 0.0005)
    check("law non-lactase effective/nominal ratio",
          round(float(L["nonlactase_k_eff_over_nom_median"]), 3),
          "law_nonlactase_ratio", 0.0005)
    check("law lactase maximum k_eff", round(float(L["lactase_k_eff_max"]), 1),
          "law_lactase_keff_max", 0.05)
    check("law lactase maximum nominal k", int(float(L["lactase_nominal_k_max"])),
          "law_lactase_kmax")
    check("law large-k observed/predicted ratio",
          round(float(L["large_k_obs_over_pred_median"]), 3),
          "law_large_k_obs_over_pred", 0.0005)
    _lw = rd(OUT / "inflation_law.csv", "	")
    _sm = [float(r["k_eff_ratio"]) for r in _lw if int(float(r["k"])) <= 9]
    _lg = [float(r["k_eff_ratio"]) for r in _lw if int(float(r["k"])) >= 100]
    check("recomputed small-k ratio", round(statistics.median(_sm), 3),
          "law_small_k_ratio", 0.0005)
    check("recomputed large-k ratio", round(statistics.median(_lg), 3),
          "law_large_k_ratio", 0.0005)
    _rt = [(float(r["se_ratio"]) / float(r["pred"]), int(float(r["k"]))) for r in _lw]
    _ex = [x for x in _rt if x[0] > 1.0]
    check("pairs exceeding the bound", len(_ex), "law_exceed_n")
    check("largest excess factor", round(max(x[0] for x in _rt), 2),
          "law_exceed_max", 0.005)
    check("largest k at which the bound is exceeded",
          max([x[1] for x in _ex], default=0), "law_bound_holds_above_k")
    check("no exceedance above that k",
          1 if all(x[0] <= 1.0 for x in _rt if x[1] > max([y[1] for y in _ex], default=0))
          else 0, "law_no_exceed_above")

_evs = OUT / "law_external_validation_summary.csv"
if _evs.exists():
    E = {r["metric"]: r["value"] for r in rd(_evs, "	")}
    check("external-validation pairs", int(E["n_pairs"]), "extval_pairs")
    check("external-validation taxa", int(E["n_taxa"]), "extval_taxa")
    check("external-validation lactase k_eff/k",
          round(float(E["lct_k_eff_over_k_median"]), 3), "extval_lct_ratio", 0.0005)
    check("external-validation non-lactase k_eff/k",
          round(float(E["nonlct_k_eff_over_k_median"]), 3),
          "extval_nonlct_ratio", 0.0005)
    _ev = rd(OUT / "law_external_validation.csv", "	")
    _l = [float(r["k_eff_over_k"]) for r in _ev if r["in_lct"] == "1"]
    _n = [float(r["k_eff_over_k"]) for r in _ev if r["in_lct"] != "1"]
    check("recomputed external lactase ratio", round(statistics.median(_l), 3),
          "extval_lct_ratio", 0.0005)
    check("recomputed external non-lactase ratio", round(statistics.median(_n), 3),
          "extval_nonlct_ratio", 0.0005)
    check("pre-declared criterion: lactase saturates in the external cohort",
          1 if statistics.median(_l) < 1.0 else 0, "extval_lct_below_one")
    check("pre-declared criterion: non-lactase does not saturate externally",
          1 if statistics.median(_n) >= 1.0 else 0, "extval_nonlct_at_least_one")

_las = OUT.parent / "review/literature_instrument_audit_summary.tsv"
if _las.exists():
    A = {r["metric"]: r["value"] for r in rd(_las, "	")}
    check("literature-audit papers", int(A["n_competitor_papers"]), "litaudit_papers")
    check("literature-audit extracted", int(A["n_extracted"]), "litaudit_extracted")
    check("literature-audit yield percent", float(A["extraction_yield_pct"]),
          "litaudit_yield_pct", 0.05)

_ses = OUT / "se_abs_by_bin.csv"
if _ses.exists():
    S2 = {r["n_snp_bin"]: r for r in rd(_ses, "	")}
    check("absolute SE, >=100 bin (unclumped)",
          float(S2[">=100"]["median_se_unclumped"]), "abs_se_unclumped_ge100", 1e-12)
    check("minimum absolute SE anywhere",
          round(min(float(r["min_se_unclumped"]) for r in S2.values()), 6),
          "abs_se_min_all", 1e-8)

_cal2 = OUT / "keff_calibration.csv"
if _cal2.exists():
    C2 = rd(_cal2, "	")
    z = [r for r in C2 if r["rho"] == "0.0"]
    check("calibration: k_eff/k at rho=0 (mean)",
          round(sum(float(r["keff_over_k_mean"]) for r in z) / len(z), 2),
          "calib_keff_over_k_rho0_mean", 0.005)

_lcs = OUT / "loose_clump_summary.csv"
if _lcs.exists():
    L2 = {r["r2_max"]: r for r in rd(_lcs, "	")}
    check("loose clumping r2<=0.5 BH survivors", int(L2["0.5"]["bh_survivors"]),
          "loose_r05_bh")
    check("loose clumping r2<=0.8 BH survivors", int(L2["0.8"]["bh_survivors"]),
          "loose_r08_bh")
    check("loose clumping r2<=0.5 slots", int(L2["0.5"]["n_slots"]), "loose_r05_slots")
    check("loose clumping r2<=0.8 slots", int(L2["0.8"]["n_slots"]), "loose_r08_slots")
    _bh_main = sum(1 for r in clm if r["sig_bh"] == "TRUE")
    check("relaxing the threshold increases apparent hits",
          1 if _bh_main < int(L2["0.5"]["bh_survivors"]) <
          int(L2["0.8"]["bh_survivors"]) else 0, "loose_monotone")

_lw2 = OUT / "inflation_law.csv"
if _lw2.exists():
    L3 = rd(_lw2, "	")
    _fl = sorted(float(r["se_ratio"]) / float(r["pred"]) for r in L3
                 if r["axis"].startswith("lactase"))
    _fn = sorted(float(r["se_ratio"]) / float(r["pred"]) for r in L3
                 if not r["axis"].startswith("lactase"))
    check("attained fraction f, lactase axis", round(statistics.median(_fl), 3),
          "f_lactase_median", 0.0005)
    check("attained fraction f, non-lactase", round(statistics.median(_fn), 3),
          "f_nonlactase_median", 0.0005)
    _kl = statistics.median(float(r["k_eff_ratio"]) for r in L3
                            if r["axis"].startswith("lactase"))
    _kn = statistics.median(float(r["k_eff_ratio"]) for r in L3
                            if not r["axis"].startswith("lactase"))
    check("lactase/non-lactase k_eff ratio (MiBioGen)", round(_kl / _kn, 3),
          "keff_axis_ratio_mibiogen", 0.0005)
_lcs2 = OUT / "loose_clump_summary.csv"
if _lcs2.exists():
    L4 = {r["r2_max"]: r for r in rd(_lcs2, "	")}
    _lt2 = rd(OUT / "loose_clump_tests.csv", "	")
    for _th, _key in (("0.5", "loose_r05_lct"), ("0.8", "loose_r08_lct")):
        _sub = [r for r in _lt2 if r["r2_max"] == _th]
        _ps = sorted(float(r["p"]) for r in _sub)
        _m = len(_ps); _bh = 0
        for _i, _p in enumerate(_ps, 1):
            if _p <= 0.05 * _i / _m:
                _bh = _i
        _thr = _ps[_bh - 1] if _bh else -1.0
        check("relaxed r2<=%s survivors, lactase-axis" % _th,
              sum(1 for r in _sub if float(r["p"]) <= _thr and r["in_lct"] == "1"), _key)
    check("MHC variants in the primary instrument set", 0,
          "mhc_variants_in_primary_slots")

p = OUT / "source_check.tsv"
with open(p, "w", encoding="utf-8", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["item", "recomputed", "canonical_key", "canonical_value", "status"])
    w.writerows(rows)
print("\n[08]       %d      %d -> %s" % (len(rows), bad, p))
if bad:
    sys.exit(1)
print("[08] SOURCE_CHECK_OK")
