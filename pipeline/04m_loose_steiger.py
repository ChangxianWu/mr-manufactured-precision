# -*- coding: utf-8 -*-
"""Steiger direction testing for the relaxed-threshold instrument sets."""
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

OUT = Path("D:/Horzontal_project/20B3/work/output")
DATA = Path("D:/Horzontal_project/20B3/work/data")
TESTS = OUT / "loose_clump_tests.csv"
NW = DATA / "finngen_R12_pheno_n.tsv"
OUTCSV = OUT / "loose_steiger.csv"
SUMCSV = OUT / "loose_steiger_summary.csv"


def main():
    d = pd.read_csv(TESTS, sep="\t")
    EP = {"SEPSIS": "AB1_OTHER_SEPSIS", "SEPTIC_STR": "AB1_STREPTO_SEPSIS",
          "SEPTIC_PNEU": "PNEUMO_SEPTHICHE", "PUERP_SEPS": "O15_PUERP_SEPSIS",
          "NEONAT_SEPS": "P16_BACTERIAL_SEPSIS_NEWBO", "PNEUMONIA": "J10_PNEUMONIA",
          "PNEUMOBACT": "J10_PNEUMOBACT", "VIRALPNEUMO": "J10_VIRALPNEUMO",
          "ARDS": "J10_ARDS", "PYELONEPHR": "N14_PYELONEPHR",
          "MENINGBACT": "G6_MENINGBACT", "ENDOCARD": "I9_ENDOCARD",
          "OSTEOMYEL": "M13_OSTEOMYELITIS", "CELLULITIS": "L12_CELLULITIS",
          "ABSCESS_CUT": "L12_ABSCESS_CUT", "ERYSIPELAS": "AB1_ERYSIPELAS",
          "INTEST_INF": "AB1_INTESTINAL_INFECTIONS", "GASTROENT": "AB1_GASTROENTERITIS_NOS",
          "BACT_OTHER": "AB1_OTHER_BACTERIAL", "BACT_NOS": "AB1_BACTINF_NOS",
          "ACTINOMYC": "AB1_ACTINOMYCOSIS", "COVID19": "U22_COVID19_CONFIRMED",
          "ZOSTER": "AB1_ZOSTER"}
    n = pd.read_csv(NW, sep="\t", encoding="utf-8-sig")
    nc = {}
    for _, r in n.iterrows():
        try:
            nc[r["phenocode"]] = (float(r["num_cases"]), float(r["num_controls"]))
        except (TypeError, ValueError):
            continue

    rows = []
    for r2 in sorted(d["r2_max"].unique()):
        s = d[d["r2_max"] == r2].copy()
        s["q"] = __import__("scipy.stats", fromlist=["false_discovery_control"]) \
            .false_discovery_control(s["p"].values, method="bh")
        surv = s[s["q"] < 0.05]
        nl = surv[surv["in_lct"] == 0]
        for _, r in nl.iterrows():
            code = r["endpoint"]
            ph = EP.get(code)
            if ph is None or ph not in nc:
                continue
            cases, ctrls = nc[ph]
            N_out = cases + ctrls
            K = cases / N_out
            zK = norm.ppf(1 - K)
            phi = norm.pdf(zK)
            eaf = float(r["eaf1"]) if pd.notna(r["eaf1"]) else np.nan
            be, se_e = float(r["beta_exp1"]), float(r["se_exp1"])
            bo, se_o = float(r["b"]), float(r["se"])
            r2_exp = 2 * eaf * (1 - eaf) * be ** 2
            r2_obs = 2 * eaf * (1 - eaf) * bo ** 2
            r2_liab = r2_obs * K * (1 - K) / phi ** 2
            se_r2_exp = 4 * eaf * (1 - eaf) * abs(be) * se_e
            se_r2_obs = 4 * eaf * (1 - eaf) * abs(bo) * se_o
            se_r2_liab = se_r2_obs * K * (1 - K) / phi ** 2
            N_exp = float(r["N_exp"])
            diff = r2_exp * N_exp - r2_liab * N_out
            se_d = math.sqrt((se_r2_exp * N_exp) ** 2 + (se_r2_liab * N_out) ** 2)
            z = diff / se_d if se_d > 0 else np.nan
            p_one = norm.sf(z) if np.isfinite(z) else np.nan
            rows.append({"r2_max": r2, "endpoint": code, "name": r["name"], "rsid": r["rsid"],
                         "n_snp": int(r["n_snp"]), "k_eff_over_k": 1.0 if r["n_snp"] == 1 else np.nan,
                         "f": 1.0 if r["n_snp"] == 1 else np.nan,
                         "b_out": bo, "se_out": se_o, "p_mr": float(r["p"]),
                         "r2_exp": r2_exp, "r2_liab": r2_liab, "N_exp": N_exp, "N_out": N_out,
                         "steiger_z": z, "steiger_p_one_sided": p_one,
                         "direction": ("tool_to_outcome" if (np.isfinite(z) and p_one < 0.05)
                                       else ("outcome_to_tool" if (np.isfinite(z) and p_one > 0.95)
                                             else "indeterminate"))})
    out = pd.DataFrame(rows)
    out.to_csv(OUTCSV, index=False)

    n_total = len(d)
    n_nl = len(out)
    n_single = int((out["n_snp"] == 1).sum()) if n_nl else 0
    n_correct = int((out["direction"] == "tool_to_outcome").sum()) if n_nl else 0
    n_reverse = int((out["direction"] == "outcome_to_tool").sum()) if n_nl else 0
    summ = pd.DataFrame([
        {"key": "loose_tests_total", "value": str(n_total), "note": "         MR     552×2 "},
        {"key": "loose_nonlct_survivors", "value": str(n_nl), "note": "      BH             "},
        {"key": "loose_nonlct_single_snp", "value": str(n_single),
         "note": "   n_snp=1    ——               k_eff/k=f=1 "},
        {"key": "loose_nonlct_direction_ok", "value": str(n_correct),
         "note": "      tool_to_outcome    p<0.05    "},
        {"key": "loose_nonlct_direction_reverse", "value": str(n_reverse),
         "note": "       outcome_to_tool    p>0.95    ——     "},
        {"key": "loose_nonlct_names",
         "value": "; ".join(sorted(out["name"].unique())) if n_nl else " ",
         "note": "           taxon"},
        {"key": "loose_coloc_status", "value": "deferred_to_04o",
         "note": "       04o      loose_coloc_summary.csv            "},
    ])
    summ.to_csv(SUMCSV, index=False)
    print("[loose-steiger]         %d       %d       %d    %d" % (n_nl, n_single, n_correct, n_reverse))
    if n_nl:
        print(out[["r2_max", "endpoint", "name", "n_snp", "b_out", "steiger_p_one_sided",
                   "direction"]].to_string(index=False))


if __name__ == "__main__":
    main()
