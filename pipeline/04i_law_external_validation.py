# -*- coding: utf-8 -*-
"""External validation of the law with exposure effects re-estimated in the independent cohort; two pre-declared criteria."""
import csv
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = Path(r"D:/Horzontal_project/20B3")
OUT = ROOT / "work/output"
DATA = ROOT / "work/data"

COMPL = {"A": "T", "T": "A", "C": "G", "G": "C"}


def is_pal(a, b):
    return a in COMPL and b in COMPL and COMPL[a] == b


def rd(path, delim=None):
    with open(path, encoding="utf-8-sig", newline="") as fh:
        head = fh.readline()
    d = delim or ("\t" if "\t" in head else ",")
    with open(path, encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter=d))


def ivw(pairs):
    """pairs: [(ratio, se_ratio)] -> (b, se)      IVW """
    w = [1.0 / (sr * sr) for _, sr in pairs if sr and sr > 0]
    if len(w) != len(pairs) or not w:
        return None
    b = sum(wi * r for wi, (r, _) in zip(w, pairs)) / sum(w)
    return b, math.sqrt(1.0 / sum(w))


def align(finn_ea, finn_oa, ref, alt, beta_out, se_out):
    """    EA     FinnGen               """
    if not all([finn_ea, finn_oa, ref, alt]):
        return None
    e, o, r, a = finn_ea.upper(), finn_oa.upper(), ref.upper(), alt.upper()
    if is_pal(e, o):
        return None
    if e == a and o == r:
        return beta_out, se_out
    if e == r and o == a:
        return -beta_out, se_out
    return None


def main():
    prim = rd(OUT / "instruments_primary.csv", "\t")
    clm = rd(OUT / "instruments_clumped.csv", "\t")
    prim_k = {r["taxon"]: len([x for x in (r["snp_list"] or "").split(";") if x.strip()])
              for r in prim}
    clm_set = {r["taxon"]: {x.split(":")[0] for x in (r["snp_list"] or "").split(";") if x.strip()}
               for r in clm}
    clm_k = {t: len(v) for t, v in clm_set.items()}

    finn = {}
    for r in rd(OUT / "replication_exposure.csv", "\t"):
        try:
            b = float(r["beta_finnish"]); s = float(r["se_finnish"])
        except (TypeError, ValueError):
            continue
        if not (s > 0):
            continue
        finn.setdefault((r["taxon"], r["rsid"]), (b, s, r["ea_finnish"], r["oa_finnish"]))

    taxa = sorted({t for t, _ in finn})
    print("[04i]        taxa %d       %d" % (len(taxa), len(finn)))

    rows = []
    for ep_file in sorted(DATA.glob("outcome_extract/outcome_extract_*.tsv")):
        key = ep_file.stem.replace("outcome_extract_", "")
        O = rd(ep_file, "\t")
        for r in O:
            t = r["taxon"]
            if t not in taxa:
                continue
            f = finn.get((t, r["rsid"]))
            if f is None:
                continue
            beta_f, se_f, ea_f, oa_f = f
            if beta_f == 0:
                continue
            try:
                beta_o = float(r["beta_out"]); se_o = float(r["se_out"])
            except (TypeError, ValueError):
                continue
            if not (se_o > 0):
                continue
            al = align(ea_f, oa_f, r["ref_out"], r["alt_out"], beta_o, se_o)
            if al is None:
                continue
            beta_o_al, se_o_al = al
            ratio = beta_o_al / beta_f
            se_ratio_snp = se_o_al / abs(beta_f)
            rows.append({"endpoint": key, "taxon": t, "name": r["name"],
                         "level": r["level"], "rsid": r["rsid"],
                         "in_lct": r["in_lct"], "ratio": ratio,
                         "se_snp": se_ratio_snp,
                         "clumped": 1 if r["rsid"] in clm_set.get(t, set()) else 0})

    res = []
    by = defaultdict(list)
    for r in rows:
        by[(r["endpoint"], r["taxon"])].append(r)
    for (ep, t), g in by.items():
        if len(g) < 2:
            continue
        u = ivw([(x["ratio"], x["se_snp"]) for x in g])
        cg = [x for x in g if x["clumped"] == 1]
        if len(cg) < 1 or u is None:
            continue
        c = ivw([(x["ratio"], x["se_snp"]) for x in cg])
        if c is None:
            continue
        b_u, se_u = u
        b_c, se_c = c
        k_fin = len(g)
        r_obs = se_c / se_u
        res.append({
            "endpoint": ep, "taxon": t, "name": g[0]["name"], "level": g[0]["level"],
            "in_lct": g[0]["in_lct"],
            "k_finnish_matched": k_fin,
            "k_primary": prim_k.get(t, ""), "k_clumped": clm_k.get(t, ""),
            "b_unclumped": b_u, "se_unclumped": se_u,
            "b_clumped": b_c, "se_clumped": se_c,
            "se_ratio_observed": r_obs,
            "pred_sqrt_k": math.sqrt(k_fin),
            "k_eff": r_obs ** 2,
            "k_eff_over_k": (r_obs ** 2) / k_fin,
        })

    if not res:
        raise SystemExit("[04i]       Finnish      ")
    res.sort(key=lambda x: (x["in_lct"] != "1", -x["k_finnish_matched"]))
    with open(OUT / "law_external_validation.csv", "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(res[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(res)

    def med(v):
        return statistics.median(v) if v else float("nan")

    lct = [x["k_eff_over_k"] for x in res if x["in_lct"] == "1"]
    non = [x["k_eff_over_k"] for x in res if x["in_lct"] != "1"]
    S = [
        ("n_pairs", len(res)),
        ("n_taxa", len({x["taxon"] for x in res})),
        ("lct_pairs", len(lct)),
        ("lct_k_eff_over_k_median", "%.3f" % med(lct)),
        ("nonlct_pairs", len(non)),
        ("nonlct_k_eff_over_k_median", "%.3f" % med(non)),
        ("lct_k_max", max(x["k_finnish_matched"] for x in res if x["in_lct"] == "1")),
        ("se_ratio_median_lct", "%.2f" % med([x["se_ratio_observed"] for x in res if x["in_lct"] == "1"])),
        ("pred_sqrt_k_median_lct", "%.2f" % med([x["pred_sqrt_k"] for x in res if x["in_lct"] == "1"])),
        ("se_ratio_median_nonlct", "%.2f" % med([x["se_ratio_observed"] for x in res if x["in_lct"] != "1"])),
        ("pred_sqrt_k_median_nonlct", "%.2f" % med([x["pred_sqrt_k"] for x in res if x["in_lct"] != "1"])),
    ]
    with open(OUT / "law_external_validation_summary.csv", "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["metric", "value"])
        w.writerows(S)
    for k, v in S:
        print("  %-34s %s" % (k, v))
    print("[04i] -> law_external_validation*.csv")


if __name__ == "__main__":
    main()
