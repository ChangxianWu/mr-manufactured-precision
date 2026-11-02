# -*- coding: utf-8 -*-
"""LD-aware GLS estimator (archived for transparency; results not reported)."""
import csv
import io
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

P_INSTR = Path("D:/Horzontal_project/20B3/work/output/instruments_primary.csv")
P_TOPHITS = Path("D:/Horzontal_project/20B3/work/data/mibiogen_tophits.tsv")
P_GENO = Path("D:/Horzontal_project/20B3/work/cache/eur_geno_cache.tsv")
P_MR = Path("D:/Horzontal_project/20B3/work/output/mr_all_endpoints.csv")
P_SENS = Path("D:/Horzontal_project/20B3/work/output/sens_inflation.csv")
D_EXTRACT = "D:/Horzontal_project/20B3/work/data/outcome_extract"
OUT_GLS = Path("D:/Horzontal_project/20B3/work/output/ld_aware_gls.csv")
OUT_SUM = Path("D:/Horzontal_project/20B3/work/output/ld_aware_gls_summary.csv")
OUT_CAL = Path("D:/Horzontal_project/20B3/work/output/keff_calibration.csv")
OUT_SE = Path("D:/Horzontal_project/20B3/work/output/se_abs_by_bin.csv")

COMPL = {"A": "T", "T": "A", "C": "G", "G": "C"}
ENDPOINTS = ["SEPSIS", "SEPTIC_STR", "SEPTIC_PNEU", "PUERP_SEPS", "NEONAT_SEPS",
             "PNEUMONIA", "PNEUMOBACT", "VIRALPNEUMO", "ARDS", "PYELONEPHR",
             "MENINGBACT", "ENDOCARD", "OSTEOMYEL", "CELLULITIS", "ABSCESS_CUT",
             "ERYSIPELAS", "INTEST_INF", "GASTROENT", "BACT_OTHER", "BACT_NOS",
             "ACTINOMYC", "COVID19", "ZOSTER"]


def rd(path, delim=None):
    text = Path(path).read_text(encoding="utf-8-sig")
    head = text.split("\n", 1)[0]
    d = delim or ("\t" if "\t" in head else ",")
    return list(csv.DictReader(io.StringIO(text), delimiter=d))


def is_pal(a, b):
    return a in COMPL and b in COMPL and COMPL[a] == b


def align(ea, oa, ref, alt, beta_out, se_out):
    """                      (beta, se)   None """
    if not all([ea, oa, ref, alt]):
        return None
    e, o, r, a = ea.upper(), oa.upper(), ref.upper(), alt.upper()
    if is_pal(e, o):
        return None
    if e == a and o == r:
        return beta_out, se_out
    if e == r and o == a:
        return -beta_out, se_out
    return None


def load_genotypes():
    """(chr,pos) -> np.array        2*MAF    """
    geno = {}
    lines = P_GENO.read_text(encoding="utf-8").splitlines()
    for ln in lines[1:]:
        f = ln.split("\t")
        if len(f) < 3:
            continue
        raw = f[2].split(",")
        vals = np.array([np.nan if x == "." else float(x) for x in raw])
        if np.isnan(vals).any():
            maf = np.nanmean(vals) / 2.0
            vals = np.where(np.isnan(vals), 2.0 * maf, vals)
        geno[(int(f[0]), str(f[1]))] = vals
    return geno


def ld_matrix(doses):
    """ =    Pearson r           0        """
    k = len(doses)
    R = np.eye(k)
    for i in range(k):
        for j in range(i + 1, k):
            a, b = doses[i], doses[j]
            if a.var() <= 0 or b.var() <= 0:
                continue
            r = float(np.corrcoef(a, b)[0, 1])
            if not math.isfinite(r):
                r = 0.0
            R[i, j] = R[j, i] = r
    return R


def naive_ivw(x, y, se):
    r = y / x
    w = (x * x) / (se * se)
    return float((w * r).sum() / w.sum()), math.sqrt(1.0 / w.sum())


def gls_ivw(x, y, se, R):
    """V = diag(se) R diag(se)                       """
    V = np.outer(se, se) * R
    Vi = np.linalg.pinv(V, rcond=1e-10)
    d = float(x @ Vi @ x)
    if not math.isfinite(d) or d <= 0:
        return None
    return float((x @ Vi @ y) / d), math.sqrt(1.0 / d)


def norm_sf2(z):
    return 2.0 * (1.0 - 0.5 * (1.0 + math.erf(abs(z) / math.sqrt(2.0))))


def main():
    prim = rd(P_INSTR, "\t")
    top = rd(P_TOPHITS, "\t")
    pos_of = {}
    for r in top:
        try:
            if float(r["P"]) < 5e-8:
                pos_of[r["rsID"]] = (int(r["chr"]), str(r["pos"]))
        except (TypeError, ValueError):
            continue
    geno = load_genotypes()
    print("[04k]       %d    tophits         %d"
          % (len(geno), len(pos_of)))

    instr = {}
    for r in prim:
        items = []
        for tok in (r["snp_list"] or "").split(";"):
            tok = tok.strip()
            if not tok:
                continue
            p = tok.split(":")
            if len(p) < 3:
                continue
            try:
                items.append((p[0], float(p[1]), float(p[2])))
            except ValueError:
                continue
        instr[r["taxon"]] = {"name": r["name"], "level": r["level"],
                             "in_lct": r["in_lct"], "items": items}

    mr = {(r["taxon"], r["endpoint"]): r for r in rd(P_MR)}
    rows = []
    for ep in ENDPOINTS:
        f = "%s/outcome_extract_%s.tsv" % (D_EXTRACT, ep)
        if not Path(f).exists():
            continue
        by_taxon = defaultdict(list)
        for r in rd(f, "\t"):
            t = r["taxon"]
            if t not in instr:
                continue
            try:
                bo, so = float(r["beta_out"]), float(r["se_out"])
            except (TypeError, ValueError):
                continue
            if not (so > 0):
                continue
            al = align(r["ea_instr"], r["oa_instr"], r["ref_out"], r["alt_out"], bo, so)
            if al is None:
                continue
            by_taxon[t].append((r["rsid"], al[0], al[1]))
        for t, hits in by_taxon.items():
            meta = instr[t]
            beta_x = {rs: bx for rs, bx, _ in meta["items"]}
            se_x = {rs: sx for rs, bx, sx in meta["items"]}
            usable = []
            for rs, bo_al, so_al in hits:
                if rs not in beta_x:
                    continue
                key = pos_of.get(rs)
                if key is None or key not in geno:
                    continue
                usable.append((rs, beta_x[rs], se_x[rs], bo_al, so_al, key))
            if len(usable) < 2:
                continue
            xs = np.array([u[1] for u in usable], dtype=float)
            ys = np.array([u[3] for u in usable], dtype=float)
            ses = np.array([u[4] for u in usable], dtype=float)
            if np.any(xs == 0):
                continue
            R = ld_matrix([geno[u[5]] for u in usable])
            nv = naive_ivw(xs, ys, ses)
            gl = gls_ivw(xs, ys, ses, R)
            if nv is None or gl is None:
                continue
            b_n, se_n = nv
            b_g, se_g = gl
            key = mr.get((t, ep))
            rows.append({
                "endpoint": ep, "taxon": t, "name": meta["name"],
                "level": meta["level"], "in_lct": meta["in_lct"],
                "k_used": len(usable),
                "b_naive": b_n, "se_naive": se_n,
                "b_gls": b_g, "se_gls": se_g,
                "inflation_factor": se_g / se_n,
                "p_naive": norm_sf2(b_n / se_n),
                "p_gls": norm_sf2(b_g / se_g),
                "naive_bh_sig": key["sig_bh"] if key else "",
                "sqrt_k": math.sqrt(len(usable)),
            })

    if not rows:
        raise SystemExit("[04k]       ")
    with io.StringIO() as buf:
        w = csv.DictWriter(buf, fieldnames=list(rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
        OUT_GLS.write_text(buf.getvalue(), encoding="utf-8", newline="")

    def binof(k):
        return "2-3" if k <= 3 else ("4-9" if k <= 9 else ("10-99" if k <= 99 else ">=100"))

    bybin = defaultdict(list)
    for r in rows:
        bybin[binof(r["k_used"])].append(r)
    summ = []
    print("\n%-8s %5s %6s %12s %12s %10s %10s"
          % ("bin", "n", "med_k", "se_naive", "se_gls", "IF_med", "sqrt(k)"))
    for b in ["2-3", "4-9", "10-99", ">=100"]:
        v = bybin.get(b, [])
        if not v:
            continue
        mk = statistics.median(x["k_used"] for x in v)
        mn = statistics.median(x["se_naive"] for x in v)
        mg = statistics.median(x["se_gls"] for x in v)
        f = statistics.median(x["inflation_factor"] for x in v)
        print("%-8s %5d %6.0f %12.4e %12.4e %10.2f %10.2f"
              % (b, len(v), mk, mn, mg, f, math.sqrt(mk)))
        key = b.replace("-", "_").replace(">=", "ge")
        summ += [("n_pairs_" + key, len(v)),
                 ("median_k_" + key, "%.0f" % mk),
                 ("median_se_naive_" + key, "%.4e" % mn),
                 ("median_se_gls_" + key, "%.4e" % mg),
                 ("median_inflation_factor_" + key, "%.2f" % f),
                 ("median_sqrt_k_" + key, "%.2f" % math.sqrt(mk))]

    sig = [x for x in rows if x["naive_bh_sig"] == "TRUE"]
    still = [x for x in sig if x["p_gls"] < 0.05]
    summ += [
        ("n_pairs_with_ld", len(rows)),
        ("n_taxa", len({x["taxon"] for x in rows})),
        ("naive_bh_sig_pairs", len(sig)),
        ("naive_bh_sig_still_p05_under_gls", len(still)),
        ("overall_median_inflation_factor",
         "%.2f" % statistics.median(x["inflation_factor"] for x in rows)),
        ("max_inflation_factor", "%.1f" % max(x["inflation_factor"] for x in rows)),
    ]
    with io.StringIO() as buf:
        w = csv.writer(buf, delimiter="\t", lineterminator="\n")
        w.writerow(["metric", "value"])
        w.writerows(summ)
        OUT_SUM.write_text(buf.getvalue(), encoding="utf-8", newline="")
    print("\n     BH    %d   LD     p<0.05     %d  "
          % (len(sig), len(still)))
    print("  IF    %.2f    %.1f"
          % (statistics.median(x["inflation_factor"] for x in rows),
             max(x["inflation_factor"] for x in rows)))

    rng = np.random.default_rng(20261013)
    cal = []
    for rho in (0.0, 0.3, 0.6, 1.0):
        for k in (3, 10, 30, 100):
            ifs, keffs = [], []
            for _ in range(150):
                n = 800
                common = rng.normal(size=n)
                Z = math.sqrt(rho) * common + math.sqrt(max(0.0, 1 - rho)) * rng.normal(size=(k, n))
                R = np.corrcoef(Z)
                bx = rng.normal(0.05, 0.01, size=k)
                by = 0.4 * bx + rng.normal(0, 0.02, size=k)
                se = rng.uniform(0.01, 0.02, size=k)
                nv = naive_ivw(bx, by, se)
                gl = gls_ivw(bx, by, se, R)
                if nv is None or gl is None:
                    continue
                ifs.append(gl[1] / nv[1])
                ratio = by / bx
                se_r = se / np.abs(bx)
                j = int(np.argmin(np.abs(ratio / se_r)))
                se_c = se_r[j]
                keffs.append((se_c / nv[1]) ** 2 / k)
            if ifs:
                cal.append({
                    "rho": rho, "k": k, "n_sim": len(ifs),
                    "IF_mean": round(float(np.mean(ifs)), 2),
                    "IF_theory_sqrt_1_plus_k1rho":
                        round(math.sqrt(1 + (k - 1) * rho), 2),
                    "keff_over_k_mean": round(float(np.mean(keffs)), 3),
                })
    if cal:
        with io.StringIO() as buf:
            w = csv.DictWriter(buf, fieldnames=list(cal[0].keys()), delimiter="\t",
                               lineterminator="\n")
            w.writeheader()
            w.writerows(cal)
            OUT_CAL.write_text(buf.getvalue(), encoding="utf-8", newline="")
        print("\n  IF   k_eff/k       ")
        print("  %-6s %-5s %8s %8s %12s" % ("rho", "k", "IF_sim", "IF_theory", "k_eff/k"))
        for c in cal:
            print("  %-6.1f %-5d %8.2f %8.2f %12.3f"
                  % (c["rho"], c["k"], c["IF_mean"],
                     c["IF_theory_sqrt_1_plus_k1rho"], c["keff_over_k_mean"]))

    si = rd(P_SENS, "\t")
    se_rows = []
    for b in ["2-3", "4-9", "10-99", ">=100"]:
        v = [r for r in si if r["n_snp_bin"] == b]
        if not v:
            continue
        se_rows.append({
            "n_snp_bin": b, "n_pairs": len(v),
            "median_se_unclumped": "%.4e" % statistics.median(float(r["se_raw"]) for r in v),
            "median_se_clumped": "%.4e" % statistics.median(float(r["se_clumped"]) for r in v),
            "min_se_unclumped": "%.4e" % min(float(r["se_raw"]) for r in v),
            "n_p_underflow": sum(1 for r in v if float(r["p_raw"]) == 0.0),
        })
    with io.StringIO() as buf:
        w = csv.DictWriter(buf, fieldnames=list(se_rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(se_rows)
        OUT_SE.write_text(buf.getvalue(), encoding="utf-8", newline="")
    print("\n     SE       SE %.2e         1e-308     "
          % min(float(r["min_se_unclumped"]) for r in se_rows))
    print("[04k]   ")


if __name__ == "__main__":
    main()
