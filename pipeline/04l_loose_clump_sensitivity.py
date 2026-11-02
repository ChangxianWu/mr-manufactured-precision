# -*- coding: utf-8 -*-
"""Relaxed clumping thresholds (r2 0.5 / 0.8): survivors by axis, Supplementary S12."""
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
D_EXTRACT = "D:/Horzontal_project/20B3/work/data/outcome_extract"
OUT_SUM = Path("D:/Horzontal_project/20B3/work/output/loose_clump_summary.csv")
OUT_TESTS = Path("D:/Horzontal_project/20B3/work/output/loose_clump_tests.csv")

THRESHOLDS = [0.5, 0.8]
KB = 10000
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


def align(ea, oa, ref, alt, beta_out, eaf_out):
    """  03b/03c           MAF   0.5             """
    if not all([ea, oa, ref, alt]):
        return None
    e, o, r, a = ea.upper(), oa.upper(), ref.upper(), alt.upper()
    if is_pal(e, o) and eaf_out is not None and abs(eaf_out - 0.5) < 0.05:
        return None
    if e == a and o == r:
        return beta_out
    if e == r and o == a:
        return -beta_out
    return None


def load_geno():
    geno = {}
    for ln in P_GENO.read_text(encoding="utf-8").splitlines()[1:]:
        f = ln.split("\t")
        if len(f) < 3:
            continue
        v = np.array([np.nan if x == "." else float(x) for x in f[2].split(",")])
        if np.isnan(v).any():
            v = np.where(np.isnan(v), np.nanmean(v), v)
        geno[(int(f[0]), str(f[1]))] = v
    return geno


def r2(a, b):
    if a.var() <= 0 or b.var() <= 0:
        return 0.0
    v = float(np.corrcoef(a, b)[0, 1])
    return 0.0 if not math.isfinite(v) else v * v


def main():
    prim = rd(P_INSTR, "\t")
    N_OF = {}
    for _r in prim:
        try:
            N_OF[_r["taxon"]] = float(_r["N"])
        except (KeyError, TypeError, ValueError):
            continue
    top = rd(P_TOPHITS, "\t")
    pos_of = {}
    for r in top:
        try:
            if float(r["P"]) < 5e-8:
                pos_of[r["rsID"]] = (int(r["chr"]), str(r["pos"]))
        except (TypeError, ValueError):
            continue
    geno = load_geno()

    # taxon -> [(rsid, chr, pos, p, beta, se)]
    cand = {}
    for r in prim:
        lst = []
        for tok in (r["snp_list"] or "").split(";"):
            tok = tok.strip()
            if not tok:
                continue
            f = tok.split(":")
            if len(f) < 6:
                continue
            key = pos_of.get(f[0])
            if key is None or key not in geno:
                continue
            try:
                lst.append((f[0], key[0], int(key[1]), float(f[3]), float(f[1]), float(f[2])))
            except ValueError:
                continue
        cand[r["taxon"]] = {"name": r["name"], "level": r["level"],
                            "in_lct": r["in_lct"], "list": lst}

    out = {}
    for ep in ENDPOINTS:
        f = "%s/outcome_extract_%s.tsv" % (D_EXTRACT, ep)
        if not Path(f).exists():
            continue
        d = {}
        for r in rd(f, "\t"):
            try:
                bo, so = float(r["beta_out"]), float(r["se_out"])
            except (TypeError, ValueError):
                continue
            if so <= 0:
                continue
            try:
                eaf = float(r["eaf_out"])
            except (TypeError, ValueError):
                eaf = None
            d[(r["taxon"], r["rsid"])] = (bo, so, r["ref_out"], r["alt_out"],
                                          r["ea_instr"], r["oa_instr"], eaf)
        out[ep] = d

    tests = []
    summary = []
    for r2max in THRESHOLDS:
        n_slot = 0
        per_taxon_kept = {}
        for t, meta in cand.items():
            lst = sorted(meta["list"], key=lambda x: x[3])
            kept = []
            for c in lst:
                bad = False
                for k in kept:
                    if k[1] != c[1] or abs(k[2] - c[2]) > KB * 1000:
                        continue
                    if r2(geno[(k[1], str(k[2]))], geno[(c[1], str(c[2]))]) > r2max:
                        bad = True
                        break
                if not bad:
                    kept.append(c)
            per_taxon_kept[t] = kept
            n_slot += len(kept)

        sig_pairs = []
        for ep, d in out.items():
            for t, kept in per_taxon_kept.items():
                used = []
                for rs in kept:
                    if (t, rs[0]) not in d:
                        continue
                    bo, so, ref, alt, eai, oai, eaf = d[(t, rs[0])]
                    ba = align(eai, oai, ref, alt, bo, eaf)
                    if ba is None:
                        continue
                    used.append((rs[0], rs[4], rs[5], ba, so))
                if not used:
                    continue
                ratios, ses = [], []
                r2_exp, used_rs, n_exp = 0.0, [], None
                bx1 = sx1 = eaf1 = None
                for rs, bx, sx, by, sy in used:
                    _eaf = d[(t, rs)][6]
                    if _eaf is not None and 0 < _eaf < 1:
                        r2_exp += 2.0 * _eaf * (1.0 - _eaf) * bx * bx
                    used_rs.append(rs)
                    n_exp = N_OF.get(t, n_exp)
                    if bx1 is None:
                        bx1, sx1, eaf1 = bx, sx, _eaf
                    if bx == 0:
                        continue
                    ratios.append(by / bx)
                    ses.append(sy / abs(bx))
                if not ratios:
                    continue
                w = np.array([1.0 / (s * s) for s in ses])
                r = np.array(ratios)
                b = float((w * r).sum() / w.sum())
                se = float(math.sqrt(1.0 / w.sum()))
                p = 2.0 * (1.0 - 0.5 * (1.0 + math.erf(abs(b / se) / math.sqrt(2.0))))
                tests.append({"r2_max": r2max, "endpoint": ep, "taxon": t,
                              "name": cand[t]["name"], "level": cand[t]["level"],
                              "in_lct": cand[t]["in_lct"], "n_snp": len(ratios),
                              "b": b, "se": se, "p": p,
                              "rsid": ";".join(used_rs),
                              "r2_exp_sum": r2_exp, "N_exp": n_exp,
                              "beta_exp1": bx1, "se_exp1": sx1, "eaf1": eaf1})
                sig_pairs.append((p, ep, t, cand[t]["name"], len(ratios)))

        # BH
        ps = sorted(x["p"] for x in tests if x["r2_max"] == r2max)
        m = len(ps)
        bh = 0
        for i, p in enumerate(ps, 1):
            if p <= 0.05 * i / m:
                bh = i
        survivors = []
        if bh:
            thr = sorted(ps)[bh - 1]
            survivors = [x for x in tests if x["r2_max"] == r2max and x["p"] <= thr]
        summary.append({
            "r2_max": r2max,
            "n_taxa_with_instruments": sum(1 for v in per_taxon_kept.values() if v),
            "n_slots": n_slot,
            "n_tests": m,
            "bh_survivors": bh,
            "min_p": "%.3e" % min(ps),
            "max_n_snp": max(len(v) for v in per_taxon_kept.values()),
            "survivor_names": ";".join(sorted({s["name"] for s in survivors})),
        })
        print("[04l] r2<=%.1f: taxa %d, slots %d, tests %d, BH %d, min p %.3e, max k %d"
              % (r2max, sum(1 for v in per_taxon_kept.values() if v), n_slot, m, bh,
                 min(ps), max(len(v) for v in per_taxon_kept.values())))
        for s in sorted(survivors, key=lambda x: x["p"])[:6]:
            print("        %-26s %-12s k=%d p=%.2e" % (s["name"][:26], s["endpoint"],
                                                       s["n_snp"], s["p"]))

    with io.StringIO() as buf:
        w = csv.DictWriter(buf, fieldnames=list(summary[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(summary)
        OUT_SUM.write_text(buf.getvalue(), encoding="utf-8", newline="")

    with io.StringIO() as buf:
        w = csv.DictWriter(buf, fieldnames=list(tests[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(tests)
        OUT_TESTS.write_text(buf.getvalue(), encoding="utf-8", newline="")
    print("[04l]   ")


if __name__ == "__main__":
    main()
