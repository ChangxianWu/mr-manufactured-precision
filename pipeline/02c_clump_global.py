# ============================================================
# ------------------------------------------------------------------
#
# ============================================================
import csv
import sys
from pathlib import Path

ROOT = "D:/Horzontal_project/20B3"
TOPHITS = ROOT + "/work/data/mibiogen_tophits.tsv"
GENO_CACHE = ROOT + "/work/cache/eur_geno_cache.tsv"
OUTDIR = ROOT + "/work/output"
PVAL = 5e-8
KB = 10000
R2_SET = [0.001, 0.01, 0.1]


def r2_of(x, y):
    n = 0
    sx = sy = sxx = syy = sxy = 0
    for a, b in zip(x, y):
        if a is None or b is None:
            continue
        n += 1
        sx += a
        sy += b
        sxx += a * a
        syy += b * b
        sxy += a * b
    if n < 10:
        return None
    mx, my = sx / n, sy / n
    vx = sxx / n - mx * mx
    vy = syy / n - my * my
    if vx <= 0 or vy <= 0:
        return None
    return ((sxy / n - mx * my) ** 2) / (vx * vy)


def main():
    pool = {}
    with open(TOPHITS, encoding="utf-8") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            try:
                p = float(r["P"])
            except (TypeError, ValueError):
                continue
            if p >= PVAL:
                continue
            key = (int(r["chr"]), r["pos"])
            e = pool.setdefault(key, {"rsid": r["rsID"], "p": p, "taxa": set()})
            e["taxa"].add(r["bac"])
            if p < e["p"]:
                e["p"] = p
                e["rsid"] = r["rsID"]
    print("[02c] P<%.0e       (chr,pos) = %d    taxa = %d"
          % (PVAL, len(pool), len({t for v in pool.values() for t in v["taxa"]})))

    geno = {}
    with open(GENO_CACHE, encoding="utf-8") as fh:
        fh.readline()
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            geno[(int(f[0]), f[1])] = [None if x == "." else int(x)
                                       for x in f[2].split(",")]
    have = [k for k in pool if k in geno]
    print("[02c]      %d / %d   " % (len(have), len(pool)))

    order = sorted(have, key=lambda k: (pool[k]["p"], k))
    res = {}
    per_kept = []
    for r2_max in R2_SET:
        kept = []
        dropped = 0
        for k in order:
            chr_, pos = k
            bad = None
            for kk in kept:
                if kk[0] != chr_:
                    continue
                if abs(int(kk[1]) - int(pos)) > KB * 1000:
                    continue
                rr = r2_of(geno[kk], geno[k])
                if rr is not None and rr > r2_max:
                    bad = (kk, rr)
                    break
            if bad is None:
                kept.append(k)
            else:
                dropped += 1
                per_kept.append({"r2_max": r2_max, "kept": pool[k]["rsid"],
                                 "chr": chr_, "pos": pos, "p": "%.4g" % pool[k]["p"],
                                 "dropped_by": pool[bad[0]]["rsid"],
                                 "r2": "%.4f" % bad[1],
                                 "n_taxa": len(pool[k]["taxa"])})
        res[r2_max] = kept
        print("[02c]    r2<=%-5s %d -> %d    %.1f× "
              % (r2_max, len(order), len(kept),
                 (len(order) / len(kept)) if kept else float("nan")))

    r2m = R2_SET[0]
    kept = res[r2m]
    rows = []
    for k in kept:
        e = pool[k]
        rows.append({"rsid": e["rsid"], "chr": k[0], "pos_grch37": k[1],
                     "p_min": "%.4g" % e["p"], "n_taxa": len(e["taxa"]),
                     "taxa": ";".join(sorted(e["taxa"]))})
    rows.sort(key=lambda r: (int(r["chr"]), int(r["pos_grch37"])))
    with open(OUTDIR + "/clump_global.csv", "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["rsid", "chr", "pos_grch37", "p_min",
                                           "n_taxa", "taxa"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    summ = []
    for r2_max in R2_SET:
        summ.append({"r2_max": r2_max, "n_variants_pooled": len(order),
                     "n_independent": len(res[r2_max]),
                     "inflation_factor": "%.1f" % (len(order) / len(res[r2_max]))})
    with open(OUTDIR + "/clump_global_summary.csv", "w",
              encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["r2_max", "n_variants_pooled",
                                           "n_independent", "inflation_factor"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(summ)

    with open(OUTDIR + "/clump_global_pairs.csv", "w",
              encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["r2_max", "kept", "chr", "pos", "p",
                                           "dropped_by", "r2", "n_taxa"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(per_kept)

    cov = {t for k in kept for t in pool[k]["taxa"]}
    print("[02c]      %d      taxa %d      clump_global*.csv"
          % (len(kept), len(cov)))


if __name__ == "__main__":
    sys.exit(main())
