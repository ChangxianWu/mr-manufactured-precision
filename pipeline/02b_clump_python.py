# -*- coding: utf-8 -*-
"""Pure-Python LD clumping against 1000 Genomes EUR (greedy, by ascending P, r2 threshold within a 10 Mb window)."""
import argparse
import collections
import csv
import gzip
import os
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

_D = "D:/Horzontal_project"
KEEP = _D + "/20B3/work/data/1000g_eur_keep2.txt"
TOPHITS = _D + "/20B3/work/data/mibiogen_tophits.tsv"
INSTR = _D + "/20B3/work/output/instruments_primary.csv"
OUTDIR = _D + "/20B3/work/output"
GENO_CACHE = _D + "/20B3/work/cache/eur_geno_cache.tsv"

VCF = {
    1: _D + "/gwas_data/1000g_eur/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    2: _D + "/gwas_data/1000g_eur/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    3: _D + "/gwas_data/1000g_eur/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    4: _D + "/gwas_data/1000g_eur/ALL.chr4.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    5: _D + "/gwas_data/1000g_eur/ALL.chr5.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    6: _D + "/gwas_data/1000g_eur/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    7: _D + "/gwas_data/1000g_eur/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    8: _D + "/gwas_data/1000g_eur/ALL.chr8.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    9: _D + "/gwas_data/1000g_eur/ALL.chr9.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    10: _D + "/gwas_data/1000g_eur/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    11: _D + "/gwas_data/1000g_eur/ALL.chr11.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    12: _D + "/gwas_data/1000g_eur/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    13: _D + "/gwas_data/1000g_eur/ALL.chr13.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    14: _D + "/gwas_data/1000g_eur/ALL.chr14.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    15: _D + "/gwas_data/1000g_eur/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    16: _D + "/gwas_data/1000g_eur/ALL.chr16.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    17: _D + "/gwas_data/1000g_eur/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    18: _D + "/gwas_data/1000g_eur/ALL.chr18.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    19: _D + "/gwas_data/1000g_eur/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    20: _D + "/gwas_data/1000g_eur/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    21: _D + "/gwas_data/1000g_eur/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    22: _D + "/gwas_data/1000g_eur/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
}

CLUMP_PAIRS = {
    0.001: _D + "/20B3/work/output/clump_pairs_0p001.csv",
    0.01: _D + "/20B3/work/output/clump_pairs_0p01.csv",
    0.1: _D + "/20B3/work/output/clump_pairs_0p1.csv",
}
CLUMP_SUMMARY = _D + "/20B3/work/output/clump_summary.csv"
INSTR_CLUMPED = _D + "/20B3/work/output/instruments_clumped.csv"


def load_keep():
    keep = []
    with open(KEEP, encoding="utf-8") as fh:
        for line in fh:
            c = line.split()
            if c:
                keep.append(c[0])
    return keep


def load_snp_coords():
    m = {}
    with open(TOPHITS, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rs = (r.get("rsID") or "").strip()
            if rs:
                m[(r["bac"], rs)] = (int(r["chr"]), int(r["pos"]))
    return m


def load_candidates(coords):
    """  taxon       SNP {taxon: [(rsid, chrom, pos, P)]} """
    out = collections.defaultdict(list)
    with open(INSTR, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            for item in (r["snp_list"] or "").split(";"):
                p = item.split(":")
                if len(p) != 6:
                    continue
                c = coords.get((r["taxon"], p[0]))
                if c is None:
                    continue
                out[r["taxon"]].append((p[0], c[0], c[1], float(p[3])))
    return out


def sample_columns(chrom):
    with gzip.open(VCF[chrom], "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]
    return []


def read_genotypes(chrom, positions, keep_idx):
    """  positions      EUR     0/1/2       None  

       ({pos_str: doses},     ) POS    **   **   
    """
    want = {str(p) for p in positions}
    got = {}
    n_seen = 0
    with gzip.open(VCF[chrom], "rt", errors="replace") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            head = line.split("\t", 2)
            if len(head) < 3:
                continue
            n_seen += 1
            key = head[1]
            if key not in want:
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            fmt = f[8].split(":")
            gi = fmt.index("GT") if "GT" in fmt else 0
            doses = []
            for idx in keep_idx:
                sample = f[9 + idx]
                gt = sample.split(":", gi + 1)[gi] if ":" in sample else sample
                if gt in (".", "./.", ".|."):
                    doses.append(None)
                    continue
                a = gt.replace("|", "/").split("/")
                if len(a) != 2 or "." in a:
                    doses.append(None)
                    continue
                doses.append(int(a[0]) + int(a[1]))
            got[key] = doses
    return got, n_seen


def build_genotypes(by_chr, keep, append=False):
    """   by_chr          

    `append=False`                   
    `append=True`        **  **             
    """
    allg = {}
    for c in sorted(by_chr):
        names = sample_columns(c)
        n2i = {n: i for i, n in enumerate(names)}
        keep_idx = [n2i[s] for s in keep if s in n2i]
        got, n_seen = read_genotypes(c, by_chr[c], keep_idx)
        allg[c] = got
        print("     chr%-2d    %d      %d / %d    | EUR %d/%d"
              % (c, n_seen, len(got), len(by_chr[c]), len(keep_idx), len(keep)),
              flush=True)
    mode = "a" if append and os.path.isfile(GENO_CACHE) else "w"
    with open(GENO_CACHE, mode, encoding="utf-8", newline="") as fh:
        if mode == "w":
            fh.write("chrom\tpos\tdoses\n")
        for c in sorted(allg):
            for p, doses in sorted(allg[c].items(), key=lambda kv: int(kv[0])):
                fh.write("%d\t%s\t%s\n" % (
                    c, p, ",".join("." if v is None else str(v) for v in doses)))
    print("[02b]      %s %s    %d    "
          % ("  " if mode == "a" else "  ", GENO_CACHE,
             sum(len(v) for v in allg.values())))
    return allg


def load_genotypes(by_chr, keep):
    """    **        **            

        12                              
                                     
    """
    out = collections.defaultdict(dict)
    if os.path.isfile(GENO_CACHE):
        with open(GENO_CACHE, encoding="utf-8") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            if hdr[:2] == ["chrom", "pos"]:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) < 3:
                        continue
                    out[int(f[0])][f[1]] = [None if x == "." else int(x)
                                            for x in f[2].split(",")]
    have = {c for c, v in out.items() if v}
    missing = {c: v for c, v in by_chr.items() if c not in have}
    if not missing:
        print("[02b]         %d     " % len(have))
        return out
    print("[02b]     %s     " % sorted(missing))
    new = build_genotypes(missing, keep, append=True)
    for c, v in new.items():
        out[c] = v
    return out


def r2_of(x, y):
    """      0/1/2   Pearson r²         """
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


def clump_taxon(cands, geno, r2_max, kb):
    """   clumping   P         SNP   r² <= r2_max   

    `geno`    **   **   VCF      ——   int        SNP   
    "no_genotype" clumping      0       
    """
    kept = []
    dropped = {}
    for rsid, _c, p, _pv in sorted(cands, key=lambda t: t[3]):
        g = geno.get(str(p))
        if g is None:
            dropped[rsid] = "no_genotype"
            continue
        bad = None
        for k_rs, k_p in kept:
            if abs(p - k_p) > kb * 1000:
                continue
            r2 = r2_of(g, geno.get(str(k_p)) or [])
            if r2 is not None and r2 > r2_max:
                bad = "r2=%.3f with %s" % (r2, k_rs)
                break
        if bad:
            dropped[rsid] = bad
        else:
            kept.append((rsid, p))
    return kept, dropped


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--r2", default="0.001,0.01,0.1")
    ap.add_argument("--kb", type=int, default=10000)
    a = ap.parse_args()
    thresholds = [float(x) for x in a.r2.split(",") if x]

    keep = load_keep()
    coords = load_snp_coords()
    cands = load_candidates(coords)
    by_chr = collections.defaultdict(set)
    for _t, lst in cands.items():
        for _rs, c, p, _pv in lst:
            by_chr[c].add(p)
    print("[02b]    SNP %d taxa     %d        %d"
          % (len(cands), len(by_chr), sum(len(v) for v in by_chr.values())))

    allg = load_genotypes(by_chr, keep)

    summary = []
    for r2_max in thresholds:
        per_taxon = []
        for t, lst in cands.items():
            genome = {str(p): allg[c].get(str(p)) for _r, c, p, _v in lst}
            kept, dropped = clump_taxon(lst, genome, r2_max, a.kb)
            per_taxon.append({
                "taxon": t, "n_before": len(lst), "n_after": len(kept),
                "kept": ";".join(k for k, _ in kept),
                "dropped": ";".join("%s(%s)" % (k, v) for k, v in dropped.items()),
            })
        tot_b = sum(r["n_before"] for r in per_taxon)
        tot_a = sum(r["n_after"] for r in per_taxon)
        summary.append({"r2_max": r2_max, "n_taxa": len(per_taxon),
                        "n_snp_before": tot_b, "n_snp_after": tot_a,
                        "inflation_factor": "%.1f" % (tot_b / tot_a) if tot_a else ""})
        with open(CLUMP_PAIRS[r2_max], "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(per_taxon[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(per_taxon)
        print("[02b] r2<=%-5s SNP %d -> %d    %.1f× "
              % (r2_max, tot_b, tot_a, (tot_b / tot_a) if tot_a else 0))

    with open(CLUMP_SUMMARY, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(summary)

    r2m = thresholds[0]
    kept_map = {}
    with open(CLUMP_PAIRS[r2m], encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            kept_map[r["taxon"]] = set((r["kept"] or "").split(";")) - {""}
    rows = []
    with open(INSTR, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            keepset = kept_map.get(r["taxon"], set())
            items = [x for x in (r["snp_list"] or "").split(";")
                     if x and x.split(":")[0] in keepset]
            if not items:
                continue
            d = dict(r)
            d["n_snp_pre_clump"] = r["n_snp"]
            d["n_snp"] = len(items)
            d["snp_list"] = ";".join(items)
            best = min(items, key=lambda x: float(x.split(":")[3])).split(":")
            d["rsid"] = best[0]
            d["beta"] = best[1]
            d["se"] = best[2]
            d["pval"] = best[3]
            d["ea"] = best[4]
            d["oa"] = best[5]
            rows.append(d)
    if not rows:
        raise SystemExit("[02b] clumping      SNP ——           ")
    with open(INSTR_CLUMPED, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print("[02b] clumped     %d taxa    27  SNP    %d"
          % (len(rows), sum(int(r["n_snp"]) for r in rows)))
    for r in sorted(rows, key=lambda x: -int(x["n_snp"])):
        print("   %-42s %4s -> %3s" % (r["taxon"], r["n_snp_pre_clump"], r["n_snp"]))


if __name__ == "__main__":
    main()
