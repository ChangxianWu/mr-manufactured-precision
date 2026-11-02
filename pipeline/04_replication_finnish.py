# -*- coding: utf-8 -*-
"""Cross-cohort replication against an independent Finnish single-cohort mGWAS."""
import collections
import csv
import gzip
import json
import os
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

BASE = os.path.realpath("D:/Horzontal_project")
FIN_DIR = BASE + "/20B3/work/data/finnish"
PLAN = BASE + "/20B3/work/data/exposure_libraries/finnish_replication_plan.tsv"
INSTR = BASE + "/20B3/work/output/instruments_primary.csv"
TOPHITS = BASE + "/20B3/work/data/mibiogen_tophits.tsv"
MAP = BASE + "/20B3/work/data/lead_snp_hg38_map.csv"
OUTDIR = BASE + "/20B3/work/output"

FIN_FILES = {}


def _fill_fin_files():
    for acc in ("GCST90032179", "GCST90032219", "GCST90032220", "GCST90032221",
                "GCST90032222", "GCST90032223", "GCST90032224", "GCST90032225",
                "GCST90032226", "GCST90032227", "GCST90032228", "GCST90032229",
                "GCST90032230", "GCST90032379", "GCST90032380", "GCST90032539",
                "GCST90032568"):
        FIN_FILES[acc] = FIN_DIR + "/" + acc + "_buildGRCh38.tsv.gz"


_fill_fin_files()

FIN_HEADER = ["chromosome", "base_pair_location", "effect_allele", "other_allele",
              "effect_allele_frequency", "beta", "standard_error", "p_value"]
COMPL = {"A": "T", "T": "A", "C": "G", "G": "C"}


def load_instruments():
    """   SNP taxon → [{rsid, chr37, pos37, ea, oa, beta, se}] 

    ⚠️      **   SNP   **  GRCh37       MiBioGen tophits  
            taxon       ——       taxon     SNP    
               283   SNP    1        
    """
    rs2p = {}
    with open(TOPHITS, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rs = (r.get("rsID") or "").strip()
            if rs:
                rs2p[(r["bac"], rs)] = (int(r["chr"]), int(r["pos"]))
    out = collections.defaultdict(list)
    miss = 0
    with open(INSTR, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            for item in (r["snp_list"] or "").split(";"):
                p = item.split(":")
                if len(p) != 6:
                    continue
                c = rs2p.get((r["taxon"], p[0]))
                if c is None:
                    miss += 1
                    continue
                out[r["taxon"]].append({
                    "rsid": p[0], "chr37": c[0], "pos37": c[1],
                    "ea": p[4].upper(), "oa": p[5].upper(),
                    "beta": float(p[1]), "se": float(p[2]),
                })
    if miss:
        print("[04] ⚠️ %d   SNP    SNP        " % miss)
    return out


def load_pos38():
    """(chr37, pos37) → (chr38, pos38) 

          `rsid / chrom / pos38`    03a   rsID   Ensembl      
        MiBioGen tophits   rsID      GRCh37    
    """
    rs2p37 = {}
    with open(TOPHITS, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rs = (r.get("rsID") or "").strip()
            if rs:
                rs2p37[rs] = (int(r["chr"]), int(r["pos"]))
    m = {}
    with open(MAP, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            k = rs2p37.get(r["rsid"])
            if k and r["pos38"]:
                m[k] = (int(r["chrom"]), int(r["pos38"]))
    return m


def load_plan():
    """taxon → [(finnish_trait, accession)] """
    out = collections.defaultdict(list)
    with open(PLAN, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["mibiogen_taxon"]].append((r["finnish_trait"], r["finnish_accession"]))
    return out


def streams_for(accessions):
    """  accession        (chr38, pos38)    """
    return accessions


def scan_file(acc, want_pos):
    """       Finnish sumstats    {pos38: row}     want_pos       """
    fp = FIN_FILES[acc]
    if not os.path.isfile(fp):
        return None
    hit = {}
    with gzip.open(fp, "rt", errors="replace") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ci = {c: i for i, c in enumerate(hdr)}
        c_chr = ci.get("chromosome", 0)
        c_pos = ci.get("base_pair_location", 1)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(hdr):
                continue
            key = f[c_pos]
            if key not in want_pos:
                continue
            hit[key] = f
    return hit


def main():
    ins = load_instruments()
    pos38 = load_pos38()
    plan = load_plan()

    need = collections.defaultdict(set)     # acc -> {pos str}
    key2taxa = collections.defaultdict(list)  # (acc, pos) -> [(taxon, rsid, ea, oa, beta, se)]
    for taxon, snps in ins.items():
        for tr, acc in plan.get(taxon, ()):
            for s in snps:
                k = (s["chr37"], s["pos37"])
                if k not in pos38:
                    continue
                c38, p38 = pos38[k]
                need[acc].add(str(p38))
                key2taxa[(acc, str(p38))].append((taxon, tr, s))
    print("[04]         taxa %d Finnish accession %d        %d"
          % (len(plan), len(need), sum(len(v) for v in need.values())))

    rows = []
    for acc in sorted(need):
        hit = scan_file(acc, need[acc])
        if hit is None:
            print("  [%s]        " % acc)
            continue
        n = 0
        for p38, f in hit.items():
            d = dict(zip(FIN_HEADER, f))
            try:
                fb = float(d["beta"]); fs = float(d["standard_error"])
                fe = d["effect_allele"].upper(); fo = d["other_allele"].upper()
                fef = float(d["effect_allele_frequency"])
            except Exception:  # noqa: BLE001
                continue
            for taxon, tr, s in key2taxa[(acc, p38)]:
                if fe == s["ea"] and fo == s["oa"]:
                    bi, st = fb, "same"
                elif fe == s["oa"] and fo == s["ea"]:
                    bi, st = -fb, "flip"
                else:
                    continue
                if fef <= 0 or fef >= 1 or fs <= 0:
                    continue
                rows.append({
                    "taxon": taxon, "finnish_trait": tr, "accession": acc,
                    "rsid": s["rsid"], "pos_grch38": p38,
                    "ea_mibiogen": s["ea"], "oa_mibiogen": s["oa"],
                    "beta_mibiogen": s["beta"], "se_mibiogen": s["se"],
                    "ea_finnish": fe, "oa_finnish": fo, "eaf_finnish": fef,
                    "beta_finnish": bi, "se_finnish": fs, "align": st,
                    "same_sign": int((s["beta"] > 0) == (bi > 0)),
                })
                n += 1
        print("  [%s]      %d    %d " % (acc, n, len(need[acc])))

    with open(OUTDIR + "/replication_exposure.csv", "w",
              encoding="utf-8", newline="") as fh:
        if rows:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
            w.writeheader(); w.writerows(rows)
    print("[04]        %d   -> replication_exposure.csv" % len(rows))

    by = collections.defaultdict(list)
    for r in rows:
        by[r["taxon"]].append(r)
    summ = []
    for taxon, rs in sorted(by.items()):
        n = len(rs)
        same = sum(r["same_sign"] for r in rs)
        import statistics as st
        corr = ""
        if n >= 3:
            xs = [r["beta_mibiogen"] for r in rs]
            ys = [r["beta_finnish"] for r in rs]
            mx, my = sum(xs) / n, sum(ys) / n
            sx = (sum((v - mx) ** 2 for v in xs)) ** 0.5
            sy = (sum((v - my) ** 2 for v in ys)) ** 0.5
            if sx > 0 and sy > 0:
                corr = "%.3f" % (sum((xs[i] - mx) * (ys[i] - my) for i in range(n)) / (sx * sy))
        summ.append({
            "taxon": taxon, "n_snp": n,
            "n_finnish_traits": len({r["finnish_trait"] for r in rs}),
            "same_sign": same, "same_sign_rate": "%.3f" % (same / n) if n else "",
            "beta_corr_pearson": corr,
        })
    with open(OUTDIR + "/replication_summary.csv", "w",
              encoding="utf-8", newline="") as fh:
        if summ:
            w = csv.DictWriter(fh, fieldnames=list(summ[0].keys()), delimiter="\t")
            w.writeheader(); w.writerows(summ)
    print("[04]   taxon    %d   -> replication_summary.csv" % len(summ))
    for s in summ:
        print("   %-40s n=%-4s same_sign=%s %s corr=%s"
              % (s["taxon"], s["n_snp"], s["same_sign"], s["same_sign_rate"],
                 s["beta_corr_pearson"]))


if __name__ == "__main__":
    main()
