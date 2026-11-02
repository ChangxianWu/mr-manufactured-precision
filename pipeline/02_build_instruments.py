# -*- coding: utf-8 -*-
"""Build exposure instruments from the MiBioGen top-hit table (P<5e-8; per-taxon variant lists; allele frequencies)."""
import argparse
import collections
import csv
import os
import re
import subprocess
import sys
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

_PROJ = "D:/Horzontal_project"
TOPHITS = _PROJ + "/20B3/work/data/mibiogen_tophits.tsv"
FIN_STUDIES = _PROJ + "/20B3/work/data/exposure_libraries/35115689_studies.tsv"
DMP_STUDIES = _PROJ + "/20B3/work/data/exposure_libraries/35115690_studies.tsv"
OUT_FUNNEL = _PROJ + "/20B3/work/output/iv_funnel.csv"
OUT_PRIMARY = _PROJ + "/20B3/work/output/instruments_primary.csv"
OUT_SENSITIVE = _PROJ + "/20B3/work/output/instruments_sensitive.csv"
OUT_TAXON = _PROJ + "/20B3/work/output/taxon_levels.csv"
PLINK = _PROJ + "/gwas_data/tools/plink.exe"
CLUMP_DIR = _PROJ + "/20B3/work/cache/clump"

_EUR_BFILES = {
    1: _PROJ + "/gwas_data/1000g_eur/eur_bfile/eur_chr1",
    2: _PROJ + "/gwas_data/1000g_eur/eur_bfile/eur_chr2",
    3: _PROJ + "/gwas_data/1000g_eur/eur_bfile/eur_chr3",
    4: _PROJ + "/gwas_data/1000g_eur/eur_bfile/eur_chr4",
    5: _PROJ + "/gwas_data/1000g_eur/eur_bfile/eur_chr5",
    6: _PROJ + "/gwas_data/1000g_eur/eur_bfile/eur_chr6",
    7: _PROJ + "/gwas_data/1000g_eur/eur_bfile/eur_chr7",
    8: _PROJ + "/gwas_data/1000g_eur/eur_bfile/eur_chr8",
    9: _PROJ + "/gwas_data/1000g_eur/eur_bfile/eur_chr9",
    10: _PROJ + "/gwas_data/1000g_eur/eur_bfile/eur_chr10",
    11: _PROJ + "/gwas_data/1000g_eur/eur_bfile/eur_chr11",
    12: _PROJ + "/gwas_data/1000g_eur/eur_bfile/eur_chr12",
    13: _PROJ + "/gwas_data/1000g_eur/eur_bfile/eur_chr13",
    14: _PROJ + "/gwas_data/1000g_eur/eur_bfile/eur_chr14",
    15: _PROJ + "/gwas_data/1000g_eur/eur_bfile/eur_chr15",
    16: _PROJ + "/gwas_data/1000g_eur/eur_bfile/eur_chr16",
    17: _PROJ + "/gwas_data/1000g_eur/eur_bfile/eur_chr17",
    18: _PROJ + "/gwas_data/1000g_eur/eur_bfile/eur_chr18",
    19: _PROJ + "/gwas_data/1000g_eur/eur_bfile/eur_chr19",
    20: _PROJ + "/gwas_data/1000g_eur/eur_bfile/eur_chr20",
    21: _PROJ + "/gwas_data/1000g_eur/eur_bfile/eur_chr21",
    22: _PROJ + "/gwas_data/1000g_eur/eur_bfile/eur_chr22",
}

PVAL_PRIMARY, PVAL_SENSITIVE = 5e-8, 1e-5
F_MIN, CLUMP_R2, CLUMP_KB = 10.0, 0.001, 10000
MHC = (6, 25_000_000, 35_000_000)
LCT = (2, 135_000_000, 138_000_000)
STRUCTURAL_LEVELS = ("phylum", "class", "order", "family", "genus")
_RSID_OK = re.compile(r"rs\d+\Z")


def parse_taxon(bac):
    """`genus.Bifidobacterium.id.436` -> ('genus', 'Bifidobacterium', 436) """
    m = re.match(r"^(%s)\.(.+?)\.id\.(\d+)\Z" % "|".join(STRUCTURAL_LEVELS), bac)
    if not m:
        return ("unknown", bac, None)
    return (m.group(1), m.group(2), int(m.group(3)))


def norm(t):
    return re.sub(r"[^a-z0-9]", "", t.lower())


def load_reference_names():
    """       Finnish Qin 2022    + DMP Lopera-Maya 2022     """
    fin_genera = set()
    if os.path.isfile(FIN_STUDIES):
        with open(FIN_STUDIES, encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                n = r["trait"].replace(" abundance in stool", "").strip()
                fin_genera.add(norm(n.split()[0]))
    dmp_names = set()
    if os.path.isfile(DMP_STUDIES):
        with open(DMP_STUDIES, encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                m = re.search(r"Gut microbiota abundance \((.+)\)", r["trait"])
                if m:
                    last = m.group(1).split(".")[-1]
                    dmp_names.add(norm(re.sub(r"^[a-z]__", "", last)))
    return fin_genera, dmp_names


def clump(rows, chrom):
    """plink --clump chrom     _EUR_BFILES      

              shell=False                    
    rsID   `rs\\d+`                    rsid    
              None         "  clump"       
    """
    if chrom not in _EUR_BFILES:
        return None
    bfile = _EUR_BFILES[chrom]
    if not os.path.isfile(bfile + ".bed") or not os.path.isfile(PLINK):
        return None
    safe = [(r["rsid"], r["P"]) for r in rows
            if r["rsid"] and _RSID_OK.match(r["rsid"])]
    if not safe:
        return None
    os.makedirs(CLUMP_DIR, exist_ok=True)
    inp = CLUMP_DIR + "/in_chr%d.txt" % chrom
    out_base = CLUMP_DIR + "/out_chr%d" % chrom
    with open(inp, "w", encoding="utf-8", newline="") as fh:
        fh.write("SNP\tP\n")
        for rs, p in safe:
            fh.write("%s\t%.6g\n" % (rs, p))
    argv = [PLINK, "--bfile", bfile, "--clump", inp,
            "--clump-p1", "1", "--clump-r2", str(CLUMP_R2),
            "--clump-kb", str(CLUMP_KB), "--out", out_base]
    proc = subprocess.run(argv, capture_output=True, text=True, shell=False,
                          check=False)
    clumped = out_base + ".clumped"
    if not os.path.isfile(clumped):
        sys.stderr.write("[02] plink clump    chr%d: %s\n"
                         % (chrom, (proc.stderr or proc.stdout or "")[-300:]))
        return None
    keep = set()
    with open(clumped, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("CHR\t") or not line.strip():
                continue
            c = line.split()
            if len(c) >= 3:
                keep.add(c[2])
    return keep


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-clump", action="store_true")
    a = ap.parse_args()

    with open(TOPHITS, encoding="utf-8") as fh:
        D = list(csv.DictReader(fh, delimiter="\t"))
    for r in D:
        r["rsid"] = (r.get("rsID") or "").strip()
        r["P"] = float(r["P"]); r["beta"] = float(r["beta"]); r["SE"] = float(r["SE"])
        r["chr"] = int(r["chr"]); r["pos"] = int(r["pos"]); r["N"] = float(r["N"] or 0)
        r["F"] = (r["beta"] / r["SE"]) ** 2
        r["z"] = float(r["metaZ_wz"]) if r["metaZ_wz"] else None
        r["level"], r["name"], r["tid"] = parse_taxon(r["bac"])

    taxa = collections.defaultdict(list)
    for r in D:
        taxa[r["bac"]].append(r)
    print("[02] tophits %d   / %d taxa API      P<=1e-4 " % (len(D), len(taxa)))

    bad = [r for r in D if r["z"] is not None and r["z"] * r["beta"] < 0]
    print("[02] beta/metaZ       %d   %.3f%% " % (len(bad), 100.0 * len(bad) / len(D)))

    fin_genera, dmp_names = load_reference_names()
    print("[02]       Finnish %d    / DMP %d    "
          % (len(fin_genera), len(dmp_names)))

    def in_lct(r):
        return r["chr"] == LCT[0] and LCT[1] <= r["pos"] <= LCT[2]

    def in_mhc(r):
        return r["chr"] == MHC[0] and MHC[1] <= r["pos"] <= MHC[2]

    funnel = [("L0_nominal_taxa", len(taxa), "MiBioGen    taxa  "),
              ("L0_rows", len(D), "tophits     P<=1e-4 ")]

    for thr, tag in ((PVAL_PRIMARY, "primary"), (PVAL_SENSITIVE, "sensitive")):
        elig = {k: [x for x in v if x["P"] < thr] for k, v in taxa.items()}
        elig = {k: v for k, v in elig.items() if v}
        funnel.append(("L1_%s_taxa_any_snp" % tag, len(elig),
                       "P<%g     1   SNP   taxa" % thr))
        funnel.append(("L1_%s_taxa_ge2snp" % tag,
                       sum(1 for v in elig.values() if len(v) >= 2),
                       "P<%g     2   SNP   taxa" % thr))

    prim = {k: v for k, v in taxa.items() if any(x["P"] < PVAL_PRIMARY for x in v)}
    lct_taxa = {k for k, v in prim.items() if in_lct(min(v, key=lambda x: x["P"]))}
    non_lct = {k: v for k, v in prim.items() if k not in lct_taxa}
    funnel.append(("L2_primary_taxa", len(prim), "P<5e-8      taxa"))
    funnel.append(("L2_lct_projected", len(lct_taxa),
                   "      LCT   chr2:135-138Mb   taxa"))
    funnel.append(("L2_non_lct_taxa", len(non_lct), "   LCT            taxa"))

    bins = collections.Counter()
    for v in non_lct.values():
        b = min(v, key=lambda x: x["P"])
        bins[(b["chr"], b["pos"] // 5_000_000)] += 1
    funnel.append(("L2_non_lct_loci_5mb", len(bins), "  LCT      5Mb     "))

    aligned = {}
    for k in prim:
        lvl, name, tid = parse_taxon(k)
        n = norm(name)
        aligned[k] = (n in fin_genera) or (n in dmp_names)
    funnel.append(("L3_aligned_taxa", sum(1 for x in aligned.values() if x),
                   "     Finnish   DMP   taxa"))
    funnel.append(("L3_non_lct_aligned",
                   sum(1 for k, x in aligned.items() if x and k in non_lct),
                   "  LCT       taxa"))

    with open(OUT_FUNNEL, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["layer", "n", "note"])
        w.writerows(funnel)
    for r in funnel:
        print("   %-28s %5s  %s" % (r[0], r[1], r[2]))

    with open(OUT_TAXON, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["bac", "level", "name", "n_tophits", "best_rsid", "chr", "pos_grch37",
                    "best_p", "best_beta", "best_se", "F", "N", "Ncohorts",
                    "n_sig_snp", "in_lct", "in_mhc", "aligned_fin_or_dmp"])
        for k, v in sorted(taxa.items(), key=lambda kv: min(x["P"] for x in kv[1])):
            b = min(v, key=lambda x: x["P"])
            lvl, nm, _ = parse_taxon(k)
            w.writerow([k, lvl, nm, len(v), b["rsid"], b["chr"], b["pos"],
                        "%.4g" % b["P"], b["beta"], b["SE"], "%.2f" % b["F"],
                        int(b["N"]), b["Ncohorts"],
                        sum(1 for x in v if x["P"] < PVAL_PRIMARY),
                        int(in_lct(b)), int(in_mhc(b)), int(aligned.get(k, False))])

    def build(thr, tag, out_path, do_clump):
        recs = []
        for k, v in taxa.items():
            sig = sorted([x for x in v if x["P"] < thr and not in_mhc(x)],
                         key=lambda x: x["P"])
            if not sig:
                continue
            if do_clump:
                by_chr = collections.defaultdict(list)
                for x in sig:
                    by_chr[x["chr"]].append(x)
                kept = set()
                for c, rs in by_chr.items():
                    kp = clump(rs, c)
                    if kp is None:
                        sys.stderr.write("[02] chr%s             clump\n" % c)
                        kp = {x["rsid"] for x in rs}
                    kept |= kp
                sig = [x for x in sig if x["rsid"] in kept]
            if not sig:
                continue
            Fs = [x["F"] for x in sig]
            if len(sig) == 1 and Fs[0] < F_MIN:
                continue
            lvl, nm, _ = parse_taxon(k)
            b = sig[0]
            recs.append({
                "taxon": k, "level": lvl, "name": nm,
                "n_snp": len(sig), "rsid": b["rsid"],
                "chr": b["chr"], "pos_grch37": b["pos"],
                "ea": b["alt"].upper(), "oa": b["ref"].upper(),
                "beta": b["beta"], "se": b["SE"], "pval": b["P"],
                "F_mean": sum(Fs) / len(Fs), "F_min": min(Fs), "F_max": max(Fs),
                "N": int(b["N"]), "Ncohorts": b["Ncohorts"],
                "in_lct": int(in_lct(b)), "aligned": int(aligned.get(k, False)),
                "snp_list": ";".join("%s:%.6g:%.6g:%.6g:%s:%s"
                                     % (x["rsid"], x["beta"], x["SE"], x["P"],
                                        x["alt"].upper(), x["ref"].upper())
                                     for x in sig),
            })
        with open(out_path, "w", encoding="utf-8", newline="") as fh:
            if recs:
                w = csv.DictWriter(fh, fieldnames=list(recs[0].keys()), delimiter="\t")
                w.writeheader(); w.writerows(recs)
        print("[02] %s %d   taxa       clump=%s " % (tag, len(recs), do_clump))
        return recs

    prim_recs = build(PVAL_PRIMARY, "primary", OUT_PRIMARY, not a.no_clump)
    build(PVAL_SENSITIVE, "sensitive", OUT_SENSITIVE, False)

    print("[02]       taxa %d    LCT %d "
          % (len(prim_recs), sum(r["in_lct"] for r in prim_recs)))
    if prim_recs:
        print("[02]      SNP    %d   SNP %d"
              % (sum(1 for r in prim_recs if r["n_snp"] > 1),
                 sum(1 for r in prim_recs if r["n_snp"] == 1)))


if __name__ == "__main__":
    main()
