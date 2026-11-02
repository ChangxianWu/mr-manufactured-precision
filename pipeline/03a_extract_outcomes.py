# -*- coding: utf-8 -*-
"""Remote tabix extraction of FinnGen R12 outcome summary statistics at the instrument windows."""
import argparse
import collections
import csv
import json
import os
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

BASE = os.path.realpath("D:/Horzontal_project")
sys.path.insert(0, BASE)
from cv_shared.bgzf_tabix import TabixReader  # noqa: E402
from cv_shared.net_guard import open_url  # noqa: E402

INSTR = BASE + "/20B3/work/output/instruments_primary.csv"
INSTR_SENS = BASE + "/20B3/work/output/instruments_sensitive.csv"
INSTR_CLUMPED = BASE + "/20B3/work/output/instruments_clumped.csv"
OUTDIR_PRIMARY = BASE + "/20B3/work/data/outcome_extract"
OUTDIR_SENSITIVE = BASE + "/20B3/work/data/outcome_extract_sensitive"
CACHE_PRIMARY = BASE + "/20B3/work/cache/tbi"
CACHE_SENSITIVE = BASE + "/20B3/work/cache/tbi_sensitive"
OUTDIR_CLUMPED = BASE + "/20B3/work/data/outcome_extract_clumped"
CACHE_CLUMPED = BASE + "/20B3/work/cache/tbi_clumped"
MAP_CLUMPED = BASE + "/20B3/work/data/lead_snp_hg38_map_clumped.csv"
MAP_PRIMARY = BASE + "/20B3/work/data/lead_snp_hg38_map.csv"
MAP_SENSITIVE = BASE + "/20B3/work/data/lead_snp_hg38_map_sensitive.csv"
for d in (OUTDIR_PRIMARY, OUTDIR_SENSITIVE, OUTDIR_CLUMPED, CACHE_PRIMARY,
          CACHE_SENSITIVE, CACHE_CLUMPED):
    os.makedirs(d, exist_ok=True)

URL_BASE = "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release"

ENDPOINTS = {
    "SEPSIS": "AB1_OTHER_SEPSIS", "SEPTIC_STR": "AB1_STREPTO_SEPSIS",
    "SEPTIC_PNEU": "PNEUMO_SEPTHICHE", "PUERP_SEPS": "O15_PUERP_SEPSIS",
    "NEONAT_SEPS": "P16_BACTERIAL_SEPSIS_NEWBO", "PNEUMONIA": "J10_PNEUMONIA",
    "PNEUMOBACT": "J10_PNEUMOBACT", "VIRALPNEUMO": "J10_VIRALPNEUMO",
    "ARDS": "J10_ARDS", "PYELONEPHR": "N14_PYELONEPHR",
    "MENINGBACT": "G6_MENINGBACT", "ENDOCARD": "I9_ENDOCARD",
    "OSTEOMYEL": "M13_OSTEOMYELITIS", "CELLULITIS": "L12_CELLULITIS",
    "ABSCESS_CUT": "L12_ABSCESS_CUT", "ERYSIPELAS": "AB1_ERYSIPELAS",
    "INTEST_INF": "AB1_INTESTINAL_INFECTIONS",
    "GASTROENT": "AB1_GASTROENTERITIS_NOS",
    "BACT_OTHER": "AB1_OTHER_BACTERIAL", "BACT_NOS": "AB1_BACTINF_NOS",
    "ACTINOMYC": "AB1_ACTINOMYCOSIS", "COVID19": "U22_COVID19_CONFIRMED",
    "ZOSTER": "AB1_ZOSTER",
}

_OUT_FILES = {
    "SEPSIS": OUTDIR_PRIMARY + "/outcome_extract_SEPSIS.tsv",
    "SEPTIC_STR": OUTDIR_PRIMARY + "/outcome_extract_SEPTIC_STR.tsv",
    "SEPTIC_PNEU": OUTDIR_PRIMARY + "/outcome_extract_SEPTIC_PNEU.tsv",
    "PUERP_SEPS": OUTDIR_PRIMARY + "/outcome_extract_PUERP_SEPS.tsv",
    "NEONAT_SEPS": OUTDIR_PRIMARY + "/outcome_extract_NEONAT_SEPS.tsv",
    "PNEUMONIA": OUTDIR_PRIMARY + "/outcome_extract_PNEUMONIA.tsv",
    "PNEUMOBACT": OUTDIR_PRIMARY + "/outcome_extract_PNEUMOBACT.tsv",
    "VIRALPNEUMO": OUTDIR_PRIMARY + "/outcome_extract_VIRALPNEUMO.tsv",
    "ARDS": OUTDIR_PRIMARY + "/outcome_extract_ARDS.tsv",
    "PYELONEPHR": OUTDIR_PRIMARY + "/outcome_extract_PYELONEPHR.tsv",
    "MENINGBACT": OUTDIR_PRIMARY + "/outcome_extract_MENINGBACT.tsv",
    "ENDOCARD": OUTDIR_PRIMARY + "/outcome_extract_ENDOCARD.tsv",
    "OSTEOMYEL": OUTDIR_PRIMARY + "/outcome_extract_OSTEOMYEL.tsv",
    "CELLULITIS": OUTDIR_PRIMARY + "/outcome_extract_CELLULITIS.tsv",
    "ABSCESS_CUT": OUTDIR_PRIMARY + "/outcome_extract_ABSCESS_CUT.tsv",
    "ERYSIPELAS": OUTDIR_PRIMARY + "/outcome_extract_ERYSIPELAS.tsv",
    "INTEST_INF": OUTDIR_PRIMARY + "/outcome_extract_INTEST_INF.tsv",
    "GASTROENT": OUTDIR_PRIMARY + "/outcome_extract_GASTROENT.tsv",
    "BACT_OTHER": OUTDIR_PRIMARY + "/outcome_extract_BACT_OTHER.tsv",
    "BACT_NOS": OUTDIR_PRIMARY + "/outcome_extract_BACT_NOS.tsv",
    "ACTINOMYC": OUTDIR_PRIMARY + "/outcome_extract_ACTINOMYC.tsv",
    "COVID19": OUTDIR_PRIMARY + "/outcome_extract_COVID19.tsv",
    "ZOSTER": OUTDIR_PRIMARY + "/outcome_extract_ZOSTER.tsv",
}
_OUT_FILES_SENS = {
    "SEPSIS": OUTDIR_SENSITIVE + "/outcome_extract_SEPSIS.tsv",
    "SEPTIC_STR": OUTDIR_SENSITIVE + "/outcome_extract_SEPTIC_STR.tsv",
    "SEPTIC_PNEU": OUTDIR_SENSITIVE + "/outcome_extract_SEPTIC_PNEU.tsv",
    "PUERP_SEPS": OUTDIR_SENSITIVE + "/outcome_extract_PUERP_SEPS.tsv",
    "NEONAT_SEPS": OUTDIR_SENSITIVE + "/outcome_extract_NEONAT_SEPS.tsv",
    "PNEUMONIA": OUTDIR_SENSITIVE + "/outcome_extract_PNEUMONIA.tsv",
    "PNEUMOBACT": OUTDIR_SENSITIVE + "/outcome_extract_PNEUMOBACT.tsv",
    "VIRALPNEUMO": OUTDIR_SENSITIVE + "/outcome_extract_VIRALPNEUMO.tsv",
    "ARDS": OUTDIR_SENSITIVE + "/outcome_extract_ARDS.tsv",
    "PYELONEPHR": OUTDIR_SENSITIVE + "/outcome_extract_PYELONEPHR.tsv",
    "MENINGBACT": OUTDIR_SENSITIVE + "/outcome_extract_MENINGBACT.tsv",
    "ENDOCARD": OUTDIR_SENSITIVE + "/outcome_extract_ENDOCARD.tsv",
    "OSTEOMYEL": OUTDIR_SENSITIVE + "/outcome_extract_OSTEOMYEL.tsv",
    "CELLULITIS": OUTDIR_SENSITIVE + "/outcome_extract_CELLULITIS.tsv",
    "ABSCESS_CUT": OUTDIR_SENSITIVE + "/outcome_extract_ABSCESS_CUT.tsv",
    "ERYSIPELAS": OUTDIR_SENSITIVE + "/outcome_extract_ERYSIPELAS.tsv",
    "INTEST_INF": OUTDIR_SENSITIVE + "/outcome_extract_INTEST_INF.tsv",
    "GASTROENT": OUTDIR_SENSITIVE + "/outcome_extract_GASTROENT.tsv",
    "BACT_OTHER": OUTDIR_SENSITIVE + "/outcome_extract_BACT_OTHER.tsv",
    "BACT_NOS": OUTDIR_SENSITIVE + "/outcome_extract_BACT_NOS.tsv",
    "ACTINOMYC": OUTDIR_SENSITIVE + "/outcome_extract_ACTINOMYC.tsv",
    "COVID19": OUTDIR_SENSITIVE + "/outcome_extract_COVID19.tsv",
    "ZOSTER": OUTDIR_SENSITIVE + "/outcome_extract_ZOSTER.tsv",
}
_OUT_FILES_CLUMPED = {
    "SEPSIS": OUTDIR_CLUMPED + "/outcome_extract_SEPSIS.tsv",
    "SEPTIC_STR": OUTDIR_CLUMPED + "/outcome_extract_SEPTIC_STR.tsv",
    "SEPTIC_PNEU": OUTDIR_CLUMPED + "/outcome_extract_SEPTIC_PNEU.tsv",
    "PUERP_SEPS": OUTDIR_CLUMPED + "/outcome_extract_PUERP_SEPS.tsv",
    "NEONAT_SEPS": OUTDIR_CLUMPED + "/outcome_extract_NEONAT_SEPS.tsv",
    "PNEUMONIA": OUTDIR_CLUMPED + "/outcome_extract_PNEUMONIA.tsv",
    "PNEUMOBACT": OUTDIR_CLUMPED + "/outcome_extract_PNEUMOBACT.tsv",
    "VIRALPNEUMO": OUTDIR_CLUMPED + "/outcome_extract_VIRALPNEUMO.tsv",
    "ARDS": OUTDIR_CLUMPED + "/outcome_extract_ARDS.tsv",
    "PYELONEPHR": OUTDIR_CLUMPED + "/outcome_extract_PYELONEPHR.tsv",
    "MENINGBACT": OUTDIR_CLUMPED + "/outcome_extract_MENINGBACT.tsv",
    "ENDOCARD": OUTDIR_CLUMPED + "/outcome_extract_ENDOCARD.tsv",
    "OSTEOMYEL": OUTDIR_CLUMPED + "/outcome_extract_OSTEOMYEL.tsv",
    "CELLULITIS": OUTDIR_CLUMPED + "/outcome_extract_CELLULITIS.tsv",
    "ABSCESS_CUT": OUTDIR_CLUMPED + "/outcome_extract_ABSCESS_CUT.tsv",
    "ERYSIPELAS": OUTDIR_CLUMPED + "/outcome_extract_ERYSIPELAS.tsv",
    "INTEST_INF": OUTDIR_CLUMPED + "/outcome_extract_INTEST_INF.tsv",
    "GASTROENT": OUTDIR_CLUMPED + "/outcome_extract_GASTROENT.tsv",
    "BACT_OTHER": OUTDIR_CLUMPED + "/outcome_extract_BACT_OTHER.tsv",
    "BACT_NOS": OUTDIR_CLUMPED + "/outcome_extract_BACT_NOS.tsv",
    "ACTINOMYC": OUTDIR_CLUMPED + "/outcome_extract_ACTINOMYC.tsv",
    "COVID19": OUTDIR_CLUMPED + "/outcome_extract_COVID19.tsv",
    "ZOSTER": OUTDIR_CLUMPED + "/outcome_extract_ZOSTER.tsv",
}

FINNGEN_HEADER = ["#chrom", "pos", "ref", "alt", "rsids", "nearest_genes",
                  "pval", "mlogp", "beta", "sebeta", "af_alt",
                  "af_alt_cases", "af_alt_controls"]

_RS_PREFIX = "rs"

TOPHITS = BASE + "/20B3/work/data/mibiogen_tophits.tsv"


def load_snp_coords():
    """(bac, rsID) -> (chr, pos_grch37)   SNP      """
    m = {}
    with open(TOPHITS, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rs = (r.get("rsID") or "").strip()
            if not rs:
                continue
            m[(r["bac"], rs)] = (int(r["chr"]), int(r["pos"]))
    return m


def load_instruments(path, coords):
    """          snp_list    SNP        SNP      """
    recs = []
    missing = 0
    with open(path, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            for item in (r["snp_list"] or "").split(";"):
                p = item.split(":")
                if len(p) != 6:
                    continue
                c = coords.get((r["taxon"], p[0]))
                if c is None:
                    missing += 1
                    continue
                recs.append({
                    "taxon": r["taxon"], "level": r["level"], "name": r["name"],
                    "rsid": p[0], "beta": float(p[1]), "se": float(p[2]),
                    "pval": float(p[3]), "ea": p[4], "oa": p[5],
                    "chr37": c[0], "pos37": c[1],
                    "in_lct": r["in_lct"], "aligned": r["aligned"],
                })
    if missing:
        print("[03a] ⚠️ %d   SNP   tophits                 " % missing)
    return recs


def map_hg38(recs, cache_path, batch=180):
    """rsID → GRCh38    Ensembl REST `/variation/human` **POST   **  

         rsID       MiBioGen   GRCh37 FinnGen R12   GRCh38 
                145/146         rsID   Ensembl   GRCh38
                     `/map`               503     
         +      
    """
    cache = {}
    if os.path.isfile(cache_path):
        with open(cache_path, encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                cache[r["rsid"]] = (r["chrom"], r["pos38"])
    todo = sorted({r["rsid"] for r in recs
                   if r["rsid"].startswith(_RS_PREFIX) and r["rsid"] not in cache})
    print("[03a] rsID→GRCh38       %d       %d      %d "
          % (len(cache), len(todo), batch))
    if todo:
        hdr = {"Content-Type": "application/json", "Accept": "application/json",
               "User-Agent": "horzontal-pipeline 20B3 coord-map/1.0"}

        def one_batch(ids):
            body = json.dumps({"ids": ids}).encode()
            for attempt in range(4):
                try:
                    with open_url("https://rest.ensembl.org/variation/human",
                                  headers=hdr, timeout=120, data=body) as resp:
                        return json.loads(resp.read())
                except Exception as exc:  # noqa: BLE001
                    sys.stderr.write("[03a] Ensembl        %d    %s\n"
                                     % (attempt + 1, str(exc)[:120]))
                    time.sleep(4 * (attempt + 1))
            return {}

        for i in range(0, len(todo), batch):
            chunk = todo[i:i + batch]
            d = one_batch(chunk)
            for rs in chunk:
                v = d.get(rs)
                if not v:
                    cache[rs] = ("", "")
                    continue
                ms = [m for m in v.get("mappings", [])
                      if m.get("assembly_name") == "GRCh38"]
                cache[rs] = (str(ms[0]["seq_region_name"]), str(ms[0]["start"])) \
                    if ms else ("", "")
            print("        %d/%d         %d "
                  % (min(i + batch, len(todo)), len(todo),
                     sum(1 for rs in chunk if cache[rs][1])), flush=True)
        with open(cache_path, "w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["rsid", "chrom", "pos38"])
            for rs in sorted(cache):
                w.writerow([rs, cache[rs][0], cache[rs][1]])
        print("[03a]        %d   " % len(cache))

    miss = 0
    for r in recs:
        c, p = cache.get(r["rsid"], ("", ""))
        if not p:
            r["pos38"] = None
            r["chr38"] = None
            miss += 1
        else:
            r["pos38"] = int(p)
            r["chr38"] = int(c)
    if miss:
        print("[03a] ⚠️ %d   SNP   GRCh38               " % miss)
    return recs


def merged_windows(recs, pad):
    bychr = collections.defaultdict(list)
    for r in recs:
        if r.get("pos38"):
            bychr[r["chr38"]].append(r["pos38"])
    out = []
    for c, ps in bychr.items():
        ps = sorted(ps)
        s = e = ps[0]
        for p in ps[1:]:
            if p - e <= 2 * pad:
                e = p
            else:
                out.append((c, max(1, s - pad), e + pad))
                s = e = p
        out.append((c, max(1, s - pad), e + pad))
    return out


def extract(endpoint, recs, pad, cache_dir, out_path, tag, workers=8):
    pheno = ENDPOINTS[endpoint]
    url = URL_BASE + "/finngen_R12_" + pheno + ".gz"
    tr = TabixReader(url, cache_dir=cache_dir, name="finngen_R12_" + pheno + ".gz")
    wins = merged_windows(recs, pad)
    print("  [%s] %s %d     " % (endpoint, pheno, len(wins)))

    want_rs = collections.defaultdict(list)
    want_pos = collections.defaultdict(list)
    for r in recs:
        if r["rsid"].startswith(_RS_PREFIX):
            want_rs[r["rsid"]].append(r)
        if r.get("pos38"):
            want_pos[(r["chr38"], r["pos38"])].append(r)

    lock = threading.Lock()
    hits, pos_hits = {}, {}

    def do_win(w):
        c, s, e = w
        try:
            _hdr, rows = tr.query(str(c), s, e, header=FINNGEN_HEADER)
        except Exception as exc:  # noqa: BLE001
            return False, "%s: %s" % (type(exc).__name__, exc), []
        return True, "", rows

    n_ok = n_fail = 0
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(do_win, w) for w in wins]
        for fut in as_completed(futs):
            ok, err, rows = fut.result()
            if not ok:
                n_fail += 1
                sys.stderr.write("[03a]      %s\n" % err)
                continue
            n_ok += 1
            with lock:
                for row in rows:
                    if len(row) < len(FINNGEN_HEADER):
                        continue
                    d = dict(zip(FINNGEN_HEADER, row))
                    for rs in (d["rsids"] or "").split(","):
                        rs = rs.strip()
                        for r in want_rs.get(rs, ()):
                            hits[(r["taxon"], r["rsid"])] = d
                    try:
                        key = (int(d["#chrom"]), int(d["pos"]))
                    except Exception:  # noqa: BLE001
                        continue
                    if key in want_pos:
                        pos_hits[key] = d
    print("  [%s]       %d /    %d" % (endpoint, n_ok, n_fail))

    out = []
    for r in recs:
        d = hits.get((r["taxon"], r["rsid"]))
        match = "rsid"
        if d is None and r.get("pos38"):
            d = pos_hits.get((r["chr38"], r["pos38"]))
            match = "pos"
        if d is None:
            continue
        out.append({
            "endpoint": endpoint, "phenocode": pheno,
            "taxon": r["taxon"], "level": r["level"], "name": r["name"],
            "rsid": r["rsid"], "match_by": match,
            "chr": r["chr37"], "pos_grch37": r["pos37"],
            "pos_grch38": r.get("pos38") or "",
            "ea_instr": r["ea"], "oa_instr": r["oa"],
            "beta_instr": r["beta"], "se_instr": r["se"], "p_instr": r["pval"],
            "ref_out": d["ref"], "alt_out": d["alt"],
            "beta_out": d["beta"], "se_out": d["sebeta"], "pval_out": d["pval"],
            "eaf_out": d["af_alt"], "in_lct": r["in_lct"], "aligned": r["aligned"],
        })
    with open(out_path, "w", encoding="utf-8", newline="") as fh:
        if out:
            w = csv.DictWriter(fh, fieldnames=list(out[0].keys()), delimiter="\t")
            w.writeheader(); w.writerows(out)
    print("  [%s][%s]    %d  " % (tag, endpoint, len(out)))
    return len(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--endpoints", default="")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--pad", type=int, default=100000)
    ap.add_argument("--tag", default="primary",
                    choices=("primary", "sensitive", "clumped"))
    a = ap.parse_args()

    if a.tag == "primary":
        instr, cache_dir, map_cache, table = (
            INSTR, CACHE_PRIMARY, MAP_PRIMARY, _OUT_FILES)
    elif a.tag == "clumped":
        instr, cache_dir, map_cache, table = (
            INSTR_CLUMPED, CACHE_CLUMPED, MAP_CLUMPED, _OUT_FILES_CLUMPED)
    else:
        instr, cache_dir, map_cache, table = (
            INSTR_SENS, CACHE_SENSITIVE, MAP_SENSITIVE, _OUT_FILES_SENS)

    if a.endpoints:
        keys = [k for k in a.endpoints.split(",") if k]
    else:
        keys = list(ENDPOINTS)
    bad = [k for k in keys if k not in ENDPOINTS or k not in table]
    if bad:
        raise SystemExit("             %s" % bad)

    recs = load_instruments(instr, load_snp_coords())
    print("[03a] [%s]    SNP %d   / %d taxa / %d    rsID"
          % (a.tag, len(recs), len({r['taxon'] for r in recs}),
             len({r['rsid'] for r in recs})))
    recs = map_hg38(recs, map_cache)

    t0 = time.time()
    for i, k in enumerate(keys, 1):
        extract(k, recs, a.pad, cache_dir, table[k], a.tag, a.workers)
        print("[03a] %d/%d %s       %.1f    "
              % (i, len(keys), k, (time.time() - t0) / 60))


if __name__ == "__main__":
    main()
