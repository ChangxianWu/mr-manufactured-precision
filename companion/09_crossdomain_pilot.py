# -*- coding: utf-8 -*-
"""Metabolite cross-domain pilot: fetch, map, remote outcome extraction, clumping, IVW, k_eff/f and pre-declared criteria."""
import argparse
import collections
import csv
import gzip
import importlib.util
import io
import json
import math
import os
import sys
import time
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

BASE = "D:/Horzontal_project"
sys.path.insert(0, BASE)
from cv_shared.net_guard import fetch, open_url  # noqa: E402
from cv_shared.bgzf_tabix import TabixReader  # noqa: E402

ROOT = Path(BASE) / "20B3" / "work" / "crossdomain"
DATA = ROOT / "data"
OUT = ROOT / "output"
DATA.mkdir(parents=True, exist_ok=True)
OUT.mkdir(parents=True, exist_ok=True)

CLUMP_MODULE = BASE + "/20B3/work/code/02b_clump_python.py"
RS_CACHE = DATA / "rsid_pos38.tsv"

FTP = ("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
       "GCST90132001-GCST90133000/%s/%s_buildGRCh37.tsv.gz")

TRAITS = [
    ("Glycoprotein acetyls", "GCST90132728"),
    ("Serum total cholesterol", "GCST90132746"),
    ("Apolipoprotein B", "GCST90132701"),
    ("Isoleucine", "GCST90132730"),
    ("Lactate", "GCST90132731"),
    ("Phenylalanine", "GCST90132740"),
    ("Pyruvate", "GCST90132743"),
    ("3-Hydroxybutyrate", "GCST90132695"),
    ("Apolipoprotein A1", "GCST90132700"),
    ("Glucose", "GCST90132724"),
]

ENDPOINTS = {
    "SEPSIS": "AB1_OTHER_SEPSIS",
    "PNEUMONIA": "J10_PNEUMONIA",
    "GASTROENT": "AB1_GASTROENTERITIS_NOS",
    "PNEUMOBACT": "J10_PNEUMOBACT",
    "VIRALPNEUMO": "J10_VIRALPNEUMO",
    "ARDS": "J10_ARDS",
    "CELLULITIS": "L12_CELLULITIS",
    "ABSCESS_CUT": "L12_ABSCESS_CUT",
    "ERYSIPELAS": "AB1_ERYSIPELAS",
    "INTEST_INF": "AB1_INTESTINAL_INFECTIONS",
    "COVID19": "U22_COVID19_CONFIRMED",
}
FINNGEN = "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_%s.gz"

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
PAD = 10000
KB_CLUMP = 10000
R2_MAX = 0.001
CLUMP_R2 = R2_MAX

csv.field_size_limit(10 ** 8)


# ---------------------------------------------------------------- helpers
def load_clump_module():
    spec = importlib.util.spec_from_file_location("clump02b", CLUMP_MODULE)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def write_rows(path, header, rows):
    buf = io.StringIO()
    w = csv.writer(buf, delimiter="\t", lineterminator="\n")
    w.writerow(header)
    w.writerows(rows)
    Path(path).write_text(buf.getvalue(), encoding="utf-8", newline="")


def read_rows(path):
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ---------------------------------------------------------------- fetch
def remote_size(url):
    """  Range: bytes=0-0   Content-Range     EBI   HEAD     """
    for _ in range(3):
        try:
            with open_url(url, headers={"Range": "bytes=0-0"}, timeout=60) as r:
                cr = r.headers.get("Content-Range") or ""
                if "/" in cr:
                    return int(cr.rsplit("/", 1)[1])
        except Exception:  # noqa: BLE001
            time.sleep(2)
    return None


def download_ranged(url, dst, chunk=4 << 20, workers=4):
    """   Range          100 MB    EBI   TLS         SSLEOFError  
                     """
    from concurrent.futures import ThreadPoolExecutor
    total = remote_size(url)
    if not total:
        body = fetch(url, timeout=300, tries=3)
        if not body:
            return False
        dst.write_bytes(body)
        return True
    part_dir = dst.with_suffix(dst.suffix + ".parts")
    part_dir.mkdir(exist_ok=True)
    starts = list(range(0, total, chunk))
    have = [s for s in starts
            if (part_dir / ("%d.bin" % s)).exists()
            and (part_dir / ("%d.bin" % s)).stat().st_size
            == min(s + chunk, total) - s]
    print("[fetch]   %s: %.0f MB    %d      %d    %d"
          % (dst.name, total / 1e6, len(starts), len(have), len(starts) - len(have)),
          flush=True)

    def one(start):
        fp = part_dir / ("%d.bin" % start)
        end = min(start + chunk, total) - 1
        if fp.exists() and fp.stat().st_size == end - start + 1:
            return True
        for attempt in range(5):
            try:
                with open_url(url, headers={"Range": "bytes=%d-%d" % (start, end)},
                              timeout=120) as r:
                    b = r.read()
                if len(b) == end - start + 1:
                    fp.write_bytes(b)
                    return True
            except Exception:  # noqa: BLE001
                time.sleep(2 * (attempt + 1))
        return False

    ok = True
    with ThreadPoolExecutor(max_workers=workers) as ex:
        for good in ex.map(one, starts):
            ok = ok and good
    if not ok:
        return False
    blob = bytearray()
    for start in starts:
        blob.extend((part_dir / ("%d.bin" % start)).read_bytes())
    dst.write_bytes(bytes(blob))
    print("[fetch]   %s      %.1f MB" % (dst.name, len(blob) / 1e6))
    for p in part_dir.iterdir():
        p.unlink()
    part_dir.rmdir()
    return True


def fetch_traits(only=None):
    traits = [t for t in TRAITS if not only or t[1] in only]
    for name, gcst in traits:
        url = FTP % (gcst, gcst)
        dst = DATA / ("%s.tsv.gz" % gcst)
        if dst.exists() and dst.stat().st_size > 200_000_000:
            print("[fetch]        %s %s" % (gcst, name))
            continue
        t0 = time.time()
        ok = download_ranged(url, dst)
        if not ok or not dst.exists():
            print("[fetch] ✗ %s     " % gcst)
            continue
        print("[fetch] %s %s  %.1f MB  %.1fs"
              % (gcst, name, dst.stat().st_size / 1e6, time.time() - t0))

    for name, gcst in traits:
        src = DATA / ("%s.tsv.gz" % gcst)
        dst = DATA / ("%s.sig.tsv" % gcst)
        if not src.exists() or dst.exists():
            continue
        rows = []
        with gzip.open(src, "rt", errors="replace") as fh:
            rd = csv.reader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)
            hdr = next(rd)
            idx = {k: i for i, k in enumerate(hdr)}
            for f in rd:
                if len(f) != len(hdr):
                    continue
                try:
                    p = float(f[idx["p_value"]])
                except ValueError:
                    continue
                if p >= 5e-8:
                    continue
                rows.append([f[idx["chromosome"]], f[idx["base_pair_location"]],
                             f[idx["rs_id"]], f[idx["effect_allele"]],
                             f[idx["other_allele"]],
                             f[idx["effect_allele_frequency"]], f[idx["beta"]],
                             f[idx["standard_error"]], f[idx["p_value"]]])
        best = {}
        for row in rows:
            rs = row[2]
            if not rs or not rs.startswith("rs"):
                continue
            if rs not in best or float(row[8]) < float(best[8]):
                best[rs] = row
        out = sorted(best.values(), key=lambda x: (int(x[0]), int(x[1])))
        write_rows(dst, ["chr", "pos37", "rsid", "ea", "oa", "eaf",
                         "beta", "se", "p"], out)
        print("[fetch] %s      %d    %d " % (gcst, len(out), len(rows)))


# ---------------------------------------------------------------- map
def map_rsids():
    cache = {}
    if RS_CACHE.exists():
        for r in read_rows(RS_CACHE):
            cache[r["rsid"]] = (r["chrom"], r["pos38"])
    todo = set()
    for _n, gcst in TRAITS:
        p = DATA / ("%s.sig.tsv" % gcst)
        if not p.exists():
            continue
        for r in read_rows(p):
            if r["rsid"] not in cache:
                todo.add(r["rsid"])
    todo = sorted(todo)
    print("[map]    %d     %d" % (len(cache), len(todo)))
    if todo:
        hdr = {"Content-Type": "application/json", "Accept": "application/json",
               "User-Agent": "horzontal-pipeline 20B3 crossdomain/1.0"}

        def one_batch(ids):
            body = json.dumps({"ids": ids}).encode()
            for attempt in range(4):
                try:
                    with open_url("https://rest.ensembl.org/variation/human",
                                  headers=hdr, timeout=180, data=body) as resp:
                        return json.loads(resp.read())
                except Exception as exc:  # noqa: BLE001
                    sys.stderr.write("[map]      %d: %s\n" % (attempt + 1, str(exc)[:100]))
                    time.sleep(4 * (attempt + 1))
            return {}

        batch = 180
        for i in range(0, len(todo), batch):
            chunk = todo[i:i + batch]
            d = one_batch(chunk)
            for rs in chunk:
                v = d.get(rs) or {}
                ms = [m for m in v.get("mappings", [])
                      if m.get("assembly_name") == "GRCh38"]
                cache[rs] = (str(ms[0]["seq_region_name"]), str(ms[0]["start"])) if ms else ("", "")
            print("      %d/%d" % (min(i + batch, len(todo)), len(todo)), flush=True)
        write_rows(RS_CACHE, ["rsid", "chrom", "pos38"],
                   [[k, v[0], v[1]] for k, v in sorted(cache.items())])
    else:
        write_rows(RS_CACHE, ["rsid", "chrom", "pos38"],
                   [[k, v[0], v[1]] for k, v in sorted(cache.items())])
    hit = sum(1 for v in cache.values() if v[1])
    print("[map]      %d/%d   GRCh38   " % (hit, len(cache)))


# ---------------------------------------------------------------- outcome
def merged_windows(positions):
    pos = sorted(set(positions))
    out = []
    for p in pos:
        if out and p - out[-1][1] <= PAD:
            out[-1][1] = p + PAD
        else:
            out.append([p - PAD, p + PAD])
    return out


def outcome_step(only=None):
    """ **  **                TabixReader           trait 

           /ApoB      trait      ~15    11     × 8 trait        
    .tbi             reader             
    """
    cache = {r["rsid"]: (r["chrom"], r["pos38"]) for r in read_rows(RS_CACHE)}
    trait_spec = []
    for name, gcst in TRAITS:
        if only and gcst not in only:
            continue
        sig = DATA / ("%s.sig.tsv" % gcst)
        if not sig.exists():
            continue
        rows = read_rows(sig)
        want = {r["rsid"] for r in rows}
        by_chr = collections.defaultdict(list)
        for r in rows:
            c, p = cache.get(r["rsid"], ("", ""))
            if p:
                by_chr[int(c)].append(int(p))
        wins = {c: merged_windows(ps) for c, ps in by_chr.items()}
        trait_spec.append((name, gcst, want, wins))
    if not trait_spec:
        print("[outcome]      trait")
        return

    def run_endpoint(key, pheno):
        tr = TabixReader(FINNGEN % pheno, cache_dir=str(DATA / "tbi"), name=pheno)
        for name, gcst, want, wins in trait_spec:
            out_p = DATA / ("outcome_%s_%s.tsv" % (key, gcst))
            if out_p.exists():
                continue
            got = {}
            for chrom in sorted(wins):
                for lo, hi in wins[chrom]:
                    try:
                        _header, lines = tr.query(str(chrom), lo, hi)
                    except Exception as exc:  # noqa: BLE001
                        sys.stderr.write("[outcome] %s %s chr%d %d-%d   : %s\n"
                                         % (key, gcst, chrom, lo, hi, str(exc)[:70]))
                        continue
                    for line in lines:
                        f = line if isinstance(line, list) else line.rstrip("\n").split("\t")
                        if len(f) < 13 or f[0].startswith("#"):
                            continue
                        for rs in (f[4] or "").split(","):
                            if rs in want and rs not in got:
                                got[rs] = [f[0], f[1], f[2], f[3], f[4], f[6],
                                           f[8], f[9], f[10]]
            write_rows(out_p, ["chrom38", "pos38", "ref", "alt", "rsids", "pval",
                               "beta", "sebeta", "af_alt"], sorted(got.values()))
            print("[outcome] %-12s %-12s    %d/%d" % (key, gcst, len(got), len(want)),
                  flush=True)

    from concurrent.futures import ThreadPoolExecutor
    workers = min(5, len(ENDPOINTS))
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futures = [ex.submit(run_endpoint, key, pheno) for key, pheno in ENDPOINTS.items()]
        for fut in futures:
            fut.result()


# ---------------------------------------------------------------- MR
def align(ea, oa, ref, alt, beta_out, eaf_out):
    ea, oa, ref, alt = (ea or "").upper(), (oa or "").upper(), (ref or "").upper(), (alt or "").upper()
    if not (ea and oa and ref and alt):
        return None, "missing"
    if COMP.get(ea) == oa and eaf_out is not None and abs(eaf_out - 0.5) < 0.05:
        return None, "palin_ambig"
    if ea == alt and oa == ref:
        return beta_out, "ok"
    if ea == ref and oa == alt:
        return -beta_out, "flip"
    if COMP.get(ea) == oa:
        return beta_out, "palin"
    return None, "mismatch"


def ivw(pairs):
    sw = sw2 = 0.0
    for b, se in pairs:
        w = 1.0 / (se * se)
        sw += w
        sw2 += w * b
    if sw <= 0:
        return None, None
    return sw2 / sw, math.sqrt(1.0 / sw)


def _scan_chrom(args):
    """              ProcessPoolExecutor       

          chr1   2–3    22       1            
    """
    chrom, positions = args
    clump = load_clump_module()
    keep = clump.load_keep()
    names = clump.sample_columns(chrom)
    n2i = {n: i for i, n in enumerate(names)}
    keep_idx = [n2i[s] for s in keep if s in n2i]
    got, n_seen = clump.read_genotypes(chrom, positions, keep_idx)
    return chrom, got, n_seen, len(positions)


def mr_step(only=None, suffix=""):
    clump = load_clump_module()
    keep = clump.load_keep()
    cache = {r["rsid"]: (r["chrom"], r["pos38"]) for r in read_rows(RS_CACHE)}

    exp_by_trait = {}
    union = collections.defaultdict(set)
    for name, gcst in TRAITS:
        if only and gcst not in only:
            continue
        sig = DATA / ("%s.sig.tsv" % gcst)
        if not sig.exists():
            continue
        exp = {r["rsid"]: r for r in read_rows(sig)}
        exp_by_trait[(name, gcst)] = exp
        for rs, r in exp.items():
            union[int(r["chr"])].add(int(r["pos37"]))
    geno = {}
    items = [(chrom, union[chrom]) for chrom in sorted(union)]
    n_workers = min(5, len(items))
    if n_workers > 1:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=n_workers) as ex:
            for chrom, got, n_seen, n_want in ex.map(_scan_chrom, items):
                geno.update(got)
                print("      chr%-2d    %d      %d / %d   "
                      % (chrom, n_seen, len(got), n_want), flush=True)
    else:
        for chrom, positions in items:
            c, got, n_seen, n_want = _scan_chrom((chrom, positions))
            geno.update(got)
            print("      chr%-2d    %d      %d / %d   "
                  % (c, n_seen, len(got), n_want), flush=True)

    pair_rows, var_rows = [], []
    for (name, gcst), exp in exp_by_trait.items():
        cands_all = [(rs, int(r["chr"]), int(r["pos37"]), float(r["p"]))
                     for rs, r in exp.items()]
        by_chr_pos = collections.defaultdict(list)
        for _rs, c, p, _pv in cands_all:
            by_chr_pos[c].append(p)
        best_share, k_all = 0.0, len(cands_all)
        for c, ps in by_chr_pos.items():
            ps = sorted(ps)
            j = 0
            for i in range(len(ps)):
                while ps[i] - ps[j] > KB_CLUMP * 1000:
                    j += 1
                best_share = max(best_share, (i - j + 1) / k_all)
        if k_all == 0:
            continue

        for key in ENDPOINTS:
            outp = DATA / ("outcome_%s_%s.tsv" % (key, gcst))
            if not outp.exists():
                continue
            obs = {r["rsids"]: r for r in read_rows(outp)}
            usable = []
            for rs, r in exp.items():
                o = obs.get(rs)
                if not o:
                    continue
                try:
                    b_out = float(o["beta"]); se_out = float(o["sebeta"])
                    eaf = float(o["af_alt"]) if o["af_alt"] not in ("", "NA") else None
                    bx = float(r["beta"])
                except (TypeError, ValueError):
                    continue
                ba, st = align(r["ea"], r["oa"], o["ref"], o["alt"], b_out, eaf)
                if ba is None or st not in ("ok", "flip") or se_out <= 0 or bx == 0:
                    continue
                b_wald = ba / bx
                se_wald = se_out / abs(bx)
                usable.append({"rsid": rs, "chr": int(r["chr"]),
                               "pos": int(r["pos37"]), "p": float(r["p"]),
                               "b": b_wald, "se": se_wald,
                               "F": (bx / float(r["se"])) ** 2})
                var_rows.append([gcst, key, rs, r["p"], b_wald, se_wald, st])
            if not usable:
                continue

            b_u, se_u = ivw([(u["b"], u["se"]) for u in usable])
            if b_u is None:
                continue
            cands = [(u["rsid"], u["chr"], u["pos"], u["p"]) for u in usable]
            kept, _dropped = clump.clump_taxon(cands, geno, CLUMP_R2, KB_CLUMP)
            kept_rs = {k[0] for k in kept}
            kept_u = [u for u in usable if u["rsid"] in kept_rs]
            if not kept_u:
                continue
            b_c, se_c = ivw([(u["b"], u["se"]) for u in kept_u])
            k = len(usable)
            kk = len(kept_u)
            keff = (se_c / se_u) ** 2 if se_u > 0 else None
            f = (se_c / se_u) / math.sqrt(k) if se_u > 0 else None
            medF = median([u["F"] for u in usable])
            nF_lt10 = sum(1 for u in usable if u["F"] < 10)

            lead = min(usable, key=lambda u: u["p"])
            lg = geno.get(str(lead["pos"]))
            block = []
            for u in usable:
                if u["chr"] != lead["chr"]:
                    continue
                if abs(u["pos"] - lead["pos"]) > KB_CLUMP * 1000:
                    continue
                if u["rsid"] == lead["rsid"]:
                    block.append(u)
                    continue
                r2 = clump.r2_of(lg or [], geno.get(str(u["pos"])) or [])
                if r2 is not None and r2 > CLUMP_R2:
                    block.append(u)
            f_block = ""
            k_block = 0
            if len(block) >= 2:
                _bb, se_b = ivw([(u["b"], u["se"]) for u in block])
                if se_b and se_b > 0:
                    k_block = len(block)
                    f_block = round((lead["se"] / se_b) / math.sqrt(k_block), 4)
            pair_rows.append([name, gcst, key, k, kk, round(se_u, 6), round(se_c, 6),
                              round(keff, 4), round(f, 4), round(best_share, 3),
                              round(b_u, 4), round(b_c, 4), round(medF, 2),
                              nF_lt10, k_block, f_block])

    write_rows(OUT / ("variants%s.csv" % suffix),
               ["trait_gcst", "endpoint", "rsid", "p_exp", "b_wald", "se_wald", "align"],
               var_rows)
    write_rows(OUT / ("pairs%s.csv" % suffix),
               ["trait", "gcst", "endpoint", "k", "k_kept", "se_unclumped",
                "se_clumped", "k_eff", "f", "max_window_share", "b_unclumped",
                "b_clumped", "median_F", "n_F_lt10", "k_block", "f_block"], pair_rows)
    print("[mr]     %d    -> %s" % (len(pair_rows), OUT / ("pairs%s.csv" % suffix)))
    summarize(suffix)


def summarize(suffix=""):
    rows = read_rows(OUT / ("pairs%s.csv" % suffix))
    if not rows:
        print("[mr]    ")
        return
    by_trait = collections.defaultdict(list)
    for r in rows:
        by_trait[(r["trait"], r["gcst"])].append(r)
    out = []
    for (t, g), rs in sorted(by_trait.items()):
        fs = sorted(float(x["f"]) for x in rs)
        ks = sorted(float(x["k_eff"]) / float(x["k"]) for x in rs)
        out.append([t, g, len(rs), sum(int(x["k"]) for x in rs) / len(rs),
                    round(median(fs), 4), round(median(ks), 4),
                    round(float(rs[0]["max_window_share"]), 3)])
    write_rows(OUT / ("summary%s.csv" % suffix),
               ["trait", "gcst", "pairs", "k_mean", "f_median", "keff_over_k_median",
                "max_window_share"], out)
    for r in out:
        print("[sum] %-24s pairs=%d k=%.1f f_med=%.3f keff/k=%.3f single_locus=%.2f"
              % (r[0], r[2], r[3], r[4], r[5], r[6]))


def spearman(x, y):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = rank(x), rank(y)
    n = len(x)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    return num / (dx * dy) if dx and dy else float("nan")


def stats_step(suffix=""):
    """       P1–P5             """
    rows = [r for r in read_rows(OUT / ("pairs%s.csv" % suffix)) if int(r["k"]) >= 2]
    if not rows:
        print("[stats]   k>=2           ")
        return
    f = [float(r["f"]) for r in rows]
    share = [float(r["max_window_share"]) for r in rows]
    k = [float(r["k"]) for r in rows]
    kk = [float(r["k_kept"]) for r in rows]
    ratio = [float(r["k_eff"]) / float(r["k"]) for r in rows]

    rho1 = spearman(share, f)
    hi = [v for v, s in zip(f, share) if s >= 0.8]
    lo = [v for v, s in zip(f, share) if s < 0.8]
    rho3 = spearman([math.log10(x) for x in k], [math.log10(max(x, 1e-9)) for x in ratio])
    rho4 = spearman(kk, f)

    pairs5 = [(float(r["f"]), float(r["f_block"])) for r in rows
              if r.get("f_block") not in ("", None)]
    wins = sum(1 for a, b in pairs5 if b > a)
    diffs = [b - a for a, b in pairs5]
    n5 = len(pairs5)
    if n5:
        from math import comb
        kk_ = min(wins, n5 - wins)
        p_sign = min(1.0, 2 * sum(comb(n5, i) for i in range(0, kk_ + 1)) / (2 ** n5))
    else:
        p_sign = float("nan")

    lines = [
        "#           Kettunen NMR     × FinnGen R12      ",
        "",
        "            P1    ↔f     P2    trait   f    "
        "P3 k↔k_eff/k     P4    block  ↔f     P5   block     f     trait ",
        "",
        "|    |     |    |    |",
        "|---|---|---|---|",
        "| P1 | Spearman(     , f) | %.3f n=%d  | %s |"
        % (rho1, len(rows), "  " if rho1 > 0 else "   "),
        "| P2 | f       %.3f n=%d  vs    %.3f n=%d  |    %.3f | %s |"
        % (median(hi) if hi else float("nan"), len(hi),
           median(lo) if lo else float("nan"), len(lo),
           (median(hi) - median(lo)) if hi and lo else float("nan"),
           ("  " if median(hi) > median(lo) else "   ") if hi and lo
           else "         trait "),
        "| P3 | Spearman(log k, log k_eff/k) | %.3f | %s |"
        % (rho3, "  " if rho3 < 0 else "   "),
        "| P4 | Spearman(k_kept, f) | %.3f | %s |"
        % (rho4, "  " if rho4 < 0 else "   "),
        "| P5 | f_block > f     | %d/%d     %+.3f      p=%.3f | %s |"
        % (wins, n5, median(diffs) if diffs else float("nan"), p_sign,
           "  " if n5 and wins > n5 / 2 else ("   " if n5 else "    ")),
        "",
        "     k>=2  %d f    %.3f k_eff/k    %.3f "
        "   F     %.1f–%.1f F<10      %d   "
        % (len(rows), median(f), median(ratio),
           min(float(r["median_F"]) for r in rows),
           max(float(r["median_F"]) for r in rows),
           sum(int(r["n_F_lt10"]) for r in rows)),
        "",
        "| trait |    | k | k_kept | SE  | k_eff | f | f_block | k_block |   F |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        lines.append("| %s | %s | %s | %s | %.2f | %.1f | %.3f | %s | %s | %.1f |"
                     % (r["trait"], r["endpoint"], r["k"], r["k_kept"],
                        float(r["se_clumped"]) / float(r["se_unclumped"]),
                        float(r["k_eff"]), float(r["f"]),
                        r["f_block"] if r["f_block"] not in ("", None) else "—",
                        r["k_block"] if r["k_block"] not in ("", None) else "—",
                        float(r["median_F"])))
    lines.append("")
    txt = "\n".join(lines)
    (OUT / ("stats%s.md" % suffix)).write_text(txt, encoding="utf-8", newline="")
    print(txt)


def median(v):
    v = sorted(v)
    n = len(v)
    if n == 0:
        return float("nan")
    return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])


def main():
    global CLUMP_R2
    ap = argparse.ArgumentParser()
    ap.add_argument("--step", required=True,
                    choices=["fetch", "map", "outcome", "mr", "stats", "all"])
    ap.add_argument("--only", default="",
                    help="      GCST         trait      ")
    ap.add_argument("--r2", type=float, default=R2_MAX,
                    help="clumping   r²              0.001 ")
    ap.add_argument("--suffix", default="",
                    help="                      ")
    a = ap.parse_args()
    CLUMP_R2 = a.r2
    only = {x.strip() for x in a.only.split(",") if x.strip()} or None
    if a.step in ("fetch", "all"):
        fetch_traits(only)
    if a.step in ("map", "all"):
        map_rsids()
    if a.step in ("outcome", "all"):
        outcome_step(only)
    if a.step in ("mr", "all"):
        mr_step(only, a.suffix)
    if a.step == "stats":
        stats_step(a.suffix)


if __name__ == "__main__":
    main()
