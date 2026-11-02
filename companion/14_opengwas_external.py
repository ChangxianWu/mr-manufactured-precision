# -*- coding: utf-8 -*-
"""External outcome replication through OpenGWAS (UK Biobank sepsis and pneumonia)."""
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
from cv_shared.net_guard import open_url  # noqa: E402

ROOT = Path(BASE) / "20B3" / "work" / "crossdomain"
DATA = ROOT / "data"
OUT = ROOT / "output"
TOKEN_FILE = Path.home() / ".zcode" / "secrets" / "opengwas_token.txt"
API = "https://api.opengwas.io/api/"

ENDPOINTS = {
    "SEPSIS_UKB": "ieu-b-4980",          # UKB, N=486,484
    "PNEUMONIA_UKB": "ebi-a-GCST90018901",  # N=480,299
}
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
CHUNK = 50
PAD = 10000
R2_MAX = 0.001
KB_CLUMP = 10000
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def token():
    t = os.environ.get("OPENGAWS_TOKEN", "").strip()
    if t:
        return t
    if TOKEN_FILE.exists():
        return TOKEN_FILE.read_text(encoding="ascii").strip()
    raise SystemExit("   OpenGWAS token      OPENGAWS_TOKEN         ")


def api_post(path, body, tries=4):
    hdr = {"Authorization": "Bearer " + token(),
           "Content-Type": "application/json", "Accept": "application/json"}
    for attempt in range(tries):
        try:
            with open_url(API + path, headers=hdr,
                          data=json.dumps(body).encode(), timeout=120) as r:
                return json.loads(r.read())
        except Exception as exc:  # noqa: BLE001
            sys.stderr.write("[og] %s   %d    : %s\n"
                             % (path, attempt + 1, str(exc)[:90]))
            time.sleep(3 * (attempt + 1))
    return None


def read_rows(path):
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_rows(path, header, rows):
    buf = io.StringIO()
    w = csv.writer(buf, delimiter="\t", lineterminator="\n")
    w.writerow(header)
    w.writerows(rows)
    Path(path).write_text(buf.getvalue(), encoding="utf-8", newline="")


# ---------------------------------------------------------------- fetch
def fetch_step(only=None):
    for key, gid in ENDPOINTS.items():
        for name, gcst in TRAITS:
            out = DATA / ("og_%s_%s.tsv" % (key, gcst))
            if out.exists():
                continue
            sig = DATA / ("%s.sig.tsv" % gcst)
            if not sig.exists():
                continue
            exp = read_rows(sig)
            by_rs = {r["rsid"]: r for r in exp}
            by_pos = {}
            for r in exp:
                by_pos[(r["chr"], r["pos37"])] = r
            got = {}
            rsids = [r["rsid"] for r in exp]
            def fetch_chunk(chunk):
                """     400/            1     """
                res = api_post("associations", {"id": [gid], "variant": chunk}, tries=4)
                if res is not None:
                    return res
                if len(chunk) == 1:
                    return []
                mid = len(chunk) // 2
                return fetch_chunk(chunk[:mid]) + fetch_chunk(chunk[mid:])

            for i in range(0, len(rsids), CHUNK):
                chunk = rsids[i:i + CHUNK]
                res = fetch_chunk(chunk)
                if not res:
                    continue
                for a in res:
                    rs = a.get("rsid")
                    if rs in by_rs and rs not in got:
                        got[rs] = a
                    else:
                        key_pos = (str(a.get("chr")), str(a.get("position")))
                        if key_pos in by_pos:
                            orig = by_pos[key_pos]["rsid"]
                            if orig not in got:
                                got[orig] = a
                time.sleep(0.3)
            rows = []
            for rs, a in got.items():
                rows.append([rs, a.get("chr"), a.get("position"), a.get("ea"),
                             a.get("nea"), a.get("eaf"), a.get("beta"), a.get("se"),
                             a.get("p"), a.get("n")])
            write_rows(out, ["rsid", "chr", "pos", "ea", "nea", "eaf", "beta",
                             "se", "p", "n"], rows)
            print("[og] %-14s %-22s    %d/%d" % (key, gcst, len(rows), len(exp)))


# ---------------------------------------------------------------- MR
def _scan_chrom(args):
    chrom, positions = args
    spec = importlib.util.spec_from_file_location(
        "pilot09", str(ROOT / "09_crossdomain_pilot.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    clump = mod.load_clump_module()
    keep = clump.load_keep()
    names = clump.sample_columns(chrom)
    n2i = {n: i for i, n in enumerate(names)}
    keep_idx = [n2i[s] for s in keep if s in n2i]
    got, n_seen = clump.read_genotypes(chrom, positions, keep_idx)
    return chrom, got, n_seen, len(positions)


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
    return (sw2 / sw, math.sqrt(1.0 / sw)) if sw > 0 else (None, None)


def median(v):
    v = sorted(v)
    n = len(v)
    if n == 0:
        return float("nan")
    return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])


def mr_step():
    exp_by_trait = {}
    union = collections.defaultdict(set)
    for name, gcst in TRAITS:
        sig = DATA / ("%s.sig.tsv" % gcst)
        if not sig.exists():
            continue
        exp = {r["rsid"]: r for r in read_rows(sig)}
        exp_by_trait[(name, gcst)] = exp
        for r in exp.values():
            union[int(r["chr"])].add(int(r["pos37"]))
    from concurrent.futures import ProcessPoolExecutor
    geno = {}
    items = [(c, union[c]) for c in sorted(union)]
    with ProcessPoolExecutor(max_workers=min(5, len(items))) as ex:
        for chrom, got, n_seen, n_want in ex.map(_scan_chrom, items):
            geno.update(got)
            print("      chr%-2d    %d / %d" % (chrom, len(got), n_want), flush=True)

    spec = importlib.util.spec_from_file_location(
        "pilot09", str(ROOT / "09_crossdomain_pilot.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    clump = mod.load_clump_module()

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

        for key in ENDPOINTS:
            ogp = DATA / ("og_%s_%s.tsv" % (key, gcst))
            if not ogp.exists():
                continue
            obs = {r["rsid"]: r for r in read_rows(ogp)}
            usable = []
            for rs, r in exp.items():
                o = obs.get(rs)
                if not o:
                    continue
                try:
                    b_out = float(o["beta"]); se_out = float(o["se"])
                    eaf = float(o["eaf"]) if o["eaf"] not in ("", "NA") else None
                    bx = float(r["beta"])
                except (TypeError, ValueError):
                    continue
                ba, st = align(r["ea"], r["oa"], o["nea"], o["ea"], b_out, eaf)
                if ba is None or st not in ("ok", "flip") or se_out <= 0 or bx == 0:
                    continue
                usable.append({"rsid": rs, "chr": int(r["chr"]),
                               "pos": int(r["pos37"]), "p": float(r["p"]),
                               "b": ba / bx, "se": se_out / abs(bx),
                               "F": (bx / float(r["se"])) ** 2})
                var_rows.append([gcst, key, rs, r["p"], ba / bx, se_out / abs(bx), st])
            if len(usable) < 2:
                continue
            b_u, se_u = ivw([(u["b"], u["se"]) for u in usable])
            if b_u is None:
                continue
            cands = [(u["rsid"], u["chr"], u["pos"], u["p"]) for u in usable]
            kept, _ = clump.clump_taxon(cands, geno, R2_MAX, KB_CLUMP)
            kept_rs = {k[0] for k in kept}
            kept_u = [u for u in usable if u["rsid"] in kept_rs]
            if not kept_u:
                continue
            b_c, se_c = ivw([(u["b"], u["se"]) for u in kept_u])
            k, kk = len(usable), len(kept_u)
            pair_rows.append([name, gcst, key, k, kk, round(se_u, 6), round(se_c, 6),
                              round((se_c / se_u) ** 2, 4),
                              round((se_c / se_u) / math.sqrt(k), 4),
                              round(best_share, 3), round(b_u, 4), round(b_c, 4),
                              round(median([u["F"] for u in usable]), 2),
                              sum(1 for u in usable if u["F"] < 10)])
    write_rows(OUT / "pairs_opengwas.csv",
               ["trait", "gcst", "endpoint", "k", "k_kept", "se_unclumped",
                "se_clumped", "k_eff", "f", "max_window_share", "b_unclumped",
                "b_clumped", "median_F", "n_F_lt10"], pair_rows)
    print("[og]     %d -> %s" % (len(pair_rows), OUT / "pairs_opengwas.csv"))


def stats_step():
    rows = [r for r in read_rows(OUT / "pairs_opengwas.csv") if int(r["k"]) >= 2]
    if not rows:
        print("[og]    ")
        return

    def rank(v):
        o = sorted(range(len(v)), key=lambda i: v[i])
        out = [0.0] * len(v)
        i = 0
        while i < len(o):
            j = i
            while j + 1 < len(o) and v[o[j + 1]] == v[o[i]]:
                j += 1
            a = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                out[o[k]] = a
            i = j + 1
        return out

    def sp(x, y):
        rx, ry = rank(x), rank(y)
        n = len(x)
        mx, my = sum(rx) / n, sum(ry) / n
        num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
        dx = (sum((a - mx) ** 2 for a in rx)) ** 0.5
        dy = (sum((b - my) ** 2 for b in ry)) ** 0.5
        return num / (dx * dy) if dx and dy else float("nan")

    f = [float(r["f"]) for r in rows]
    share = [float(r["max_window_share"]) for r in rows]
    p1 = sp(share, f)
    p4 = sp([float(r["k_kept"]) for r in rows], f)
    lines = ["# OpenGWAS          FinnGen ", "",
             "    k>=2  %d f    %.3f " % (len(rows), median(f)),
             "P1    ↔f %.3f P4 block  ↔f %.3f " % (p1, p4), ""]
    txt = "\n".join(lines)
    (OUT / "stats_opengwas.md").write_text(txt, encoding="utf-8", newline="")
    print(txt)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--step", required=True, choices=["fetch", "mr", "stats"])
    a = ap.parse_args()
    if a.step == "fetch":
        fetch_step()
    elif a.step == "mr":
        mr_step()
    else:
        stats_step()


if __name__ == "__main__":
    main()
