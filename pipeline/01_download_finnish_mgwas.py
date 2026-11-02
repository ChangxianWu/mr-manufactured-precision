# -*- coding: utf-8 -*-
"""Download the Finnish single-cohort microbiome GWAS summary statistics (GWAS Catalog, resumable chunked HTTP)."""
import argparse
import csv
import hashlib
import os
import re
import sys
import threading
import time
from pathlib import Path
from urllib.parse import urlparse

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

BASE = os.path.realpath("D:/Horzontal_project")
sys.path.insert(0, BASE)
from cv_shared.net_guard import fetch, open_url  # noqa: E402

ROOT = Path(BASE) / "20B3"
LIBDIR = ROOT / "work" / "data" / "exposure_libraries"
DST = ROOT / "work" / "data" / "finnish"
MAN = ROOT / "work" / "data" / "_finnish_manifest.tsv"
DST.mkdir(parents=True, exist_ok=True)

PLAN = LIBDIR / "finnish_replication_plan.tsv"
_lock = threading.Lock()


def md5_of(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


_MD5_CACHE = {}


def expected_md5(url_dir, fname):
    """GWAS Catalog    md5sum.txt          md5          """
    key = url_dir.rstrip("/")
    if key not in _MD5_CACHE:
        b = fetch(key + "/md5sum.txt", timeout=60, tries=2)
        _MD5_CACHE[key] = b.decode("utf-8", "replace") if b else ""
    for line in _MD5_CACHE[key].split("\n"):
        parts = line.split()
        if len(parts) >= 2 and parts[1] == fname:
            return parts[0]
    return None


def httpsify(u):
    """GWAS Catalog       sumstats URL   http://       net_guard     

           https    `assert_safe_url`     https   EBI FTP   https
        Stage 0     200         https              
    """
    if u.startswith("http://"):
        return "https://" + u[len("http://"):]
    return u


def one(acc, url_base, fname):
    """     Finnish sumstats    (acc, status, size, md5, note) """
    url_base = httpsify(url_base)
    url = url_base.rstrip("/") + "/" + fname
    dst = DST / fname
    part = DST / (fname + ".part")

    if dst.exists() and dst.stat().st_size > 1_000_000:
        exp = expected_md5(url_base, fname)
        got = md5_of(dst)
        if exp and got != exp:
            with _lock:
                print("  [md5     ] %s got=%s exp=%s" % (acc, got[:8], exp[:8]))
            dst.unlink()
        else:
            with _lock:
                print("  [  ] %s (%d B, md5=%s)" % (acc, dst.stat().st_size, got[:8]))
            return (acc, "skip_md5ok", dst.stat().st_size, got, "")

    body = fetch(url, timeout=900, tries=3, backoff=(5, 15, 30))
    if body is None:
        return (acc, "fail_download", 0, "", url)

    part.write_bytes(body)
    got = md5_of(part)
    exp = expected_md5(url_base, fname)
    if exp and got != exp:
        part.unlink()
        return (acc, "fail_md5", len(body), got,
                "expected=%s url=%s" % (exp, url))
    part.replace(dst)
    return (acc, "ok", dst.stat().st_size, got, "")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()

    plan = list(csv.DictReader(open(PLAN, encoding="utf-8"), delimiter="\t"))
    accs = {}
    for r in plan:
        u = httpsify(r["sumstats_url"].strip()).rstrip("/")
        acc = r["finnish_accession"].strip()
        if u.endswith("/" + acc):
            accs.setdefault(acc, u)
        else:
            accs.setdefault(acc, u + "/" + acc)
    items = sorted(accs.items())
    if a.limit:
        items = items[:a.limit]
    print("[01]      %d   Finnish sumstats workers=%d " % (len(items), a.workers))

    rows, lock = [], threading.Lock()
    sem = threading.Semaphore(a.workers)
    t0 = time.time()

    def work(acc, d):
        with sem:
            fn = "%s_buildGRCh38.tsv.gz" % acc
            try:
                res = one(acc, d, fn)
            except Exception as exc:  # noqa: BLE001
                res = (acc, "fail_exc", 0, "", "%s: %s" % (type(exc).__name__, exc))
            with lock:
                rows.append(res)
                print("  %-14s %-16s %8d B  %.0fs" % (res[0], res[1], res[2], time.time() - t0))

    ths = [threading.Thread(target=work, args=(acc, d)) for acc, d in items]
    for t in ths:
        t.start()
    for t in ths:
        t.join()

    ok = sum(1 for r in rows if r[1] in ("ok", "skip_md5ok"))
    with open(MAN, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["accession", "status", "size", "md5", "note"])
        for r in sorted(rows):
            w.writerow(r)
    print("\n[01]    %d/%d -> %s" % (ok, len(items), MAN))
    print("[01]    %.1f         %.2f GB" % (
        (time.time() - t0) / 60,
        sum(f.stat().st_size for f in DST.glob("*.gz")) / 1073741824))


if __name__ == "__main__":
    main()
