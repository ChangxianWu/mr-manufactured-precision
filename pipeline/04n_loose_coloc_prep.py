# -*- coding: utf-8 -*-
"""Prepare colocalization inputs for the relaxed-threshold sets."""
import argparse
import io
import sys
import urllib.parse
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

BASE = Path("D:/Horzontal_project")
sys.path.insert(0, str(BASE))
from cv_shared.bgzf_tabix import TabixReader  # noqa: E402
from cv_shared.net_guard import open_url  # noqa: E402

FG_BASE = "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release"
CAT = ("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
       "GCST90017001-GCST90018000/GCST90017038/GCST90017038_buildGRCh37.tsv")
CO = BASE / "20B3" / "work" / "data" / "coloc"
CACHE = BASE / "20B3" / "work" / "cache" / "tabix"

MAP38 = BASE / "20B3" / "work" / "data" / "lead_snp_hg38_map.csv"


def pos38_of(rsid):
    import csv as _csv
    with open(MAP38, encoding="utf-8") as fh:
        for r in _csv.DictReader(fh, delimiter="\t"):
            if r["rsid"] == rsid:
                return int(r["chrom"]), int(r["pos38"])
    raise SystemExit("[prep] %s    %s        " % (rsid, MAP38))


PAD = 100_000
FF_HEADER = ["#chrom", "pos", "ref", "alt", "rsids", "nearest_genes", "pval", "mlogp",
             "beta", "sebeta", "af_alt", "af_alt_cases", "af_alt_controls"]


def check_host(url):
    """     https +         localhost/  /  /      """
    u = urllib.parse.urlparse(url)
    if u.scheme != "https":
        raise SystemExit("[prep]     https %s" % url)
    host = (u.hostname or "").lower()
    if host in ("localhost", "127.0.0.1", "::1") or host.endswith(".localhost"):
        raise SystemExit("[prep]        %s" % host)
    import ipaddress
    try:
        ip = ipaddress.ip_address(host)
    except ValueError:
        return host
    if ip.is_private or ip.is_loopback or ip.is_reserved or ip.is_link_local:
        raise SystemExit("[prep]     /     %s" % host)
    return host


def extract_finngen(tag, pheno, chrom, pos, out_name):
    check_host(FG_BASE)
    url = "%s/finngen_R12_%s.gz" % (FG_BASE, pheno)
    tr = TabixReader(url, cache_dir=str(CACHE), name="finngen_R12_%s.gz" % pheno)
    _hdr, rows = tr.query(str(chrom), pos - PAD, pos + PAD, header=FF_HEADER)
    ncol = len(FF_HEADER)
    i_rs = FF_HEADER.index("rsids")
    good = [r for r in rows if len(r) >= ncol]
    n_bad = len(rows) - len(good)
    n_nors = sum(1 for r in good if not r[i_rs].strip())
    path = CO / out_name
    with io.open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(FF_HEADER) + "\n")
        for r in good:
            fh.write("\t".join(r[:ncol]) + "\n")
    print("[prep] %s %d       %d     rsID %d -> %s"
          % (tag, len(good), n_bad, n_nors, path.name))


def extract_exposure(chrom, pos, out_name, hi=88800000, lo=86200000, tries=3):
    """       TSV GRCh37                    p    

           coloc   **      **       top-hit/p<0.05     
                  ——                      

         EBI      TSV      EOF        ~10%      
    10%                    **    **    Range   
              Range              Content-Length    complete 
    """
    import time
    check_host(CAT)
    path = CO / out_name
    hdr = None
    n_total = n_chr = n_pgt05 = n_win = 0
    maxp_win = 0.0
    keep = []
    complete = False
    nbytes = total_len = 0
    for attempt in range(tries):
        n_total = n_chr = n_pgt05 = 0
        maxp_win = 0.0
        keep = []
        nbytes = 0
        try:
            with open_url(CAT, timeout=600) as resp:
                total_len = int(resp.headers.get("Content-Length") or 0)
                raw_h = resp.readline()
                nbytes += len(raw_h)
                hdr = raw_h.decode("utf-8", "replace").rstrip("\n").split("\t")
                idx_c = hdr.index("chromosome")
                idx_p = hdr.index("base_pair_location")
                idx_pv = hdr.index("p_value")
                idx_n = hdr.index("n")
                for raw in resp:
                    nbytes += len(raw)
                    line = raw.decode("utf-8", "replace")
                    f = line.rstrip("\n").split("\t")
                    if len(f) <= max(idx_c, idx_p, idx_pv):
                        continue
                    n_total += 1
                    try:
                        pv = float(f[idx_pv])
                    except ValueError:
                        pv = None
                    if pv is not None and pv > 0.05:
                        n_pgt05 += 1
                    c = f[idx_c]
                    if c not in (str(chrom), "chr%d" % chrom):
                        continue
                    n_chr += 1
                    try:
                        p = int(float(f[idx_p]))
                    except ValueError:
                        continue
                    if lo <= p <= hi:
                        keep.append(line)
                        if pv is not None:
                            maxp_win = max(maxp_win, pv)
            if total_len and nbytes >= total_len:
                complete = True
                break
            sys.stderr.write("[prep]    %d/%d     EOF        %d   \n"
                             % (nbytes, total_len, attempt + 1))
        except Exception as exc:  # noqa: BLE001
            sys.stderr.write("[prep]        %d    %s\n" % (attempt + 1, str(exc)[:120]))
        time.sleep(5)
    n_win = len(keep)
    with io.open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(hdr) + "\n")
        for line in keep:
            if not line.endswith("\n"):
                line += "\n"
            fh.write(line)
    audit = CO / "loose_exposure_audit.csv"
    with io.open(audit, "w", encoding="utf-8", newline="") as fh:
        fh.write("source\tbytes_read\ttotal_bytes\tcomplete\trows_total\trows_chr%d\t"
                 "rows_window\tmax_p_window\trows_p_gt_0p05_total\tn_col\n" % chrom)
        fh.write("%s\t%d\t%d\t%d\t%s\t%d\t%d\t%.6g\t%s\t%s\n"
                 % (CAT.rsplit("/", 1)[-1], nbytes, total_len, int(complete),
                    n_total if complete else "", n_chr, n_win, maxp_win,
                    n_pgt05 if complete else "",
                    keep[0].split("\t")[idx_n].strip() if n_win else ""))
    print("[prep]    chr9:%d-%d    %d     %d/%d      =%s chr%d %d   -> %s"
          % (lo, hi, n_win, nbytes, total_len, complete, chrom, n_chr, path.name))
    if not complete:
        print("[prep] ⚠️      ⇒   **  **  '    p>0.05   '   ")
    elif n_pgt05 == 0:
        print("[prep] ⚠️    %d     p>0.05    ⇒      p<0.05     "
              "   coloc             " % n_total)
    return dict(total_rows=(n_total if complete else None), rows_chr=n_chr,
                rows_window=n_win, max_p_window=maxp_win,
                rows_p_gt_0p05=(n_pgt05 if complete else None), complete=complete)



def pos37_of(rsid):
    """GRCh37      MiBioGen                       """
    import csv as _csv
    p_instr = BASE / "20B3" / "work" / "output" / "instruments_primary.csv"
    with open(p_instr, encoding="utf-8") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > 6 and f[4] == rsid:
                return int(f[5]), int(f[6])
    raise SystemExit("[prep] %s    instruments_primary.csv  " % rsid)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-exposure", action="store_true")
    a = ap.parse_args()
    CO.mkdir(parents=True, exist_ok=True)
    c1, p1 = pos38_of("rs35866622")
    c2, p2 = pos38_of("rs736744")
    _c2, e2 = pos37_of("rs736744")
    print("[prep] rs35866622 GRCh38 %d:%d | rs736744 GRCh38 %d:%d / GRCh37 %d:%d"
          % (c1, p1, c2, p2, _c2, e2))
    extract_finngen("PNEUMOBACT×rs35866622", "J10_PNEUMOBACT", c1, p1,
                    "outcome_PNEUMOBACT_chr19_region.tsv")
    extract_finngen("BACT_OTHER×rs736744", "AB1_OTHER_BACTERIAL", c2, p2,
                    "outcome_BACT_OTHER_chr9_region.tsv")
    if not a.skip_exposure:
        extract_exposure(_c2, e2, "GCST90017038_chr9_region.tsv",
                         lo=e2 - 250_000, hi=e2 + 250_000)


if __name__ == "__main__":
    main()
