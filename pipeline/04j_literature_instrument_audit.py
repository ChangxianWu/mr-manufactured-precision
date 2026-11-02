# -*- coding: utf-8 -*-
"""Abstract-level audit of instrument-count reporting in the competing literature."""
import csv
import io
import re
import statistics
import sys
import time
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

sys.path.insert(0, "D:/Horzontal_project")
from cv_shared.probe_kit import epmc_search  # noqa: E402

COMP_TSV = Path("D:/Horzontal_project/20B3/work/tools/b3_stage0_competitors_2026-10-11.tsv")
OUT_TSV = Path("D:/Horzontal_project/20B3/work/review/literature_instrument_audit.tsv")
SUM_TSV = Path("D:/Horzontal_project/20B3/work/review/literature_instrument_audit_summary.tsv")

PATTERNS = [
    re.compile(r"(?P<n>\d{1,4})\s+(?:independent\s+)?(?:SNPs?|single[- ]nucleotide polymorphisms?)"
               r"\s+(?:as|were used as|served as|as instrumental)", re.I),
    re.compile(r"(?P<n>\d{1,4})\s+(?:genetic\s+)?instruments?\b", re.I),
    re.compile(r"(?P<n>\d{1,4})\s+(?:independent\s+)?instrumental variables?\b", re.I),
    re.compile(r"(?:used|using|selected|identified|retained)\s+(?P<n>\d{1,4})\s+"
               r"(?:independent\s+)?(?:SNPs?|variants?)\b", re.I),
    re.compile(r"(?P<n>\d{1,4})\s+(?:independent\s+)?(?:SNPs?|variants?)\s+(?:were|was)\s+"
               r"(?:used|selected|retained|identified|included)", re.I),
]


def sentences(text):
    return re.split(r"(?<=[.;])\s+", text or "")


def extract(abstract):
    """   (n, quoted)   (None,   )           """
    for s in sentences(abstract):
        for pat in PATTERNS:
            m = pat.search(s)
            if m:
                n = int(m.group("n"))
                if 1 <= n <= 5000:
                    return n, s.strip()[:240]
    return None, ""


def read_competitors():
    text = COMP_TSV.read_text(encoding="utf-8-sig")
    head = text.split("\n", 1)[0]
    delim = "\t" if "\t" in head else ","
    return list(csv.DictReader(io.StringIO(text), delimiter=delim))


def main():
    comp = read_competitors()
    seen = set()
    keep = []
    for r in comp:
        pmid = (r.get("pmid") or "").strip()
        title = (r.get("title") or "").strip()
        key = title.lower()
        if not pmid or key in seen:
            continue
        seen.add(key)
        if "skin" in key:
            continue
        keep.append(r)
    print("[04j]    %d     PMID           " % len(keep))

    rows = []
    for r in keep:
        pmid = r["pmid"].strip()
        ab = ""
        try:
            res = epmc_search("EXT_ID:%s" % pmid, page_size=1)
            hits = res.get("results", []) if isinstance(res, dict) else []
            if hits:
                ab = hits[0].get("abstractText") or ""
        except Exception as e:  # noqa: BLE001
            print("   [warn] %s fetch failed: %s" % (pmid, type(e).__name__))
        time.sleep(0.35)
        n, quote = extract(ab)
        rows.append({
            "pmid": pmid,
            "first_pub": r.get("first_pub", ""),
            "journal": (r.get("journal") or "").strip(),
            "title": (r.get("title") or "").strip(),
            "abstract_chars": len(ab),
            "reported_instruments": "" if n is None else n,
            "quoted_evidence": quote,
            "extraction": "extracted" if n is not None else "not_extractable",
        })

    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=list(rows[0].keys()), delimiter="\t",
                       lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    OUT_TSV.write_text(buf.getvalue(), encoding="utf-8", newline="")

    got = [x for x in rows if x["extraction"] == "extracted"]
    nums = sorted(int(x["reported_instruments"]) for x in got)
    summ = [
        ("n_competitor_papers", len(rows)),
        ("n_abstract_text_available", sum(1 for x in rows if x["abstract_chars"] > 0)),
        ("n_extracted", len(got)),
        ("n_not_extractable", len(rows) - len(got)),
        ("extraction_yield_pct", "%.1f" % (100.0 * len(got) / len(rows))),
    ]
    if nums:
        summ += [
            ("reported_instruments_min", nums[0]),
            ("reported_instruments_median", "%.0f" % statistics.median(nums)),
            ("reported_instruments_max", nums[-1]),
            ("n_reporting_ge_20", sum(1 for x in nums if x >= 20)),
            ("n_reporting_ge_100", sum(1 for x in nums if x >= 100)),
        ]
    sbuf = io.StringIO()
    sw = csv.writer(sbuf, delimiter="\t", lineterminator="\n")
    sw.writerow(["metric", "value"])
    sw.writerows(summ)
    SUM_TSV.write_text(sbuf.getvalue(), encoding="utf-8", newline="")
    for k, v in summ:
        print("  %-32s %s" % (k, v))
    print("[04j] -> %s" % OUT_TSV)


if __name__ == "__main__":
    main()
