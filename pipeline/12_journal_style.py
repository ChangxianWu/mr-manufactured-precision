# -*- coding: utf-8 -*-
"""Journal style comparison (structured abstract, word limits)."""
import csv
import json
import sys
import time
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

sys.path.insert(0, "D:/Horzontal_project")
from cv_shared.probe_kit import epmc_search  # noqa: E402

DATE = "2026-10-12"
OUT = Path(r"D:/Horzontal_project/20B3/work/review") / ("journal_style_%s.tsv" % DATE)

TARGETS = {
    "Gut Microbes": 'JOURNAL:"Gut Microbes" AND (TITLE_ABS:"Mendelian randomization" OR TITLE_ABS:"Mendelian randomisation") AND FIRST_PDATE:[2025-01-01 TO 2026-12-31]',
    "Frontiers in Cellular and Infection Microbiology": 'JOURNAL:"Frontiers in Cellular and Infection Microbiology" AND (TITLE_ABS:"Mendelian randomization" OR TITLE_ABS:"Mendelian randomisation") AND FIRST_PDATE:[2025-01-01 TO 2026-12-31]',
    "Frontiers in Microbiology": 'JOURNAL:"Frontiers in Microbiology" AND (TITLE_ABS:"Mendelian randomization" OR TITLE_ABS:"Mendelian randomisation") AND FIRST_PDATE:[2025-01-01 TO 2026-12-31]',
    "Microbiology Spectrum": 'JOURNAL:"Microbiology Spectrum" AND FIRST_PDATE:[2025-01-01 TO 2026-12-31]',
    "PLOS ONE": 'JOURNAL:"PLOS ONE" AND (TITLE_ABS:"Mendelian randomization" OR TITLE_ABS:"Mendelian randomisation") AND FIRST_PDATE:[2025-01-01 TO 2026-12-31]',
    "Scientific Reports": 'JOURNAL:"Scientific Reports" AND (TITLE_ABS:"Mendelian randomization" OR TITLE_ABS:"Mendelian randomisation") AND FIRST_PDATE:[2025-01-01 TO 2026-12-31]',
}

rows = []
for j, q in TARGETS.items():
    try:
        r = epmc_search(q, page_size=6)
    except Exception as e:  # noqa: BLE001
        print("  %-52s SEARCH_FAIL %s" % (j, type(e).__name__))
        continue
    res = r.get("results", []) if isinstance(r, dict) else []
    for a in res[:4]:
        ab = (a.get("abstractText") or "").strip()
        secs = []
        for k in ("Background", "Methods", "Results", "Conclusion", "Conclusions",
                  "Introduction", "Objective", "Objective:", "Aim", "Discussion"):
            if k.lower() in ab.lower():
                secs.append(k)
        rows.append({
            "journal": j,
            "pmid": a.get("pmid", ""),
            "year": (a.get("firstPublicationDate") or "")[:4],
            "title_head": (a.get("title") or "")[:70],
            "abstract_words": len(ab.split()),
            "abstract_structured": "yes" if len(set(secs)) >= 2 else "no",
            "abstract_sections": ";".join(dict.fromkeys(secs)),
            "has_keywords": "yes" if a.get("keywordList", {}).get("keyword") else "no",
            "n_authors": len(a.get("authorList", {}).get("author", []) or []),
        })
        print("  %-52s %s %s w=%d struct=%s kw=%s" % (
            j[:52], a.get("pmid", ""), (a.get("firstPublicationDate") or "")[:4],
            len(ab.split()), rows[-1]["abstract_structured"], rows[-1]["has_keywords"]))
    time.sleep(0.4)

OUT.parent.mkdir(parents=True, exist_ok=True)
cols = ["journal", "pmid", "year", "title_head", "abstract_words",
        "abstract_structured", "abstract_sections", "has_keywords", "n_authors"]
with open(OUT, "w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
    w.writeheader()
    w.writerows(rows)
print("\n[12] %d       -> %s" % (len(rows), OUT))
