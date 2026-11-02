# -*- coding: utf-8 -*-
"""Verify every reference against Europe PMC and write the reference pool."""
import csv
import io
import json
import sys
import time
import urllib.parse
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

sys.path.insert(0, "D:/Horzontal_project")
from cv_shared.probe_kit import _get  # noqa: E402

DATE = "2026-10-11"
TOOLS = Path(r"D:/Horzontal_project/20B3/work/tools")
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

REFS = [
    ("Kurilshikov2021", "33462485", "MiBioGen      211 taxa / 24    meta ",
     "host factors influencing human gut microbiome composition"),
    ("Qin2022", "35115689", "Finnish       GWAS 473 taxa       ",
     "Combined effects of host genetics and diet on human gut microbiota"),
    ("LoperaMaya2022", "35115690", "Dutch Microbiome Project 207 taxa        ",
     "host genetics on the gut microbiome in 7,738 participants"),
    ("Kurki2023", "36653562", "FinnGen R12         ",
     "FinnGen provides genetic insights"),
    ("DaveySmith2014", "15075143", "MR         -  -      ",
     "Mendelian randomization: prospects, potentials, and limitations"),
    ("Bowden2015", "26050253", "MR-Egger             ",
     "Mendelian randomization with invalid instruments"),
    ("Bowden2016", "27061298", "            ",
     "Consistent Estimation in Mendelian Randomization"),
    ("Verbanck2018", "29686387", "MR-PRESSO        ",
     "Detection of widespread horizontal pleiotropy"),
    ("Hemani2018", "29846171", "TwoSampleMR    ",
     "The MR-Base platform"),
    ("Burgess2015", "24114802", "    IVW      ",
     "Mendelian randomization analysis with multiple genetic variants"),
    ("Skrivankova2021", "34698778", "STROBE-MR     ",
     "Strengthening the Reporting of Observational Studies"),
    ("Hemani2017", "29149188", "Steiger      ",
     "Orienting the causal relationship"),
    ("Giambartolomei2014", "24830394", "coloc       PP.H4 ",
     "Bayesian test for colocalisation"),
    ("Auton2015", "26432245", "1000 Genomes phase 3 LD      ",
     "A global reference for human genetic variation"),
    ("Blekhman2015", "26374288", "     -       LCT   ",
     "Host genetic variation impacts microbiome composition"),
    ("Wacklin2011", "21625510", "FUT2            ",
     "Secretor genotype (FUT2 gene) is strongly associated"),
    ("Sanna2019", "30778224", "  -   -   MR SCFA   ",
     "Causal relationships among the gut microbiome"),
    ("Rautanen2015", "25533491", "      GWAS      ",
     "survival from sepsis due to pneumonia"),
    ("Singer2016", "26903338", "Sepsis-3   ",
     "The Third International Consensus Definitions for Sepsis"),
]

FIELDS = ["key", "pmid", "use", "fingerprint", "title", "journal", "year",
          "volume", "issue", "pages", "authors", "doi"]


def main():
    TOOLS.mkdir(parents=True, exist_ok=True)
    rows, evidence = [], []
    for key, pmid, use, fp in REFS:
        u = EPMC + "?" + urllib.parse.urlencode(
            {"query": "EXT_ID:%s" % pmid, "format": "json",
             "resultType": "core", "pageSize": "1"})
        rec = None
        for _ in range(3):
            body = _get(u, tries=2)
            if body:
                rl = json.loads(body).get("resultList", {}).get("result", [])
                if rl:
                    rec = rl[0]
                    break
            time.sleep(1.5)
        if rec is None:
            print("[06] %-18s PMID %s     " % (key, pmid))
            evidence.append({"key": key, "pmid": pmid, "status": "fetch_failed"})
            continue
        title = rec.get("title", "") or ""
        ok = fp.lower() in title.lower()
        jin = rec.get("journalInfo") or {}
        evidence.append({
            "key": key, "pmid": pmid, "status": "ok" if ok else "title_mismatch",
            "title": title, "fingerprint": fp,
            "journal": (jin.get("journal") or {}).get("title", ""),
            "year": rec.get("pubYear", ""), "volume": jin.get("volume", ""),
            "issue": jin.get("issue", ""), "pages": rec.get("pageInfo", ""),
            "doi": rec.get("doi", ""), "authors": rec.get("authorString", ""),
        })
        if ok:
            rows.append({
                "key": key, "pmid": pmid, "use": use, "fingerprint": fp,
                "title": title,
                "journal": (jin.get("journal") or {}).get("title", ""),
                "year": rec.get("pubYear", ""), "volume": jin.get("volume", ""),
                "issue": jin.get("issue", ""), "pages": rec.get("pageInfo", ""),
                "authors": rec.get("authorString", ""), "doi": rec.get("doi", ""),
            })
            print("[06] %-18s PMID %-9s OK  %s" % (key, pmid, title[:66]))
        else:
            print("[06] %-18s PMID %-9s        =%s |    =%s"
                  % (key, pmid, title[:60], fp))
        time.sleep(0.35)

    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=FIELDS, delimiter="\t")
    w.writeheader()
    w.writerows(rows)
    (TOOLS / ("refs_pool_%s.tsv" % DATE)).write_text(
        buf.getvalue(), encoding="utf-8-sig", newline="")
    (TOOLS / ("refs_verify_%s.json" % DATE)).write_text(
        json.dumps(evidence, ensure_ascii=False, indent=1),
        encoding="utf-8", newline="")

    print("\n[06]    %d /    %d" % (len(rows), len(REFS)))
    for e in evidence:
        if e["status"] != "ok":
            print("   ✗ %-18s %s" % (e["key"], e["status"]))


if __name__ == "__main__":
    main()
