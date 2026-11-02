# -*- coding: utf-8 -*-
"""Journal candidate probing against Europe PMC (genre and audience counts)."""
import csv
import sys
import time
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

sys.path.insert(0, "D:/Horzontal_project")
from cv_shared.probe_kit import epmc_count, phrase, title_abs  # noqa: E402


def journal_exact(full_name):
    """          `"name"[Journal]`   /           "Genetic Epidemiology"[Journal]=3967
      JOURNAL:"Genetic Epidemiology"=101    [[probe-epmc-parenthesis-or-trap]] 
    """
    return 'JOURNAL:"%s"' % full_name

DATE = "2026-10-12"
OUT = Path(r"D:/Horzontal_project/20B3/work/review") / ("journal_probe_%s.tsv" % DATE)

CANDIDATES = [
    "BMC Medical Research Methodology",
    "BMC Infectious Diseases",
    "BMC Bioinformatics",
    "Frontiers in Microbiology",
    "Frontiers in Cellular and Infection Microbiology",
    "Frontiers in Genetics",
    "Frontiers in Immunology",
    "Microbial Genomics",
    "Microbiology Spectrum",
    "Scientific Reports",
    "PLOS ONE",
    "Medicine",
    "Microorganisms",
    "Journal of Translational Medicine",
    "International Journal of Epidemiology",
    "European Journal of Epidemiology",
    "Genetic Epidemiology",
    "Human Genetics",
    "Gut Microbes",
    "mSystems",
]

TOPICS = {
    "H1_microbiome_MR": (
        title_abs("gut microbiome", "gut microbiota", "microbiome", "microbiota")
        + " AND " + title_abs("Mendelian randomization", "Mendelian randomisation")
    ),
    "H2_microbiome_infection": (
        title_abs("gut microbiome", "gut microbiota", "microbiome", "microbiota")
        + " AND "
        + title_abs("infection", "infectious disease", "sepsis", "pneumonia",
                    "gastroenteritis", "diarrhoea", "diarrhea")
    ),
    "H3_ld_clumping_methods": (
        title_abs("LD clumping", "linkage disequilibrium clumping", "instrument independence",
                  "instrument selection", "correlated instruments")
        + " AND " + title_abs("Mendelian randomization", "Mendelian randomisation")
    ),
    "H4_reproducibility_audit": (
        title_abs("reproducibility", "reproducibility audit", "instrument structure",
                  "instrument inflation", "type 1 error inflation", "winner's curse")
        + " AND " + title_abs("Mendelian randomization", "Mendelian randomisation",
                              "genetic epidemiology", "instrumental variable")
    ),
}

Y3 = 'FIRST_PDATE:[2024-01-01 TO 2026-12-31]'

rows = []
print("[11]        %d   × %d    +   " % (len(CANDIDATES), len(TOPICS)))

for j in CANDIDATES:
    jf = journal_exact(j)
    total = epmc_count(jf + " AND " + Y3)
    row = {"journal": j, "total_3y": total,
           "volume_gate": "PASS" if (total or 0) >= 150 else "FAIL"}
    print("  %-52s 3y=%6s %s" % (j, total, row["volume_gate"]))
    for tname, tq in TOPICS.items():
        n = epmc_count(jf + " AND " + Y3 + " AND (" + tq + ")")
        row[tname] = n
        print("      %-26s %6s" % (tname, n))
        time.sleep(0.35)
    rows.append(row)
    time.sleep(0.35)

OUT.parent.mkdir(parents=True, exist_ok=True)
cols = ["journal", "total_3y", "volume_gate"] + list(TOPICS)
with open(OUT, "w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
    w.writeheader()
    w.writerows(rows)
print("\n[11] -> %s" % OUT)
