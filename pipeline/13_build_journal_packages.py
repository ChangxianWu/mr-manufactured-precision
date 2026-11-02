# -*- coding: utf-8 -*-
"""Build clean per-journal submission packages."""
import json
import shutil
import sys
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = Path(r"D:/Horzontal_project/20B3")
PKG = ROOT / "     "
JP = ROOT / "   "
FIGSRC = PKG / "03_       " / "figures"
SUPP = PKG / "02_     "
MANUSCRIPT = "20B3_Manuscript_v0.2.docx"
TITLE = ("LD clumping removes almost all apparent gut microbiome-infection Mendelian "
         "randomization signals: manufactured precision follows a perfect-redundancy bound")

FIG_ORDER = ["fig1_funnel.png", "fig2_volcano.png", "fig3_forest.png",
             "fig4_law.png", "fig5_replication.png"]

JOURNALS = [
    {
        "name": "International_Journal_of_Epidemiology",
        "tag": "IJE",
        "note": ("T1 primary. The journal publishes bias, assumption and reporting-guidance papers on "
                 "Mendelian randomization, which is the genre of this manuscript. Recent examples "
                 "include a structured quality and risk-of-bias tool for MR studies and a paper on "
                 "exposure misidentification as a recurring pitfall."),
        "limits": ("Confirm the current article type, word limit and figure limit on the journal's "
                   "author-information page before uploading; the journal's site blocked automated "
                   "retrieval, so limits could not be measured here."),
        "keywords": ("Mendelian randomization; instrument independence; linkage disequilibrium; "
                     "reproducibility; bias; gut microbiome"),
    },
    {
        "name": "European_Journal_of_Epidemiology",
        "tag": "EJE",
        "note": ("T2 first alternative. The journal recently published a guide to risk-of-bias "
                 "assessment for Mendelian randomisation studies, so reporting-standards work is "
                 "within its scope."),
        "limits": ("Confirm the current article type and limits on the journal's author-information "
                   "page before uploading."),
        "keywords": ("Mendelian randomization; reproducibility; instrument selection; "
                     "linkage disequilibrium; infectious disease"),
    },
]

README_TMPL = """# Submission package

**Title**: {title}

## Contents (numbered directory = upload slot)
- `01_Manuscript/{manuscript}` - main text with {nfig} embedded figures, three main-text tables and
  {nref} references; structured abstract of {nabsw} words
- `02_Figures/Figure1..{nfig}.png` - 300 dpi, numbered in order of first citation
- `03_Supplementary/Supplementary_Tables_20B3.xlsx` - supplementary tables S1-S12 (12 sheets)
- `04_Cover_Letter/Cover_Letter.docx`
- `05_Reporting_Checklist/STROBE_MR_checklist.docx` - STROBE-MR items, one row pending author metadata

**Keywords** (enter in the submission system): {keywords}

## Journal-specific notes
{note}

Limits currently observed: {limits}

## Author actions required before upload
1. Complete the author block, affiliations, corresponding author, email and ORCID in the
   submission system.
2. Complete the reporting-checklist row marked PENDING (funding and competing interests) and the
   author-contributions placeholder in the manuscript Declarations section.
3. Enter the keywords listed above.
4. Confirm the article type, current limits and article-processing charge on the journal website.
5. **Record the code deposition route in the decision letter.** The manuscript states that analysis
   scripts, intermediate result tables and the keyed tabulation of reported values are
   self-archived with the authors and available from the corresponding author on request.
   Public deposition is strongly encouraged but is not a condition of acceptance; if a public
   archive with a DOI is not possible, deposit at least the LD clumping implementation, because
   its correctness cannot be verified from the supplementary tables.

## Reproducibility
Manuscript values are generated from a single frozen keyed table ({nkeys} keys) produced by the archived
analysis scripts; an independent recomputation script re-derives 126 reported quantities from
the raw result files, and a numeric-difference check confirms that writing-stage edits changed
no numeric value.
"""


def abstract_words(docx_path):
    """                            """
    from docx import Document
    texts = [p.text.strip() for p in Document(str(docx_path)).paragraphs]
    i = next(i for i, t in enumerate(texts) if t.lower().startswith("abstract"))
    j = next(i for i, t in enumerate(texts) if t.startswith("1. Introduction"))
    return len(" ".join(texts[i + 1:j]).split())


def canonical_keys():
    """canonical                   """
    rows = (ROOT / "work/output/canonical_numbers.tsv").read_text(encoding="utf-8").splitlines()
    return sum(1 for r in rows[1:] if r.strip())


def main():
    if not FIGSRC.exists():
        raise SystemExit("    figures    %s    make_submission_package.py " % FIGSRC)
    for f in FIG_ORDER:
        if not (FIGSRC / f).exists():
            raise SystemExit("    %s     04_sensitivity.R       " % f)
    if JP.exists():
        shutil.rmtree(JP)

    counts = {}
    ms_src = PKG / "01_     " / MANUSCRIPT
    nabsw = abstract_words(ms_src)
    nkeys = canonical_keys()
    for j in JOURNALS:
        base = JP / j["name"]
        for sub in ("01_Manuscript", "02_Figures", "03_Supplementary",
                    "04_Cover_Letter", "05_Reporting_Checklist"):
            (base / sub).mkdir(parents=True, exist_ok=True)

        shutil.copy2(ms_src, base / "01_Manuscript" / MANUSCRIPT)

        for i, f in enumerate(FIG_ORDER, 1):
            shutil.copy2(FIGSRC / f, base / "02_Figures" / ("Figure%d.png" % i))

        shutil.copy2(SUPP / "Supplementary_Tables_20B3.xlsx",
                     base / "03_Supplementary" / "Supplementary_Tables_20B3.xlsx")
        shutil.copy2(SUPP / "Cover_Letter.docx",
                     base / "04_Cover_Letter" / "Cover_Letter.docx")
        shutil.copy2(SUPP / "STROBE_MR_checklist.docx",
                     base / "05_Reporting_Checklist" / "STROBE_MR_checklist.docx")

        (base / "README_submission_checklist.md").write_text(
            README_TMPL.format(title=TITLE, manuscript=MANUSCRIPT,
                               nfig=len(FIG_ORDER), nref=19, nabsw=nabsw, nkeys=nkeys,
                               keywords=j["keywords"], note=j["note"], limits=j["limits"]),
            encoding="utf-8", newline="")
        counts[j["name"]] = sum(1 for p in base.rglob("*") if p.is_file())
        print("[13]    %-16s %d   " % (j["name"], counts[j["name"]]))

    (JP / "package_checklist.json").write_text(json.dumps(
        {"packages": [
            {"package": j["name"],
             "manuscript_glob": "01_Manuscript/%s" % MANUSCRIPT,
             "figures_glob": "02_Figures/Figure*.png",
             "required_files": [
                 "README_submission_checklist.md",
                 "03_Supplementary/Supplementary_Tables_20B3.xlsx",
                 "04_Cover_Letter/Cover_Letter.docx",
                 "05_Reporting_Checklist/STROBE_MR_checklist.docx",
             ]}
            for j in JOURNALS]}, ensure_ascii=False, indent=2), encoding="utf-8", newline="")
    total = sum(counts.values())
    print("[13] package_checklist.json     %d    /   %d    " % (len(JOURNALS), total))


if __name__ == "__main__":
    main()
