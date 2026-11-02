# -*- coding: utf-8 -*-
"""Render the manuscript markdown to docx (single source)."""
import re
import sys
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

from docx import Document
from docx.shared import Pt, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT

MD = Path("D:/Horzontal_project/20B3/   /20B3_companion_manuscript_v0.3.md")
FIG = Path("D:/Horzontal_project/20B3/work/crossdomain/output/fig1_crossdomain.png")
OUT = Path("D:/Horzontal_project/20B3/work/crossdomain/output/20B3_companion_manuscript_v0.3.docx")


def add_runs(par, text, size=11, bold=False):
    """   **bold**   `code`         """
    text = text.replace("`", "")
    parts = re.split(r"(\*\*[^*]+\*\*)", text)
    for part in parts:
        if not part:
            continue
        if part.startswith("**") and part.endswith("**"):
            r = par.add_run(part[2:-2])
            r.bold = True
        else:
            r = par.add_run(part)
            r.bold = bold
        r.font.size = Pt(size)


def flush_table(doc, block):
    rows = [r for r in block if not re.match(r"^\|[\s\-:|]+\|$", r)]
    cells = [[c.strip() for c in r.strip().strip("|").split("|")] for r in rows]
    ncol = max(len(r) for r in cells)
    table = doc.add_table(rows=len(cells), cols=ncol)
    table.style = "Light Grid Accent 1"
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    for i, row in enumerate(cells):
        for j in range(ncol):
            txt = row[j] if j < len(row) else ""
            cell = table.cell(i, j)
            cell.text = ""
            p = cell.paragraphs[0]
            add_runs(p, txt, size=8.5, bold=(i == 0))
    doc.add_paragraph()


def main():
    lines = MD.read_text(encoding="utf-8").splitlines()
    doc = Document()
    style = doc.styles["Normal"]
    style.font.name = "Times New Roman"
    style.font.size = Pt(11)
    for s in doc.sections:
        s.left_margin = s.right_margin = Inches(1.0)

    i = 0
    first_h1 = True
    while i < len(lines):
        line = lines[i].rstrip()
        if not line.strip():
            i += 1
            continue
        if line.startswith("|"):
            block = []
            while i < len(lines) and lines[i].startswith("|"):
                block.append(lines[i])
                i += 1
            flush_table(doc, block)
            continue
        if line.startswith("# "):
            if first_h1:
                p = doc.add_paragraph()
                p.alignment = WD_ALIGN_PARAGRAPH.CENTER
                r = p.add_run(line[2:].strip())
                r.bold = True
                r.font.size = Pt(16)
                first_h1 = False
            else:
                doc.add_heading(line[2:].strip(), level=1)
            i += 1
            continue
        if line.startswith("## "):
            doc.add_heading(line[3:].strip(), level=1)
            i += 1
            continue
        if line.startswith("### "):
            doc.add_heading(line[4:].strip(), level=2)
            i += 1
            continue
        if line.startswith("- "):
            p = doc.add_paragraph(style="List Bullet")
            add_runs(p, line[2:].strip(), size=10.5)
            i += 1
            continue
        if line.strip() == "---":
            i += 1
            continue
        p = doc.add_paragraph()
        add_runs(p, line.strip(), size=11)
        p.paragraph_format.space_after = Pt(6)

        if line.strip().startswith("**Figure 1.**") and FIG.exists():
            pic = doc.add_paragraph()
            pic.alignment = WD_ALIGN_PARAGRAPH.CENTER
            pic.add_run().add_picture(str(FIG), width=Inches(6.0))
            cap = doc.add_paragraph()
            cap.alignment = WD_ALIGN_PARAGRAPH.CENTER
            r = cap.add_run("Figure 1. Attained fraction f against the share of variants in "
                            "the largest 10 Mb window, coloured by the number of blocks "
                            "retained after clumping (r²≤0.001). Dashed references: "
                            "microbiome lactase axis (0.583) and non-lactase axis (0.996).")
            r.italic = True
            r.font.size = Pt(9)
        i += 1

    OUT.parent.mkdir(parents=True, exist_ok=True)
    doc.save(OUT)
    print("[docx] -> %s %d   / %d   / %d   "
          % (OUT, len(doc.paragraphs), len(doc.tables), len(doc.inline_shapes)))


if __name__ == "__main__":
    main()
