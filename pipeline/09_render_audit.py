# -*- coding: utf-8 -*-
"""Word-to-PDF render audit: geometry, overlaps, page images."""
import csv
import os
import re
import sys
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = Path(r"D:/Horzontal_project/20B3")
DOCX = ROOT / "work/manuscript/20B3_Manuscript_v0.2.docx"
RENDER = ROOT / "work/render"
OUT_TSV = ROOT / "work/output/render_audit.tsv"
RENDER.mkdir(parents=True, exist_ok=True)

TOL = 2.0            # pt
OVERLAP_FRAC = 0.35


def to_pdf():
    import win32com.client
    pdf = str(DOCX.with_suffix(".pdf"))
    if os.path.exists(pdf):
        os.remove(pdf)
    word = win32com.client.Dispatch("Word.Application")
    word.Visible = False
    try:
        doc = word.Documents.Open(str(DOCX), ReadOnly=True)
        doc.SaveAs(pdf, FileFormat=17)   # 17 = wdFormatPDF
        doc.Close(False)
    finally:
        word.Quit()
    return pdf


def main():
    import fitz  # pymupdf

    pdf = to_pdf()
    d = fitz.open(pdf)
    rows = []
    n_bad = 0
    for i, page in enumerate(d, 1):
        W, H = page.rect.width, page.rect.height
        L, R, T, B = 20.0, W - 20.0, 20.0, H - 20.0
        issues = []

        text_blocks = page.get_text("blocks")
        img_blocks = [b for b in page.get_text("dict")["blocks"] if b.get("type") == 1]
        img_rects = [fitz.Rect(b["bbox"]) for b in img_blocks]

        for b in text_blocks:
            r = fitz.Rect(b[:4])
            txt = (b[4] or "").strip()
            if not txt:
                continue
            if r.x0 < L - TOL or r.x1 > R + TOL or r.y0 < T - TOL or r.y1 > B + TOL:
                issues.append("text out of print area: %r" % txt[:40])
            for ir in img_rects:
                inter = r & ir
                if not inter.is_empty:
                    a = inter.get_area()
                    if a / max(r.get_area(), 1e-6) > OVERLAP_FRAC:
                        issues.append("text overlaps image: %r" % txt[:40])
                        break

        for ir in img_rects:
            if ir.x0 < L - TOL or ir.x1 > R + TOL or ir.y0 < T - TOL or ir.y1 > B + TOL:
                issues.append("image out of print area")

        if not text_blocks and not img_rects:
            issues.append("blank page")
        if not re.search(r"Page\s+%d\b" % i, page.get_text()):
            issues.append("footer page number missing")

        pix = page.get_pixmap(dpi=110)
        pix.save(str(RENDER / ("page_%02d.png" % i)))
        status = "OK" if not issues else "ISSUE"
        if issues:
            n_bad += 1
        rows.append([i, W, H, len(text_blocks), len(img_rects), status,
                     "; ".join(issues)[:300]])

    with open(OUT_TSV, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["page", "width_pt", "height_pt", "text_blocks", "images",
                    "status", "issues"])
        w.writerows(rows)

    print("[09]    %d   -> %s" % (len(rows), RENDER))
    for r in rows:
        print("   p%-3d %-6s blocks=%-4d imgs=%-2d %s" % (r[0], r[5], r[3], r[4], r[6]))
    print("[09]    -> %s" % OUT_TSV)
    print("[09]       %d / %d" % (n_bad, len(rows)))
    print("[09]                          /  /   "
          "          184       ")
    if n_bad:
        sys.exit(1)


if __name__ == "__main__":
    main()
