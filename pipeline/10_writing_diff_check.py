# -*- coding: utf-8 -*-
"""Numeric-difference check between manuscript versions."""
import re
import sys
from pathlib import Path

import docx

ALLOW_DELETE = set()


def numbers(path):
    d = docx.Document(str(path))
    txt = " ".join(p.text for p in d.paragraphs)
    for t in d.tables:
        for r in t.rows:
            txt += " " + " ".join(c.text for c in r.cells)
    out = set()
    for m in re.finditer(r"\d[\d,]*\.?\d*(?:[eE][-+]?\d+)?", txt):
        out.add(m.group(0).replace(",", ""))
    return out, len(txt.split())


def main():
    old, new = Path(sys.argv[1]), Path(sys.argv[2])
    o, ow = numbers(old)
    n, nw = numbers(new)
    added = sorted(n - o)
    removed = sorted(o - n)
    bad = []
    if added:
        bad.append("              : " + ", ".join(added))
    unreg = sorted(set(removed) - ALLOW_DELETE)
    if unreg:
        bad.append("        : " + ", ".join(unreg))
    out = Path(r"D:/Horzontal_project/20B3/work/output/writing_diff.txt")
    lines = [
        "old=%s words=%d numbers=%d" % (old.name, ow, len(o)),
        "new=%s words=%d numbers=%d" % (new.name, nw, len(n)),
        "added=%s" % (added or "[]"),
        "removed=%s" % (removed or "[]"),
        "allowed_delete=%s" % (sorted(ALLOW_DELETE) or "[]"),
    ]
    out.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="")
    for line in lines:
        print("  " + line)
    if bad:
        print("[10] WRITING_DIFF_FAIL")
        for b in bad:
            print("   - " + b)
        sys.exit(1)
    print("[10] WRITING_DIFF_OK")


if __name__ == "__main__":
    main()
