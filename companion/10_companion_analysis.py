# -*- coding: utf-8 -*-
"""Robustness (trait permutation, leave-one-trait-out), consequence analysis (BH counts), f_block bootstrap, Figure 1."""
import csv
import json
import math
import statistics
import sys
from pathlib import Path

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = Path("D:/Horzontal_project/20B3/work/crossdomain")
OUT = ROOT / "output"
PAIRS = OUT / "pairs.csv"

LACTASE_F = 0.583
NONLACTASE_F = 0.996
N_PERM = 20000
RNG = np.random.default_rng(20261030)


def read_rows():
    with open(PAIRS, encoding="utf-8", newline="") as fh:
        return [r for r in csv.DictReader(fh, delimiter="\t") if int(r["k"]) >= 2]


def norm_p(b, se):
    """   P   erfc       scipy    """
    if not se or se <= 0:
        return 1.0
    z = abs(b / se)
    return math.erfc(z / math.sqrt(2.0))


def bh(pvals, alpha=0.05):
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    keep = [False] * n
    thresh = 0
    for rank, i in enumerate(order, 1):
        if pvals[i] <= alpha * rank / n:
            thresh = rank
    for rank, i in enumerate(order, 1):
        if rank <= thresh:
            keep[i] = True
    return keep


def rank(v):
    o = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(o):
        j = i
        while j + 1 < len(o) and v[o[j + 1]] == v[o[i]]:
            j += 1
        a = (i + j) / 2.0 + 1
        for k in range(i, j + 1):
            r[o[k]] = a
        i = j + 1
    return r


def spearman(x, y):
    rx, ry = rank(x), rank(y)
    n = len(x)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    return num / (dx * dy) if dx and dy else float("nan")


def perm_test(rows, x_field, y_field, more_negative=False):
    """  trait    x    trait           p   """
    by_trait = {}
    for r in rows:
        by_trait.setdefault(r["gcst"], []).append(r)
    traits = list(by_trait)
    obs = spearman([float(r[x_field]) for r in rows], [float(r[y_field]) for r in rows])
    counts = 0
    for _ in range(N_PERM):
        vals = [statistics.median([float(r[x_field]) for r in by_trait[t]]) for t in traits]
        vals = list(RNG.permutation(vals))
        x = []
        for v, t in zip(vals, traits):
            x.extend([v] * len(by_trait[t]))
        s = spearman(x, [float(r[y_field]) for r in rows])
        if more_negative:
            counts += (s <= obs)
        else:
            counts += (s >= obs)
    return obs, (counts + 1) / (N_PERM + 1)


def main():
    rows = read_rows()
    by_trait = {}
    for r in rows:
        by_trait.setdefault(r["trait"], []).append(r)
    traits = sorted(by_trait)

    p1_obs, p1_p = perm_test(rows, "max_window_share", "f")
    p4_obs, p4_p = perm_test(rows, "k_kept", "f", more_negative=True)
    loo = []
    for t in traits:
        sub = [r for r in rows if r["trait"] != t]
        loo.append((t, spearman([float(r["max_window_share"]) for r in sub],
                                [float(r["f"]) for r in sub])))
    loo.sort(key=lambda x: x[1])

    tf = [statistics.median([float(r["f"]) for r in by_trait[t]]) for t in traits]
    tc = [statistics.median([float(r["max_window_share"]) for r in by_trait[t]]) for t in traits]
    trait_rho = spearman(tc, tf)

    p_un = [norm_p(float(r["b_unclumped"]), float(r["se_unclumped"])) for r in rows]
    p_cl = [norm_p(float(r["b_clumped"]), float(r["se_clumped"])) for r in rows]
    sig_un = bh(p_un)
    sig_cl = bh(p_cl)
    n = len(rows)
    rate_un = sum(sig_un) / n * 100
    rate_cl = sum(sig_cl) / n * 100
    fold = (sum(sig_un) / n) / (sum(sig_cl) / n) if sum(sig_cl) else float("inf")

    fb = [float(r["f_block"]) for r in rows if r["f_block"] not in ("", None)]
    med_fb = statistics.median(fb)
    boots = []
    for _ in range(5000):
        s = list(RNG.choice(fb, size=len(fb), replace=True))
        boots.append(statistics.median(s))
    boots.sort()
    ci = (boots[int(0.025 * len(boots))], boots[int(0.975 * len(boots)) - 1])

    make_figure(rows)

    lines = [
        "# companion       20B3         110   ",
        "",
        "## A.             ",
        "",
        "110     10       pair   n       **trait      **"
        "     trait          %d    " % N_PERM,
        "",
        "|    |    rho |      p |",
        "|---|---|---|",
        "| P1    ↔f | %.3f | **%.4f** |" % (p1_obs, p1_p),
        "| P4 block  ↔f | %.3f | **%s** |"
        % (p4_obs, ("<0.0001" if p4_p <= 1.0 / (N_PERM + 1) else "%.4f" % p4_p)),
        "",
        "trait   n=%d Spearman(     , f   ) = **%.3f** " % (len(traits), trait_rho),
        "",
        "   trait   P1    **%.3f–%.3f**    %s      "
        % (loo[0][1], loo[-1][1], loo[0][0]),
        "",
        "## B.      clump            182/575 → 2/548        ",
        "",
        "  110            P    BH α=0.05  ",
        "",
        "|    |    BH |    |",
        "|---|---|---|",
        "|   clump | **%d / %d** | **%.2f%%** |" % (sum(sig_un), n, rate_un),
        "| clumped | **%d / %d** | **%.2f%%** |" % (sum(sig_cl), n, rate_cl),
        "",
        ("         **    ** clumped   0            "
         "   31.6%% → 0.36%% 86.7        **%.1f%% → 0%%**——    "
         "        clumped    p<0.05    %d      p=%.3g "
         "       BH ")
        % (rate_un, sum(1 for x in p_cl if x < 0.05), min(p_cl)),
        "",
        "## C.   block             0.583     ",
        "",
        "f_block    **%.3f** bootstrap 95%% CI %.3f–%.3f n=%d  "
        % (med_fb, ci[0], ci[1], len(fb)),
        "       f = %.3f %s 95%% CI   "
        % (LACTASE_F, "  " if ci[0] <= LACTASE_F <= ci[1] else "**  **"),
        "",
        "## D.  ",
        "",
        "`output/fig1_crossdomain.png` f     10 Mb       =k_kept "
        "          0.583        0.996 ",
    ]
    txt = "\n".join(lines)
    (OUT / "companion_stats.md").write_text(txt, encoding="utf-8", newline="")
    (OUT / "companion.json").write_text(json.dumps({
        "n_pairs": n, "n_traits": len(traits),
        "p1_obs": p1_obs, "p1_perm_p": p1_p, "p4_obs": p4_obs, "p4_perm_p": p4_p,
        "trait_level_rho": trait_rho,
        "loo_min": {"trait": loo[0][0], "rho": loo[0][1]},
        "loo_max": {"trait": loo[-1][0], "rho": loo[-1][1]},
        "sig_unclumped": int(sum(sig_un)), "sig_clumped": int(sum(sig_cl)),
        "rate_unclumped": rate_un, "rate_clumped": rate_cl, "rate_fold": fold,
        "f_block_median": med_fb, "f_block_ci": ci,
        "lactase_reference": LACTASE_F,
    }, ensure_ascii=False, indent=1), encoding="utf-8", newline="")
    print(txt)


def make_figure(rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.2, 5.0), dpi=300)
    xs = [float(r["max_window_share"]) for r in rows]
    ys = [float(r["f"]) for r in rows]
    ks = [int(r["k_kept"]) for r in rows]
    sc = ax.scatter(xs, ys, c=ks, cmap="viridis", s=38, edgecolor="white", linewidth=0.5,
                    vmin=min(ks), vmax=max(ks))
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label("independent blocks after clumping (k_kept)", fontsize=10)
    ax.axhline(LACTASE_F, ls="--", lw=1.2, color="#c0392b")
    ax.text(0.02, LACTASE_F + 0.02, "lactase axis (microbiome), f = 0.583",
            color="#c0392b", fontsize=8)
    ax.axhline(NONLACTASE_F, ls=":", lw=1.2, color="#7f8c8d")
    ax.text(0.02, NONLACTASE_F - 0.045, "non-lactase axis (microbiome), f = 0.996",
            color="#7f8c8d", fontsize=8)
    ax.set_xlabel("share of variants in the largest 10 Mb window", fontsize=11)
    ax.set_ylabel("attained fraction of the perfect-redundancy deflation (f)", fontsize=11)
    ax.set_ylim(0.1, 1.05)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(alpha=0.25, lw=0.5)
    fig.tight_layout()
    fig.savefig(OUT / "fig1_crossdomain.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
