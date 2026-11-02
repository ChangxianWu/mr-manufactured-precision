# ============================================================
# ------------------------------------------------------------------
#
# ============================================================
import gzip
import sys

ROOT = "D:/Horzontal_project"
VCF = {
    2: ROOT + "/gwas_data/1000g_eur/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    19: ROOT + "/gwas_data/1000g_eur/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
}
KEEP = ROOT + "/20B3/work/data/1000g_eur_keep2.txt"

POSID = {
    2: ROOT + "/20B3/work/data/coloc/GCST90016970_chr2_region.tsv",
    19: ROOT + "/20B3/work/data/coloc/GCST90017066_chr19_region.tsv",
}


def pos_to_rsid(path):
    """Extract 1000 Genomes EUR dosages for a list of positions (LD matrices for coloc)."""
    import csv
    m = {}
    with open(path, encoding="utf-8") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            v = (r.get("variant_id") or "").strip()
            if v.startswith("rs"):
                m[r["base_pair_location"].strip()] = v
    return m


def main():
    chrom = int(sys.argv[1])
    pos_file = sys.argv[2]
    out_tsv = sys.argv[3]

    with open(pos_file, encoding="utf-8") as fh:
        want = {ln.strip() for ln in fh if ln.strip()}
    with open(KEEP, encoding="utf-8") as fh:
        keep = [ln.split()[0].strip() for ln in fh if ln.strip()]
    p2r = pos_to_rsid(POSID[chrom])
    print("[04g] chr%d      %d / EUR    %d /   ->rsID    %d"
          % (chrom, len(want), len(keep), len(p2r)))

    n_hit = 0
    n_named = 0
    with gzip.open(VCF[chrom], "rt", errors="replace") as fh, \
            open(out_tsv, "w", encoding="utf-8", newline="") as out:
        out.write("pos\trsid\ta1\ta2\t" + "\t".join(keep) + "\n")
        idx = None
        for line in fh:
            if line[0] == "#":
                if line.startswith("#CHROM"):
                    names = line.rstrip("\n").split("\t")[9:]
                    n2i = {n: i for i, n in enumerate(names)}
                    idx = [n2i[s] for s in keep if s in n2i]
                continue
            head = line.split("\t", 4)
            pos = head[1]
            if pos not in want:
                continue
            f = line.rstrip("\n").split("\t")
            rsid = p2r.get(pos, ".")
            if rsid != ".":
                n_named += 1
            a1, a2 = f[3], f[4]
            fmt = f[8].split(":")
            gi = fmt.index("GT") if "GT" in fmt else 0
            doses = []
            for i in idx:
                sample = f[9 + i]
                gt = sample.split(":", gi + 1)[gi] if ":" in sample else sample
                gt = gt.replace("|", "/")
                a = gt.split("/")
                if len(a) != 2 or "." in a:
                    doses.append(".")
                else:
                    doses.append(str(int(a[0]) + int(a[1])))
            out.write("%s\t%s\t%s\t%s\t%s\n" % (pos, rsid, a1, a2,
                                                "\t".join(doses)))
            n_hit += 1
    print("[04g]    %d / %d    %d   rsID -> %s"
          % (n_hit, len(want), n_named, out_tsv))


if __name__ == "__main__":
    sys.exit(main())
