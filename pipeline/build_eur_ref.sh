#!/usr/bin/env bash
# ===========================================================================
#
#
# ===========================================================================
set -u
ROOT="${HORZ_ROOT:-D:/Horzontal_project}"
SRC="$ROOT/gwas_data/1000g_eur"
DST="$SRC/eur_bfile"
PLINK="$ROOT/gwas_data/tools/plink.exe"
KEEP="$ROOT/20B3/work/data/1000g_eur_keep2.txt"
BEDBASE="ALL.chr%d.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"

mkdir -p "$DST"
[ -f "$KEEP" ] || { echo "[  ]        $KEEP"; exit 1; }
[ -f "$PLINK" ] || { echo "[  ]    plink: $PLINK"; exit 1; }

CHRS="${1:-1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22}"

for c in $CHRS; do
  vcf="$SRC/$(printf "$BEDBASE" "$c")"
  out="$DST/eur_chr$c"
  log="$DST/_build_$c.log"

  if [ -s "$out.bed" ] && [ "$(stat -c %s "$out.bed" 2>/dev/null || echo 0)" -gt 1000000 ]; then
    echo "[  ] chr$c     ($(stat -c %s "$out.bed") B)"
    continue
  fi
  if [ ! -s "$vcf" ]; then
    echo "[  ] chr$c   VCF     $vcf"
    continue
  fi

  echo "[  ] chr$c ..."
  "$PLINK" --vcf "$vcf" --keep "$KEEP" --maf 0.01 --geno 0.05 \
           --set-missing-var-ids @:#:\$r:\$a --make-bed --out "$out" --memory 6000 > "$log" 2>&1
  if [ -s "$out.bed" ]; then
    nvar=$(wc -l < "$out.bim")
    nsmp=$(wc -l < "$out.fam")
    echo "[  ] chr$c   =$nvar   =$nsmp   =$(stat -c %s "$out.bed") B"
  else
    echo "[  ] chr$c   $log"
  fi
done

echo "[  ]     "
ls -1 "$DST"/eur_chr*.bed 2>/dev/null | wc -l
