# Manufactured precision in Mendelian randomization: microbiome and metabolite pipeline

Code for two companion manuscripts:

1. **Microbiome paper** — "LD clumping removes almost all apparent gut microbiome-infection
   Mendelian randomization signals: manufactured precision follows a perfect-redundancy bound"
   (MiBioGen 211 taxa x 23 FinnGen R12 infection endpoints).
2. **Metabolite paper** — "The diagnostic transfers: manufactured precision in metabolite
   Mendelian randomization" (10 NMR metabolite traits x 11 infection endpoints, with an
   external replication in two UK Biobank outcomes).

Both papers introduce the same two diagnostics, computable from any analysis that reports both
instrument definitions:

- `k_eff = (SE_clumped / SE_unclumped)^2` — the effective number of independent instruments;
- `f = (SE_clumped / SE_unclumped) / sqrt(k)` — the attained fraction of the maximal
  perfect-redundancy deflation.

## Repository layout

```
pipeline/    main microbiome pipeline + journal-quality tools
companion/   metabolite cross-domain scripts
```

This repository is **code and README only**. Small static inputs (exposure-library study
tables, the 1000 Genomes EUR keep-list and panel, FinnGen phenotype counts) and the frozen
result tables are not included; the scripts document where each input is expected, and
`frozen_results`-style tables can be regenerated with the reproduction order below.

The GWAS Catalog studies index (126 MB, used only for planning queries) is not included;
download the current version from the GWAS Catalog downloads page if you need it.

## Data sources (all public)

- MiBioGen top-hit table (Molgenis REST API)
- GWAS Catalog summary statistics (Kettunen 2016 NMR metabolomics, GCST90132xxx)
- FinnGen release 12 (remote tabix over the public bucket)
- 1000 Genomes phase 3 EUR panel (genotypes)
- OpenGWAS (UK Biobank outcomes; optional, requires a free JWT token)

## Requirements

- Python 3.10+ — standard library plus `openpyxl`, `numpy`, `scipy`, `matplotlib`, `python-docx`
- R 4.x — `data.table`, `ggplot2`, `coloc`, `susieR`, `MRPRESSO`

## Reproduction order (microbiome paper)

1. `02_build_instruments.py` -> `02b_clump_python.py` -> `02c_clump_global.py`
2. `03a_extract_outcomes.py` -> `03b_mr_main.R` / `03c_mr_clumped.R`
3. `04_replication_finnish.py`, `04_sensitivity.R`, `04b_robust_preclump.R`,
   `04e_coloc.R`, `04f_steiger_real.R`
4. `04h_inflation_law.R`, `04i_law_external_validation.py`, `04l_loose_clump_sensitivity.py`
5. `05_canonical_numbers.py` (freezes every reported value), `06_verify_refs.py`,
   `07_supp_tables.py`
6. `pipeline/08_source_check.py` (126 independent recomputation checks)

Metabolite companion: `companion/09_crossdomain_pilot.py` (steps `fetch/map/outcome/mr/stats`),
`10_companion_analysis.py`, `11_companion_tables.py`, `12_companion_supp.py`,
`14_opengwas_external.py`, `15_companion_check.py`.

## Notes

- Scripts contain **absolute paths** from the original run (a deliberate provenance choice);
  replace the `ROOT`/path constants with your own locations to rerun. Everything the pipeline
  downloads is public; the GWAS Catalog FTP is slow (~0.06 MB/s observed), and the
  cross-domain fetcher supports resumable chunked downloads.
- Working comments and runtime messages were written in Chinese and have been **removed**
  (not machine-translated) in this public copy; each file carries a concise English docstring.
  The original annotated scripts are available from the authors on reasonable request.
- The OpenGWAS token is read from the environment variable `OPENGAWS_TOKEN` (fallback:
  `~/.zcode/secrets/opengwas_token.txt`). It is **never** stored in this repository.
- The LD-aware GLS script (`04k`) is archived for transparency but its results are **not**
  reported in the manuscripts (shrinkage-unstable).
