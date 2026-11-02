# ============================================================
# ------------------------------------------------------------------
#   B) rs736744   × BACT_OTHER     —— Oxalobacter chr9
#
#
#
# ============================================================
.libPaths("D:/Rlibs/4.6")
source("D:/Horzontal_project/20B3/work/course_config.R")
suppressPackageStartupMessages({ library(data.table); library(coloc) })

WINDOW  <- 100000
P12_SET <- c(1e-6, 1e-5, 1e-4)
MIN_SHARED <- 50

LEAD37 <- list(
  chr19 = list(chr = 19L, pos = 49218060L, rs = "rs35866622"),
  chr9  = list(chr = 9L,  pos = 87514407L, rs = "rs736744"))
LEAD38 <- list(
  chr19 = 48714803L,
  chr9  = 84899492L)

PAIRS <- list(
  list(tag = "chr19", key = "PNEUMOBACT", taxon = "genus..Ruminococcustorquesgroup.id.14377",
       exp_file = file.path(DATA, "coloc/GCST90017066_chr19_region.tsv"),
       out_file = file.path(DATA, "coloc/outcome_PNEUMOBACT_chr19_region.tsv"),
       ld_file  = file.path(CACHE, "ld/ld_chr19.tsv")),
  list(tag = "chr9", key = "BACT_OTHER", taxon = "genus.Oxalobacter.id.2978",
       exp_file = file.path(DATA, "coloc/GCST90017038_chr9_region.tsv"),
       out_file = file.path(DATA, "coloc/outcome_BACT_OTHER_chr9_region.tsv"),
       ld_file  = NULL))

build_exp <- function(path, tag) {
  E <- data.table::fread(path, data.table = FALSE, showProgress = FALSE)
  E <- E[!is.na(E$beta) & !is.na(E$standard_error) & E$standard_error > 0 &
           E$variant_id != "" & E$variant_id != ".", ]
  L <- LEAD37[[tag]]
  E <- E[E$chromosome == L$chr & abs(E$base_pair_location - L$pos) <= WINDOW, ]
  E <- E[!duplicated(E$variant_id), ]
  cat(sprintf("[04o] %s      ±%d kb   SNP %d %s \n", tag, WINDOW / 1000,
              nrow(E), unique(E$bac)[1]))
  E
}

expand_outcome <- function(O, tag) {
  pos_out <- LEAD38[[tag]]
  rows <- list()
  for (i in seq_len(nrow(O))) {
    if (abs(as.integer(O$pos[i]) - pos_out) > WINDOW) next
    rs <- trimws(strsplit(O$rsids[i], ",")[[1]])
    rs <- rs[grepl("^rs[0-9]+$", rs)]
    for (r in rs) {
      rows[[length(rows) + 1]] <- data.frame(
        rsid = r, beta = as.numeric(O$beta[i]), se = as.numeric(O$sebeta[i]),
        eaf = as.numeric(O$af_alt[i]), pval = as.numeric(O$pval[i]),
        stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(NULL)
  D <- do.call(rbind, rows)
  D <- D[!is.na(D$beta) & !is.na(D$se) & D$se > 0 &
           !is.na(D$eaf) & D$eaf > 0 & D$eaf < 1, ]
  D[!duplicated(D$rsid), ]
}

load_ld <- function(path) {
  L <- data.table::fread(path, data.table = FALSE, showProgress = FALSE)
  ids <- L$rsid
  M <- as.matrix(L[, -(1:4)])
  M[is.na(M)] <- 0
  G <- scale(t(M))
  R <- suppressWarnings(cor(G))
  R[is.na(R)] <- 0
  diag(R) <- 1
  rownames(R) <- colnames(R) <- ids
  cat(sprintf("[04o] LD    %d × %d %s \n", nrow(R), ncol(R), basename(path)))
  R
}

run_pair <- function(P) {
  st <- list(locus = P$tag, endpoint = P$key, taxon = P$taxon,
             n_exp_region = NA_integer_, n_out_region = NA_integer_,
             n_shared = NA_integer_, status = "", note = "")
  if (!file.exists(P$exp_file)) {
    st$status <- "exposure_region_unavailable"
    st$note <- paste("             p          tabix    ",
                     "            '  '  '    '")
    if (file.exists(P$out_file)) {
      st$n_out_region <- nrow(data.table::fread(P$out_file, data.table = FALSE,
                                                showProgress = FALSE))
    }
    return(st)
  }
  if (!file.exists(P$out_file)) {
    st$status <- "input_missing"
    st$note <- "        "
    return(st)
  }
  E <- build_exp(P$exp_file, P$tag)
  O <- data.table::fread(P$out_file, data.table = FALSE, showProgress = FALSE)
  Oe <- expand_outcome(O, P$tag)
  st$n_exp_region <- nrow(E)
  if (is.null(Oe) || !nrow(Oe)) {
    st$status <- "outcome_region_empty"; return(st)
  }
  st$n_out_region <- nrow(Oe)

  M <- merge(
    data.frame(rsid = E$variant_id, beta_exp = as.numeric(E$beta),
               se_exp = as.numeric(E$standard_error), n_exp = as.numeric(E$n),
               stringsAsFactors = FALSE),
    data.frame(rsid = Oe$rsid, beta_out = Oe$beta, se_out = Oe$se,
               eaf_out = Oe$eaf, stringsAsFactors = FALSE),
    by = "rsid")
  M$maf <- pmin(M$eaf_out, 1 - M$eaf_out)
  M <- M[M$maf > 0.01, ]
  st$n_shared <- nrow(M)
  if (nrow(M) < MIN_SHARED) {
    st$status <- "insufficient_shared_variants"
    st$note <- sprintf("   SNP %d < %d    abf                 ",
                       nrow(M), MIN_SHARED)
    return(st)
  }

  M$varbeta_exp <- M$se_exp^2
  M$varbeta_out <- M$se_out^2
  ntot <- ENDPOINT_CASES[[P$key]] + ENDPOINT_CONTROLS[[P$key]]
  d2 <- list(beta = M$beta_out, varbeta = M$varbeta_out, snp = M$rsid,
             MAF = M$maf, N = ntot,
             s = ENDPOINT_CASES[[P$key]] / ntot, type = "cc")

  best <- list()
  for (p12 in P12_SET) {
    d1 <- list(beta = M$beta_exp, varbeta = M$varbeta_exp, snp = M$rsid,
               MAF = M$maf, N = max(M$n_exp, na.rm = TRUE), type = "quant")
    r <- tryCatch(coloc.abf(d1, d2, p12 = p12), error = function(e) NULL)
    if (is.null(r)) next
    s <- r$summary
    best[[length(best) + 1]] <- list(
      method = "coloc.abf", p12 = p12, n = nrow(M),
      h1 = s[["PP.H1.abf"]], h3 = s[["PP.H3.abf"]], h4 = s[["PP.H4.abf"]],
      verdict = ifelse(s[["PP.H4.abf"]] >= 0.8, "colocalised",
                       ifelse(s[["PP.H4.abf"]] >= 0.5, "suggestive", "not_supported")))
  }
  if (is.null(P$ld_file) || !file.exists(P$ld_file)) {
    st$status <- "abf_only_no_ld"
    st$note <- "  1000G EUR LD       coloc.susie"
  } else {
    Rld <- load_ld(P$ld_file)
    keep <- intersect(M$rsid, rownames(Rld))
    if (length(keep) >= MIN_SHARED) {
      Ms <- M[match(keep, M$rsid), ]
      R <- Rld[keep, keep, drop = FALSE]
      s1 <- list(beta = Ms$beta_exp, varbeta = Ms$varbeta_exp, snp = Ms$rsid,
                 MAF = Ms$maf, N = max(Ms$n_exp, na.rm = TRUE), type = "quant", LD = R)
      s2 <- list(beta = Ms$beta_out, varbeta = Ms$varbeta_out, snp = Ms$rsid,
                 MAF = Ms$maf, N = ntot,
                 s = ENDPOINT_CASES[[P$key]] / ntot, type = "cc", LD = R)
      rs <- tryCatch(suppressWarnings(coloc.susie(s1, s2, p12 = 1e-5)),
                     error = function(e) NULL)
      if (!is.null(rs) && !is.null(rs$summary) && nrow(rs$summary)) {
        j <- which.max(rs$summary$PP.H4)
        best[[length(best) + 1]] <- list(
          method = "coloc.susie", p12 = 1e-5, n = length(keep),
          h1 = NA_real_, h3 = rs$summary$PP.H3[j], h4 = rs$summary$PP.H4[j],
          verdict = ifelse(rs$summary$PP.H4[j] >= 0.8, "colocalised",
                           ifelse(rs$summary$PP.H4[j] >= 0.5, "suggestive",
                                  "not_supported")))
        st$status <- "ran_abf_susie"
      } else {
        st$status <- "ran_abf_susie_not_converged"
      }
    } else {
      st$status <- "ran_abf_ld_coverage_low"
      st$note <- sprintf("LD       SNP   %d    susie", length(keep))
    }
  }
  if (!length(best)) { st$status <- "abf_failed"; return(st) }
  idx <- which(vapply(best, function(x) x$method == "coloc.abf" && x$p12 == 1e-5, TRUE))
  a <- best[[if (length(idx)) idx[1] else 1]]
  cat(sprintf("[04o] %-5s %-12s    SNP %4d | abf p12=1e-5 PP.H4=%.3f (%s)\n",
              P$tag, P$key, a$n, a$h4, a$verdict))
  attr(st, "rows") <- best
  st
}

res <- lapply(PAIRS, run_pair)
ROWS <- do.call(rbind, lapply(res, function(s) {
  b <- attr(s, "rows")
  if (is.null(b)) return(NULL)
  do.call(rbind, lapply(b, function(x) data.frame(
    locus = s$locus, endpoint = s$endpoint, taxon = s$taxon, method = x$method,
    p12 = x$p12, n_snp_shared = x$n, PP_H1 = x$h1, PP_H3 = x$h3, PP_H4 = x$h4,
    verdict = x$verdict, stringsAsFactors = FALSE)))
}))
SUM <- do.call(rbind, lapply(res, function(s) data.frame(
  locus = s$locus, endpoint = s$endpoint, taxon = s$taxon,
  n_exp_region = s$n_exp_region, n_out_region = s$n_out_region,
  n_snp_shared = s$n_shared, status = s$status, note = s$note,
  stringsAsFactors = FALSE)))

data.table::fwrite(SUM, file.path(OUT, "loose_coloc_summary.csv"))
if (!is.null(ROWS)) data.table::fwrite(ROWS, file.path(OUT, "loose_coloc.csv"))
print(SUM, row.names = FALSE)
cat("[04o] -> loose_coloc_summary.csv", if (!is.null(ROWS)) " + loose_coloc.csv" else "",
    "\n")
