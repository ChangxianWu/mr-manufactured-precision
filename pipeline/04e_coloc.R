# ============================================================
# ------------------------------------------------------------------
#
#
#
# ============================================================
.libPaths("D:/Rlibs/4.6")
source("D:/Horzontal_project/20B3/work/course_config.R")
suppressPackageStartupMessages({ library(data.table); library(coloc) })

WINDOW  <- 100000
P12_SET <- c(1e-6, 1e-5, 1e-4)
E19_F  <- file.path(DATA, "coloc/GCST90017066_chr19_region.tsv")
E2_F   <- file.path(DATA, "coloc/GCST90016970_chr2_region.tsv")
LD19_F <- file.path(CACHE, "ld/ld_chr19.tsv")
LD2_F  <- file.path(CACHE, "ld/ld_chr2.tsv")

LEAD <- list(
  chr19 = list(chr = 19L, pos_exp = 49218060L, pos_out = 48714803L, rs = "rs35866622"),
  chr2  = list(chr = 2L,  pos_exp = 136616754L, pos_out = 135859184L, rs = "rs182549"))

ENDPOINTS_FOCAL <- c("INTEST_INF", "GASTROENT", "SEPSIS", "BACT_OTHER")

build_exp <- function(path, tag) {
  E <- data.table::fread(path, data.table = FALSE, showProgress = FALSE)
  E <- E[!is.na(E$beta) & !is.na(E$standard_error) & E$standard_error > 0 &
           E$variant_id != "" & E$variant_id != ".", ]
  L <- LEAD[[tag]]
  E <- E[E$chromosome == L$chr & abs(E$base_pair_location - L$pos_exp) <= WINDOW, ]
  E <- E[!duplicated(E$variant_id), ]
  cat(sprintf("[04e] %s      ±%d kb   SNP %d %s \n", tag, WINDOW / 1000,
              nrow(E), unique(E$bac)[1]))
  E
}

expand_outcome <- function(O, tag) {
  L <- LEAD[[tag]]
  rows <- list()
  for (i in seq_len(nrow(O))) {
    if (abs(as.integer(O$pos[i]) - L$pos_out) > WINDOW) next
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
  cat(sprintf("[04e] LD    %d × %d %s \n", nrow(R), ncol(R), basename(path)))
  R
}

run_locus <- function(E, key, out_file, tag, Rld) {
  if (!file.exists(out_file)) {
    cat(sprintf("[04e] %s %-12s      \n", tag, key)); return(NULL) }
  O <- data.table::fread(out_file, data.table = FALSE, showProgress = FALSE)
  Oe <- expand_outcome(O, tag)
  if (is.null(Oe) || !nrow(Oe)) {
    cat(sprintf("[04e] %s %-12s        SNP\n", tag, key)); return(NULL) }

  M <- merge(
    data.frame(rsid = E$variant_id, beta_exp = as.numeric(E$beta),
               se_exp = as.numeric(E$standard_error), n_exp = as.numeric(E$n),
               stringsAsFactors = FALSE),
    data.frame(rsid = Oe$rsid, beta_out = Oe$beta, se_out = Oe$se,
               eaf_out = Oe$eaf, p_out = Oe$pval, stringsAsFactors = FALSE),
    by = "rsid")
  M$maf <- pmin(M$eaf_out, 1 - M$eaf_out)
  M <- M[M$maf > 0.01, ]
  if (nrow(M) < 50) {
    cat(sprintf("[04e] %s %-12s    SNP   %d   \n", tag, key, nrow(M))); return(NULL) }
  M$varbeta_exp <- M$se_exp^2
  M$varbeta_out <- M$se_out^2

  ntot <- ENDPOINT_CASES[[key]] + ENDPOINT_CONTROLS[[key]]
  d2 <- list(beta = M$beta_out, varbeta = M$varbeta_out, snp = M$rsid,
             MAF = M$maf, N = ntot,
             s = ENDPOINT_CASES[[key]] / ntot, type = "cc")

  out <- list()
  for (p12 in P12_SET) {
    d1 <- list(beta = M$beta_exp, varbeta = M$varbeta_exp, snp = M$rsid,
               MAF = M$maf, N = max(M$n_exp, na.rm = TRUE), type = "quant")
    r <- tryCatch(coloc.abf(d1, d2, p12 = p12), error = function(e) NULL)
    if (is.null(r)) next
    s <- r$summary
    out[[length(out) + 1]] <- data.frame(
      locus = tag, endpoint = key, endpoint_label = ENDPOINT_LABEL[[key]],
      method = "coloc.abf", p12 = p12, n_snp_shared = nrow(M),
      PP_H0 = s[["PP.H0.abf"]], PP_H1 = s[["PP.H1.abf"]], PP_H2 = s[["PP.H2.abf"]],
      PP_H3 = s[["PP.H3.abf"]], PP_H4 = s[["PP.H4.abf"]],
      ratio_H4_H3 = s[["PP.H4.abf"]] / max(s[["PP.H3.abf"]], 1e-12),
      verdict = ifelse(s[["PP.H4.abf"]] >= 0.8, "colocalised",
                       ifelse(s[["PP.H4.abf"]] >= 0.5, "suggestive", "not_supported")),
      stringsAsFactors = FALSE)
  }

  keep <- intersect(M$rsid, rownames(Rld))
  if (length(keep) >= 50) {
    R <- Rld[keep, keep, drop = FALSE]
    Ms <- M[match(keep, M$rsid), ]
    s1 <- list(beta = Ms$beta_exp, varbeta = Ms$varbeta_exp, snp = Ms$rsid,
               MAF = Ms$maf, N = max(Ms$n_exp, na.rm = TRUE), type = "quant", LD = R)
    s2 <- list(beta = Ms$beta_out, varbeta = Ms$varbeta_out, snp = Ms$rsid,
               MAF = Ms$maf, N = ntot,
               s = ENDPOINT_CASES[[key]] / ntot, type = "cc", LD = R)
    rs <- tryCatch(suppressWarnings(coloc.susie(s1, s2, p12 = 1e-5)),
                   error = function(e) NULL)
    if (!is.null(rs) && !is.null(rs$summary) && nrow(rs$summary)) {
      ss <- rs$summary
      j <- which.max(ss$PP.H4)
      out[[length(out) + 1]] <- data.frame(
        locus = tag, endpoint = key, endpoint_label = ENDPOINT_LABEL[[key]],
        method = "coloc.susie", p12 = 1e-5, n_snp_shared = length(keep),
        PP_H0 = NA_real_, PP_H1 = NA_real_, PP_H2 = NA_real_,
        PP_H3 = ss$PP.H3[j], PP_H4 = ss$PP.H4[j],
        ratio_H4_H3 = ss$PP.H4[j] / max(ss$PP.H3[j], 1e-12),
        verdict = ifelse(ss$PP.H4[j] >= 0.8, "colocalised",
                         ifelse(ss$PP.H4[j] >= 0.5, "suggestive", "not_supported")),
        stringsAsFactors = FALSE)
    } else {
      out[[length(out) + 1]] <- data.frame(
        locus = tag, endpoint = key, endpoint_label = ENDPOINT_LABEL[[key]],
        method = "coloc.susie", p12 = 1e-5, n_snp_shared = length(keep),
        PP_H0 = NA_real_, PP_H1 = NA_real_, PP_H2 = NA_real_,
        PP_H3 = NA_real_, PP_H4 = NA_real_, ratio_H4_H3 = NA_real_,
        verdict = "susie_not_converged", stringsAsFactors = FALSE)
    }
  }
  D <- do.call(rbind, out)
  if (!is.null(D)) {
    a <- D$PP_H4[D$method == "coloc.abf" & D$p12 == 1e-5][1]
    v <- D$verdict[D$method == "coloc.susie"][1]
    cat(sprintf("[04e] %s %-12s    SNP %5d | abf PP.H4=%.3f | susie: %s\n",
                tag, key, nrow(M), a, v))
  }
  D
}

LD19 <- if (file.exists(LD19_F)) load_ld(LD19_F) else NULL
LD2  <- if (file.exists(LD2_F))  load_ld(LD2_F)  else NULL

res <- list()
if (file.exists(E19_F) && !is.null(LD19)) {
  E19 <- build_exp(E19_F, "chr19")
  for (key in ENDPOINTS_FOCAL) {
    r <- run_locus(E19, key,
                   file.path(DATA, sprintf("coloc/outcome_%s_chr19_region.tsv", key)),
                   "chr19", LD19)
    if (!is.null(r)) res[[paste0("chr19_", key)]] <- r
  }
}
if (file.exists(E2_F) && !is.null(LD2)) {
  E2 <- build_exp(E2_F, "chr2")
  for (key in ENDPOINTS_FOCAL) {
    r <- run_locus(E2, key,
                   file.path(DATA, sprintf("coloc/outcome_%s_chr2_region.tsv", key)),
                   "chr2", LD2)
    if (!is.null(r)) res[[paste0("chr2_", key)]] <- r
  }
}
if (!length(res)) stop("[04e]       ")
COL <- do.call(rbind, res)
COL <- COL[order(COL$locus, COL$endpoint, COL$method, COL$p12), ]
data.table::fwrite(COL, file.path(OUT, "coloc_focal.csv"))

BEST <- do.call(rbind, lapply(split(COL, paste(COL$locus, COL$endpoint)), function(g) {
  a <- g[g$method == "coloc.abf" & g$p12 == 1e-5, ]
  s <- g[g$method == "coloc.susie", ]
  data.frame(locus = g$locus[1], endpoint = g$endpoint[1],
             abf_PP_H3 = a$PP_H3[1], abf_PP_H4 = a$PP_H4[1],
             abf_PP_H1 = a$PP_H1[1],
             susie_PP_H4 = if (nrow(s)) s$PP_H4[1] else NA_real_,
             susie_verdict = if (nrow(s)) s$verdict[1] else "not_run",
             stringsAsFactors = FALSE)
}))
data.table::fwrite(BEST, file.path(OUT, "coloc_focal_best.csv"))

cat(sprintf("\n[04e]      %d    ±%d kb    %s p12 %s\n",
            nrow(COL), WINDOW / 1000, paste(unique(COL$method), collapse = "/"),
            paste(P12_SET, collapse = "/")))
cat(sprintf("[04e] chr19 abf    PP.H4=%.3f chr2 abf    PP.H4=%.3f\n",
            max(COL$PP_H4[COL$locus == "chr19" & COL$method == "coloc.abf"], na.rm = TRUE),
            max(COL$PP_H4[COL$locus == "chr2" & COL$method == "coloc.abf"], na.rm = TRUE)))
print(BEST, row.names = FALSE)
cat("[04e] -> coloc_focal.csv / coloc_focal_best.csv\n")
