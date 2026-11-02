# ============================================================
# ------------------------------------------------------------------
#
#       work/output/mr_all_endpoints_clumped.csv
# ============================================================
.libPaths("D:/Rlibs/4.6")
source("D:/Horzontal_project/20B3/work/course_config.R")
suppressPackageStartupMessages({ library(data.table) })

TAG <- "clumped"
INSTR_F <- file.path(OUT, "instruments_clumped.csv")
if (!file.exists(INSTR_F)) stop("[03c]    clumped      ", INSTR_F)
INS <- data.table::fread(INSTR_F, data.table = FALSE, showProgress = FALSE)
EXTRACT_DIR <- file.path(DATA, "outcome_extract_clumped")

COMPL <- c(A = "T", T = "A", C = "G", G = "C")
is_pal <- function(a, b) !is.na(a) && !is.na(b) && a == COMPL[[b]]

align_one <- function(ea, oa, ref, alt, beta_out, eaf_out) {
  ea <- toupper(ea); oa <- toupper(oa); ref <- toupper(ref); alt <- toupper(alt)
  if (is.na(ea) || is.na(oa) || is.na(ref) || is.na(alt))
    return(list(beta = NA_real_, status = "missing"))
  if (is_pal(ea, oa) && !is.na(eaf_out) && abs(eaf_out - 0.5) < 0.05)
    return(list(beta = NA_real_, status = "palin_ambig"))
  if (ea == alt && oa == ref) return(list(beta = beta_out, status = "ok"))
  if (ea == ref && oa == alt) return(list(beta = -beta_out, status = "flip"))
  if (is_pal(ea, oa)) return(list(beta = beta_out, status = "palin"))
  return(list(beta = NA_real_, status = "mismatch"))
}

run_endpoint <- function(key) {
  f <- file.path(EXTRACT_DIR, sprintf("outcome_extract_%s.tsv", key))
  if (!file.exists(f)) { cat(sprintf("[03c] %-12s         \n", key)); return(NULL) }
  O <- data.table::fread(f, data.table = FALSE, showProgress = FALSE)
  O$beta_out <- suppressWarnings(as.numeric(O$beta_out))
  O$se_out   <- suppressWarnings(as.numeric(O$se_out))
  O$eaf_out  <- suppressWarnings(as.numeric(O$eaf_out))
  al <- lapply(seq_len(nrow(O)), function(i)
    align_one(O$ea_instr[i], O$oa_instr[i], O$ref_out[i], O$alt_out[i],
              O$beta_out[i], O$eaf_out[i]))
  O$beta_al <- vapply(al, function(x) x$beta, numeric(1))
  O$align   <- vapply(al, function(x) x$status, character(1))
  M <- O[O$align %in% c("ok", "flip") & !is.na(O$beta_al) &
           !is.na(O$se_out) & O$se_out > 0, ]
  if (!nrow(M)) return(NULL)

  res <- lapply(split(M, M$taxon), function(g) {
    bx <- as.numeric(g$beta_instr); by <- as.numeric(g$beta_al)
    sx <- as.numeric(g$se_instr);  sy <- as.numeric(g$se_out)
    k <- length(bx)
    if (k == 1) {
      b <- by / bx
      se <- sqrt(sy^2 / bx^2 + (by^2 * sx^2) / bx^4)
      method <- "Wald ratio"
      q <- NA_real_; q_df <- NA_real_; q_p <- NA_real_; i2 <- NA_real_
      b_re <- b; se_re <- se
    } else {
      w <- 1 / (sy^2 / bx^2)
      b <- sum(w * (by / bx)) / sum(w)
      se <- sqrt(1 / sum(w))
      q <- sum(w * ((by / bx) - b)^2)
      q_df <- k - 1
      q_p <- pchisq(q, df = q_df, lower.tail = FALSE)
      i2 <- max(0, (q - q_df) / q) * 100
      c0 <- sum(w) - sum(w^2) / sum(w)
      tau2 <- max(0, (q - q_df) / c0)
      wr <- 1 / (sy^2 / bx^2 + tau2)
      b_re <- sum(wr * (by / bx)) / sum(wr)
      se_re <- sqrt(1 / sum(wr))
      method <- "IVW (fixed)"
    }
    data.frame(taxon = g$taxon[1], level = g$level[1], name = g$name[1],
               n_snp = k, method = method, b = b, se = se,
               lo = b - 1.96 * se, hi = b + 1.96 * se,
               p = 2 * pnorm(-abs(b / se)), Q = q, Q_df = q_df, Q_p = q_p, I2 = i2,
               b_re = b_re, se_re = se_re,
               p_re = 2 * pnorm(-abs(b_re / se_re)),
               n_flip = sum(g$align == "flip"),
               in_lct = g$in_lct[1], aligned = g$aligned[1],
               stringsAsFactors = FALSE)
  })
  D <- do.call(rbind, res)
  if (is.null(D) || !nrow(D)) return(NULL)
  Fm <- setNames(INS$F_min, INS$taxon)
  Fm2 <- setNames(INS$F_mean, INS$taxon)
  D$F_min <- as.numeric(Fm[D$taxon]); D$F_mean <- as.numeric(Fm2[D$taxon])
  D$endpoint <- key
  D$endpoint_label <- ENDPOINT_LABEL[[key]]
  D$axis <- ENDPOINT_AXIS[[key]]
  D$pathogen <- ENDPOINT_PATHOGEN[[key]]
  D$n_cases <- ENDPOINT_CASES[[key]]
  D$n_controls <- ENDPOINT_CONTROLS[[key]]
  D$OR <- exp(D$b); D$OR_lo <- exp(D$lo); D$OR_hi <- exp(D$hi)
  D <- D[order(D$p), ]
  cnt <- table(O$align)
  ALIGN_LOG[[key]] <<- data.frame(
    endpoint = key, n_rows = nrow(O), n_kept = nrow(M),
    n_palin_ambig = ifelse(is.na(cnt["palin_ambig"]), 0L, as.integer(cnt["palin_ambig"])),
    n_mismatch = ifelse(is.na(cnt["mismatch"]), 0L, as.integer(cnt["mismatch"])),
    stringsAsFactors = FALSE)
  data.table::fwrite(D, file.path(OUT, sprintf("mr_clumped_%s.csv", key)))
  cat(sprintf("[03c] %-12s      %d SNP    palin_ambig %d  taxa %d    p=%.3g %s \n",
              key, nrow(M), ALIGN_LOG[[key]]$n_palin_ambig[1], nrow(D), min(D$p), D$name[1]))
  D
}

ALIGN_LOG <- list()
all_res <- list()
for (k in names(ENDPOINTS)) {
  r <- run_endpoint(k)
  if (!is.null(r)) all_res[[k]] <- r
}
BIG <- do.call(rbind, all_res)
if (is.null(BIG) || !nrow(BIG)) stop("[03c]      ")

BIG$p_bonf <- pmin(1, BIG$p * nrow(BIG))
BIG$p_bh <- p.adjust(BIG$p, method = "BH")
BIG$sig_bonf <- BIG$p_bonf < 0.05
BIG$sig_bh <- BIG$p_bh < 0.05
BIG$p_bh_within_endpoint <- NA_real_
for (k in unique(BIG$endpoint)) {
  idx <- BIG$endpoint == k
  BIG$p_bh_within_endpoint[idx] <- p.adjust(BIG$p[idx], method = "BH")
}
data.table::fwrite(BIG, file.path(OUT, "mr_all_endpoints_clumped.csv"))

ALIGN <- do.call(rbind, ALIGN_LOG)
ALIGN$definition <- "clumped"
data.table::fwrite(ALIGN, file.path(OUT, "align_audit_clumped.csv"))

cat(sprintf("\n[03c] clumped        %d %d taxa × %d    \n",
            nrow(BIG), length(unique(BIG$taxon)), length(all_res)))
cat(sprintf("[03c]    palin_ambig    %d  \n", sum(ALIGN$n_palin_ambig)))
cat(sprintf("[03c] Bonferroni    %d BH FDR    %d\n",
            sum(BIG$sig_bonf), sum(BIG$sig_bh)))
cat(sprintf("[03c]       SEPSIS   BH    %d\n",
            sum(BIG$sig_bh & BIG$endpoint == PRIMARY_ENDPOINT)))
if (sum(BIG$sig_bh)) {
  s <- BIG[BIG$sig_bh, c("name", "endpoint", "n_snp", "b", "se", "OR", "p", "p_bh")]
  print(head(s[order(s$p_bh), ], 40), row.names = FALSE)
}

P1 <- file.path(OUT, "mr_all_endpoints.csv")
if (file.exists(P1)) {
  A <- data.table::fread(P1, data.table = FALSE, showProgress = FALSE)
  A <- A[, c("taxon", "endpoint", "b", "se", "p", "n_snp")]
  B <- BIG[, c("taxon", "endpoint", "b", "se", "p", "n_snp")]
  M <- merge(A, B, by = c("taxon", "endpoint"), suffixes = c("_raw", "_clumped"))
  M$d_b <- M$b_clumped - M$b_raw
  M$ratio_se <- M$se_clumped / M$se_raw
  data.table::fwrite(M, file.path(OUT, "period_compare.csv"))
  cat(sprintf("\n[03c]          %d SE      %.2f |Δb|    %.4f\n",
              nrow(M), median(M$ratio_se, na.rm = TRUE),
              median(abs(M$d_b), na.rm = TRUE)))
}
