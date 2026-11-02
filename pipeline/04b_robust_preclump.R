# ============================================================
# ------------------------------------------------------------------
#
#       work/output/robust_preclump_summary.csv
# ============================================================
.libPaths("D:/Rlibs/4.6")
source("D:/Horzontal_project/20B3/work/course_config.R")
suppressPackageStartupMessages({ library(data.table) })

args <- commandArgs(trailingOnly = TRUE)
PRESSO_ONLY <- "--presso-only" %in% args
EXTRACT_DIR <- file.path(DATA, "outcome_extract")
COMPL <- c(A = "T", T = "A", C = "G", G = "C")
is_pal <- function(a, b) !is.na(a) && !is.na(b) && a == COMPL[[b]]

align_status <- function(ea, oa, ref, alt, beta_out, eaf_out) {
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

wmedian <- function(rat, w, nboot = 1000, seed = 20261012) {
  est <- function(idx) {
    r <- rat[idx]; ww <- w[idx]
    o <- order(r); r <- r[o]; ww <- ww[o]
    cw <- cumsum(ww) / sum(ww)
    k <- which(cw >= 0.5)[1]
    if (k == 1) return(r[1])
    if (abs(cw[k] - 0.5) < 1e-12) return((r[k] + r[k - 1]) / 2)
    r[k - 1] + (r[k] - r[k - 1]) * (0.5 - cw[k - 1]) / (cw[k] - cw[k - 1])
  }
  b <- est(seq_along(rat))
  set.seed(seed)
  bs <- replicate(nboot, est(sample(seq_along(rat), replace = TRUE)))
  list(b = b, se = sd(bs))
}

aligned_taxa <- function(key) {
  f <- file.path(EXTRACT_DIR, sprintf("outcome_extract_%s.tsv", key))
  if (!file.exists(f)) return(NULL)
  O <- data.table::fread(f, data.table = FALSE, showProgress = FALSE)
  O$beta_out <- suppressWarnings(as.numeric(O$beta_out))
  O$se_out <- suppressWarnings(as.numeric(O$se_out))
  O$eaf_out <- suppressWarnings(as.numeric(O$eaf_out))
  st <- mapply(function(a, b, r, al, bo, eo)
    align_status(a, b, r, al, bo, eo), O$ea_instr, O$oa_instr,
    O$ref_out, O$alt_out, O$beta_out, O$eaf_out, SIMPLIFY = FALSE)
  O$beta_al <- vapply(st, function(x) x$beta, numeric(1))
  O$align <- vapply(st, function(x) x$status, character(1))
  O[O$align %in% c("ok", "flip") & !is.na(O$beta_al) &
      !is.na(O$se_out) & O$se_out > 0 &
      !is.na(O$beta_instr) & O$beta_instr != 0, ]
}

if (!PRESSO_ONLY) {
  res <- list()
  for (key in names(ENDPOINTS)) {
    M <- aligned_taxa(key)
    if (is.null(M)) next
    for (g in split(M, M$taxon)) {
      k <- nrow(g)
      if (k < 3) next
      bx <- as.numeric(g$beta_instr); by <- as.numeric(g$beta_al)
      sy <- as.numeric(g$se_out);     sx <- as.numeric(g$se_instr)
      rat <- by / bx
      w <- 1 / (sy^2 / bx^2)
      b_ivw <- sum(w * rat) / sum(w); se_ivw <- sqrt(1 / sum(w))
      q <- sum(w * (rat - b_ivw)^2); qdf <- k - 1
      fit <- lm(by ~ bx, weights = 1 / sy^2)
      co <- coef(fit); vc <- vcov(fit)
      wm <- wmedian(rat, w)
      res[[paste(key, g$taxon[1], sep = "|")]] <- data.frame(
        endpoint = key, taxon = g$taxon[1], name = g$name[1], level = g$level[1],
        n_snp = k, in_lct = g$in_lct[1],
        b_ivw = b_ivw, se_ivw = se_ivw, p_ivw = 2 * pnorm(-abs(b_ivw / se_ivw)),
        Q = q, I2 = max(0, (q - qdf) / q) * 100,
        Q_p = pchisq(q, qdf, lower.tail = FALSE),
        b_egger = co[2], se_egger = sqrt(vc[2, 2]),
        p_egger = 2 * pnorm(-abs(co[2] / sqrt(vc[2, 2]))),
        egger_intercept = co[1], p_intercept = 2 * pnorm(-abs(co[1] / sqrt(vc[1, 1]))),
        b_wmed = wm$b, se_wmed = wm$se,
        p_wmed = 2 * pnorm(-abs(wm$b / wm$se)),
        stringsAsFactors = FALSE)
    }
  }
  D <- do.call(rbind, res)
  D <- D[order(D$p_ivw), ]
  data.table::fwrite(D, file.path(OUT, "robust_preclump.csv"))
  robust_sig <- function(col) sum(D[[col]] < 0.05, na.rm = TRUE)
  both <- sum(D$p_egger < 0.05 & D$p_wmed < 0.05, na.rm = TRUE)
  cat(sprintf("\n[04b]     taxa >=3       %d taxa %d\n",
              nrow(D), length(unique(D$taxon))))
  cat(sprintf("[04b] IVW    %d | Egger %d |       %d |        %d\n",
              robust_sig("p_ivw"), robust_sig("p_egger"), robust_sig("p_wmed"), both))
  cat(sprintf("[04b] Egger      %d / %d\n",
              robust_sig("p_intercept"), nrow(D)))
  cat(sprintf("[04b] IVW     Egger            %d\n",
              sum(D$p_ivw < 0.05 & D$p_egger >= 0.05 & D$p_wmed >= 0.05, na.rm = TRUE)))
  cat(sprintf("[04b]    I² %.1f%%\n", median(D$I2, na.rm = TRUE)))
  S <- data.frame(
    metric = c("n_multi_variant_tests", "n_taxa", "ivw_sig", "egger_sig",
               "wmed_sig", "both_sig", "egger_intercept_sig",
               "ivw_sig_but_not_robust", "median_I2"),
    value = c(nrow(D), length(unique(D$taxon)), robust_sig("p_ivw"),
              robust_sig("p_egger"), robust_sig("p_wmed"), both,
              robust_sig("p_intercept"),
              sum(D$p_ivw < 0.05 & D$p_egger >= 0.05 & D$p_wmed >= 0.05, na.rm = TRUE),
              round(median(D$I2, na.rm = TRUE), 1)))
  data.table::fwrite(S, file.path(OUT, "robust_preclump_summary.csv"), sep = "\t")
}

M <- aligned_taxa(PRIMARY_ENDPOINT)
if (!is.null(M)) {
  pr <- list()
  for (g in split(M, M$taxon)) {
    k <- nrow(g)
    if (k < 3) next
    d <- data.frame(by = as.numeric(g$beta_al), bx = as.numeric(g$beta_instr),
                    sy = as.numeric(g$se_out), sx = as.numeric(g$se_instr))
    r <- tryCatch(suppressWarnings(
      MRPRESSO::mr_presso(BetaOutcome = "by", BetaExposure = "bx",
                          SdOutcome = "sy", SdExposure = "sx",
                          OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                          data = d, NbDistribution = 1000, SignifThreshold = 0.05)),
      error = function(e) NULL)
    if (is.null(r)) next
    main <- r$`Main MR results`
    gt <- r$`MR-PRESSO results`$`Global Test`
    pr[[g$taxon[1]]] <- data.frame(
      endpoint = PRIMARY_ENDPOINT, taxon = g$taxon[1], name = g$name[1],
      n_snp = k, in_lct = g$in_lct[1],
      presso_raw_b = main$`Causal Estimate`[1], presso_raw_p = main$`P-value`[1],
      presso_corrected_b = main$`Causal Estimate`[2],
      presso_corrected_p = main$`P-value`[2],
      global_p = if (!is.null(gt)) gt$Pvalue else NA_real_,
      n_outlier = if (!is.null(gt)) length(gt$`Outliers Indices`) else NA_integer_,
      stringsAsFactors = FALSE)
  }
  P <- do.call(rbind, pr)
  P <- P[order(P$presso_raw_p), ]
  data.table::fwrite(P, file.path(OUT, "robust_preclump_presso.csv"))
  cat(sprintf("\n[04b] MR-PRESSO %s taxa %d raw p<0.05 %d     p<0.05 %d\n",
              PRIMARY_ENDPOINT, nrow(P), sum(P$presso_raw_p < 0.05, na.rm = TRUE),
              sum(P$presso_corrected_p < 0.05, na.rm = TRUE)))
  cat(sprintf("[04b]        %d        SNP   %.0f\n",
              sum(P$global_p < 0.05, na.rm = TRUE), median(P$n_outlier, na.rm = TRUE)))
}
cat("[04b]   \n")
