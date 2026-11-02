# ============================================================
# ------------------------------------------------------------------
#
#     se(R2) = 4 p (1-p) |beta| se(beta) z = diff/sqrt(se_exp^2*N_exp^2 + se_liab^2*N_out^2)
#
#       work/output/sens_steiger_summary.csv
# ============================================================
.libPaths("D:/Rlibs/4.6")
source("D:/Horzontal_project/20B3/work/course_config.R")
suppressPackageStartupMessages({ library(data.table) })

EXTRACT_DIR <- file.path(DATA, "outcome_extract_clumped")
INS <- data.table::fread(file.path(OUT, "instruments_clumped.csv"),
                         data.table = FALSE, showProgress = FALSE)
COMPL <- c(A = "T", T = "A", C = "G", G = "C")
is_pal <- function(a, b) !is.na(a) && !is.na(b) && a == COMPL[[b]]

rows <- list()
for (key in names(ENDPOINTS)) {
  f <- file.path(EXTRACT_DIR, sprintf("outcome_extract_%s.tsv", key))
  if (!file.exists(f)) next
  O <- data.table::fread(f, data.table = FALSE, showProgress = FALSE)
  O$beta_out <- suppressWarnings(as.numeric(O$beta_out))
  O$se_out <- suppressWarnings(as.numeric(O$se_out))
  O$eaf_out <- suppressWarnings(as.numeric(O$eaf_out))
  ncase <- ENDPOINT_CASES[[key]]; ncont <- ENDPOINT_CONTROLS[[key]]
  ntot <- ncase + ncont
  K <- ncase / ntot
  zK <- qnorm(1 - K)
  phi <- dnorm(zK)
  conv <- K * (1 - K) / phi^2

  for (i in seq_len(nrow(O))) {
    ea <- toupper(O$ea_instr[i]); oa <- toupper(O$oa_instr[i])
    ref <- toupper(O$ref_out[i]); alt <- toupper(O$alt_out[i])
    if (any(is.na(c(ea, oa, ref, alt)))) next
    if (is_pal(ea, oa) && !is.na(O$eaf_out[i]) && abs(O$eaf_out[i] - 0.5) < 0.05) next
    ok <- (ea == alt && oa == ref) || (ea == ref && oa == alt) ||
      (is_pal(ea, oa))
    if (!ok) next
    b_out <- if (ea == ref && oa == alt) -O$beta_out[i] else O$beta_out[i]
    p <- O$eaf_out[i]
    if (is.na(p) || p <= 0 || p >= 1) next

    b_exp <- as.numeric(O$beta_instr[i]); se_exp <- as.numeric(O$se_instr[i])
    se_out <- O$se_out[i]
    if (is.na(b_exp) || is.na(se_exp) || se_exp <= 0 || is.na(se_out) || se_out <= 0) next
    N_exp <- NA_real_
    rows[[length(rows) + 1]] <- data.frame(
      endpoint = key, taxon = O$taxon[i], name = O$name[i], rsid = O$rsid[i],
      eaf = p, beta_exp = b_exp, se_exp = se_exp,
      beta_out = b_out, se_out = se_out, N_exp = N_exp, N_out = ntot, K = K,
      R2_exp = 2 * p * (1 - p) * b_exp^2,
      R2_out_obs = 2 * p * (1 - p) * b_out^2,
      conv = conv, stringsAsFactors = FALSE)
  }
}
D <- do.call(rbind, rows)
if (is.null(D) || !nrow(D)) stop("[04f]       ")

Nmap <- setNames(as.numeric(INS$N), INS$taxon)
D$N_exp <- as.numeric(Nmap[D$taxon])
D <- D[!is.na(D$N_exp) & !is.na(D$N_out), ]
D$R2_out_liab <- D$R2_out_obs * D$conv
D$se_R2_exp <- 4 * D$eaf * (1 - D$eaf) * abs(D$beta_exp) * D$se_exp
D$se_R2_out_liab <- 4 * D$eaf * (1 - D$eaf) * abs(D$beta_out) * D$se_out * D$conv
D$diff <- D$R2_exp * D$N_exp - D$R2_out_liab * D$N_out
D$se_diff <- sqrt((D$se_R2_exp * D$N_exp)^2 + (D$se_R2_out_liab * D$N_out)^2)
D$z <- D$diff / D$se_diff
D$p_steiger <- pnorm(-abs(D$z))
D$p_steiger_one_sided <- ifelse(D$diff > 0, pnorm(-D$z), pnorm(D$z))
D$direction <- ifelse(D$diff > 0 & D$p_steiger_one_sided < 0.05,
                      "tool_to_outcome",
                      ifelse(D$diff < 0 & D$p_steiger_one_sided < 0.05,
                             "outcome_to_tool", "indeterminate"))
D <- D[order(D$p_steiger_one_sided), ]
data.table::fwrite(D, file.path(OUT, "sens_steiger_real.csv"))

cat(sprintf("\n[04f] Steiger   %d taxa %d ×    %d \n",
            nrow(D), length(unique(D$taxon)), length(unique(D$endpoint))))
print(table(D$direction))
S <- data.frame(
  metric = c("n_steiger_rows", "n_taxa",
             "tool_to_outcome", "outcome_to_tool", "indeterminate",
             "share_tool_to_outcome"),
  value = c(nrow(D), length(unique(D$taxon)),
            sum(D$direction == "tool_to_outcome"),
            sum(D$direction == "outcome_to_tool"),
            sum(D$direction == "indeterminate"),
            sprintf("%.3f", mean(D$direction == "tool_to_outcome"))))
data.table::fwrite(S, file.path(OUT, "sens_steiger_summary.csv"), sep = "\t")
print(S, row.names = FALSE)

B <- data.table::fread(file.path(OUT, "mr_all_endpoints_clumped.csv"),
                       data.table = FALSE, showProgress = FALSE)
S2 <- B[B$sig_bh, c("taxon", "endpoint", "name", "b", "p", "p_bh")]
K2 <- merge(S2, D[, c("taxon", "endpoint", "rsid", "R2_exp", "R2_out_liab",
                      "N_exp", "N_out", "z", "p_steiger_one_sided", "direction")],
            by = c("taxon", "endpoint"), all.x = TRUE)
data.table::fwrite(K2, file.path(OUT, "sens_steiger_retained.csv"))
cat("\n[04f]    BH    pair        \n")
print(K2[, c("name", "endpoint", "rsid", "R2_exp", "R2_out_liab",
             "p_steiger_one_sided", "direction")], row.names = FALSE)
cat("[04f] -> sens_steiger_real.csv / sens_steiger_summary.csv / sens_steiger_retained.csv\n")
