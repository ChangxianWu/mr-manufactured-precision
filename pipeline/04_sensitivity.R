# ============================================================
# ------------------------------------------------------------------
#
# ============================================================
.libPaths("D:/Rlibs/4.6")
source("D:/Horzontal_project/20B3/work/course_config.R")
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

RAW <- file.path(OUT, "mr_all_endpoints.csv")
CLM <- file.path(OUT, "mr_all_endpoints_clumped.csv")
A <- data.table::fread(RAW, data.table = FALSE, showProgress = FALSE)
B <- data.table::fread(CLM, data.table = FALSE, showProgress = FALSE)

M <- merge(A[, c("taxon", "name", "endpoint", "axis", "pathogen", "n_snp", "b", "se", "p", "OR", "in_lct")],
           B[, c("taxon", "endpoint", "n_snp", "b", "se", "p", "OR")],
           by = c("taxon", "endpoint"), suffixes = c("_raw", "_clumped"))
M$se_ratio <- M$se_clumped / M$se_raw
M$d_logp <- -log10(M$p_clumped) - (-log10(M$p_raw))
M$n_snp_bin <- cut(M$n_snp_raw, breaks = c(0, 1, 3, 9, 99, 1e6),
                   labels = c("1", "2-3", "4-9", "10-99", ">=100"))
data.table::fwrite(M, file.path(OUT, "sens_inflation.csv"), sep = "\t")

agg <- aggregate(cbind(n_pairs = rep(1, nrow(M)), med_se_ratio = M$se_ratio,
                       med_p_raw = M$p_raw, med_p_clm = M$p_clumped),
                 by = list(n_snp_bin = M$n_snp_bin), FUN = median)
agg$n_pairs <- as.integer(table(M$n_snp_bin)[as.character(agg$n_snp_bin)])
print(agg, row.names = FALSE)

cl <- data.table::fread(file.path(OUT, "instruments_clumped.csv"),
                        data.table = FALSE, showProgress = FALSE)
steig <- do.call(rbind, lapply(seq_len(nrow(cl)), function(i) {
  r <- cl[i, ]
  data.frame(taxon = r$taxon, name = r$name, rsid = r$rsid, n_snp = r$n_snp,
             beta = as.numeric(r$beta), se = as.numeric(r$se),
             F = as.numeric(r$F_min), N_exp = as.numeric(r$N),
             in_lct = as.integer(r$in_lct), stringsAsFactors = FALSE)
}))
sp <- B[B$endpoint == PRIMARY_ENDPOINT, c("taxon", "b", "se", "p")]
steig <- merge(steig, sp, by = "taxon", all.x = TRUE, suffixes = c("", "_out"))
steig$direction <- ifelse(is.na(steig$b), "no_outcome",
                          ifelse(steig$F > 10, "strong_tool", "weak_tool"))
data.table::fwrite(steig, file.path(OUT, "sens_steiger.csv"), sep = "\t")

n_bh_raw <- sum(A$sig_bh); n_bh_clm <- sum(B$sig_bh)
sub_a <- A[A$in_lct == 0, ]; sub_b <- B[B$in_lct == 0, ]
lct <- data.frame(
  layer = c("BH_survivors_all", "BH_survivors_non_LCT",
            "n_tests_all", "n_tests_non_LCT"),
  unclumped = c(n_bh_raw, sum(sub_a$sig_bh), nrow(A), nrow(sub_a)),
  clumped = c(n_bh_clm, sum(sub_b$sig_bh), nrow(B), nrow(sub_b)))
data.table::fwrite(lct, file.path(OUT, "sens_lct_drop.csv"), sep = "\t")
print(lct, row.names = FALSE)

funnel <- data.table::fread(file.path(OUT, "iv_funnel.csv"), data.table = FALSE)
clumpsum <- data.table::fread(file.path(OUT, "clump_summary.csv"), data.table = FALSE)
rep_sum <- file.path(OUT, "replication_summary.csv")
rep_n <- if (file.exists(rep_sum)) nrow(data.table::fread(rep_sum, data.table = FALSE)) else 0

sink(file.path(OUT, "sens_summary.csv"))
cat("metric\tvalue\n")
cat(sprintf("nominal_taxa\t%d\n", funnel$n[funnel$layer == "L0_nominal_taxa"]))
cat(sprintf("instrumentable_taxa_P5e8\t%d\n", funnel$n[funnel$layer == "L2_primary_taxa"]))
cat(sprintf("lct_projected_taxa\t%d\n", funnel$n[funnel$layer == "L2_lct_projected"]))
cat(sprintf("non_lct_taxa\t%d\n", funnel$n[funnel$layer == "L2_non_lct_taxa"]))
cat(sprintf("aligned_to_finnish_or_dmp\t%d\n", funnel$n[funnel$layer == "L3_aligned_taxa"]))
cat(sprintf("snp_before_clump\t%d\n", clumpsum$n_snp_before[1]))
cat(sprintf("snp_after_clump\t%d\n", clumpsum$n_snp_after[1]))
cat(sprintf("clump_inflation_factor\t%.1f\n", as.numeric(clumpsum$inflation_factor[1])))
cat(sprintf("med_se_ratio_clump_vs_raw\t%.3f\n", median(M$se_ratio, na.rm = TRUE)))
cat(sprintf("med_se_ratio_only_multiSNP\t%.3f\n",
            median(M$se_ratio[M$n_snp_raw > 1], na.rm = TRUE)))
cat(sprintf("bh_survivors_raw\t%d\n", n_bh_raw))
cat(sprintf("bh_survivors_clumped\t%d\n", n_bh_clm))
cat(sprintf("bh_survivors_clumped_non_lct\t%d\n", sum(sub_b$sig_bh)))
cat(sprintf("replication_taxa_assessed\t%d\n", rep_n))
cat(sprintf("max_abs_log10p_raw\t%.1f\n", max(-log10(A$p))))
cat(sprintf("max_abs_log10p_clumped\t%.1f\n", max(-log10(B$p))))
COLOC_F <- file.path(OUT, "coloc_focal.csv")
if (file.exists(COLOC_F)) {
  CO <- data.table::fread(COLOC_F, data.table = FALSE)
  cat(sprintf("coloc_loci_assessed\t%d\n", length(unique(CO$locus))))
  cat(sprintf("coloc_pairs_assessed\t%d\n", nrow(CO)))
  cat(sprintf("coloc_colocalised_PPH4ge0.8\t%d\n", sum(CO$verdict == "colocalised")))
  cat(sprintf("coloc_suggestive_PPH4ge0.5\t%d\n", sum(CO$verdict == "suggestive")))
  cat(sprintf("coloc_max_PPH4\t%.3f\n", max(CO$PP_H4)))
  cat(sprintf("coloc_max_PPH4_pair\t%s_%s\n",
              CO$locus[which.max(CO$PP_H4)], CO$endpoint[which.max(CO$PP_H4)]))
  cat(sprintf("coloc_chr2_max_PPH4\t%.3f\n", max(CO$PP_H4[CO$locus == "chr2"])))
  cat(sprintf("coloc_chr19_max_PPH4\t%.3f\n", max(CO$PP_H4[CO$locus == "chr19"])))
}
sink()
cat("\n[04]    -> sens_summary.csv\n")

f1 <- as.data.frame(funnel)
f1$label <- c("Nominal taxa", "Top-hit rows (P<=1e-4)", "Instrumentable, P<5e-8 (>=1 SNP)",
              "Instrumentable, P<5e-8 (>=2 SNP)", "Instrumentable, P<1e-5 (>=1 SNP)",
              "Instrumentable, P<1e-5 (>=2 SNP)", "Passing P<5e-8 (orig)",
              "Primary signal in LCT window", "Independent non-LCT taxa",
              "Non-LCT 5 Mb loci", "Aligned to Finnish/DMP", "Non-LCT & aligned")[seq_len(nrow(f1))]
keep <- c("L0_nominal_taxa", "L2_primary_taxa", "L2_non_lct_taxa", "L3_non_lct_aligned")
f1b <- f1[f1$layer %in% keep, ]
f1b$label <- c("Nominal taxa\n(211)", "Instrumentable P<5e-8\n(27)",
               "Independent non-LCT\n(22)", "Cross-cohort alignable\n(3)")
f1b$label <- factor(f1b$label, levels = f1b$label[order(-f1b$n)])
p1 <- ggplot(f1b, aes(label, n, fill = label)) +
  geom_col(width = 0.6) + geom_text(aes(label = n), vjust = -0.4, size = 4) +
  scale_fill_manual(values = c("#2c7fb8", "#41b6c4", "#7fcdbb", "#c7e9b4")) +
  labs(x = NULL, y = "Number of taxa", title = "Exposure-side funnel (MiBioGen)") +
  theme_bw(base_size = 12) + theme(legend.position = "none",
                                    plot.title = element_text(face = "bold"))
ggsave(file.path(FIG, "fig1_funnel.png"), p1, width = 7, height = 4.5, dpi = 300)


AXLAB <- c("lactase-axis projection", "non-lactase taxon")
CM <- rbind(
  data.frame(set = "All variants at P<5e-8", b = A$b, p = A$p, endpoint = A$endpoint,
             axis = ifelse(A$in_lct == 1, AXLAB[1], AXLAB[2]),
             taxon = A$taxon, name = A$name),
  data.frame(set = "After LD clumping", b = B$b, p = B$p, endpoint = B$endpoint,
             axis = ifelse(B$in_lct == 1, AXLAB[1], AXLAB[2]),
             taxon = B$taxon, name = B$name))
CM$set <- factor(CM$set, levels = c("All variants at P<5e-8", "After LD clumping"))
HL <- CM[CM$set == "After LD clumping" &
           grepl("Ruminococcustorques", CM$name) &
           CM$endpoint %in% c("INTEST_INF", "GASTROENT"), ]
p3 <- ggplot(CM, aes(b, -log10(p), colour = axis)) +
  geom_point(alpha = 0.4, size = 1.1) +
  geom_hline(yintercept = -log10(0.05 / nrow(A)), linetype = 2, colour = "grey40") +
  geom_point(data = HL, shape = 1, size = 3.2, stroke = 1.1, colour = "black") +
  scale_colour_manual(values = setNames(c("#cb181d", "#08519c"), AXLAB), name = NULL) +
  facet_wrap(~ set, ncol = 2) +
  labs(x = "MR effect (log OR per SD)", y = "-log10 P",
       title = "Clumping removes the apparent signal from the lactase axis") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")
ggsave(file.path(FIG, "fig2_volcano.png"), p3, width = 9, height = 4.8, dpi = 300)

if (file.exists(rep_sum)) {
  R <- data.table::fread(rep_sum, data.table = FALSE)
  R$same_sign_rate <- as.numeric(R$same_sign_rate)
  R$name <- sub("^.*\\.", "", R$taxon)
  R <- R[!is.na(R$same_sign_rate) & R$n_snp >= 5, ]
  p4 <- ggplot(R, aes(reorder(name, same_sign_rate), same_sign_rate)) +
    geom_col(fill = "#41ab5d", width = 0.6) +
    geom_hline(yintercept = 0.5, linetype = 2, colour = "grey30") +
    coord_flip(ylim = c(0, 1.05)) +
    labs(x = NULL, y = "Direction agreement with Finnish cohort",
         title = "Cross-cohort exposure effect agreement (MiBioGen vs Finnish)") +
    theme_bw(base_size = 12) + theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(FIG, "fig5_replication.png"), p4, width = 7.5, height = 4.5, dpi = 300)
}

S <- B[B$sig_bh, ]
if (nrow(S)) {
  S$lab <- sprintf("%s\n(%s)", S$name, ENDPOINT_LABEL[S$endpoint])
  S$lab <- factor(S$lab, levels = S$lab[order(S$b)])
  p5 <- ggplot(S, aes(b, lab)) +
    geom_vline(xintercept = 0, colour = "grey50") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15, colour = "#08519c") +
    geom_point(size = 2.6, colour = "#08519c") +
    labs(x = "MR effect (log OR per SD, 95% CI)", y = NULL,
         title = "Pairs surviving BH FDR after clumping") +
    theme_bw(base_size = 12) + theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(FIG, "fig3_forest.png"), p5, width = 8, height = 3.2, dpi = 300)
}
cat("[04]    -> work/figures/ (fig1..fig5)\n")
