# ============================================================
# ------------------------------------------------------------------
#
#
#
# ============================================================
.libPaths("D:/Rlibs/4.6")
source("D:/Horzontal_project/20B3/work/course_config.R")
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

M <- data.table::fread(file.path(OUT, "sens_inflation.csv"), data.table = FALSE,
                       showProgress = FALSE)
M <- M[M$n_snp_raw >= 2, ]
M$k <- M$n_snp_raw
M$pred <- sqrt(M$k)
M$k_eff <- M$se_ratio^2
M$k_eff_ratio <- M$k_eff / M$k
M$axis <- ifelse(M$in_lct == 1, "lactase-axis projection", "non-lactase taxon")

data.table::fwrite(M[, c("taxon", "name", "endpoint", "axis", "k",
                         "se_ratio", "pred", "k_eff", "k_eff_ratio", "n_snp_bin")],
                   file.path(OUT, "inflation_law.csv"), sep = "\t")

bins <- do.call(rbind, lapply(split(M, M$n_snp_bin), function(g) {
  data.frame(bin = g$n_snp_bin[1], n_pairs = nrow(g),
             nominal_k_median = median(g$k),
             observed_ratio_median = median(g$se_ratio),
             predicted_ratio_median = median(g$pred),
             obs_over_pred_median = median(g$se_ratio / g$pred),
             k_eff_median = median(g$k_eff),
             k_eff_over_nominal = median(g$k_eff / g$k))
}))
bins <- bins[order(match(bins$bin, c("2-3", "4-9", "10-99", ">=100"))), ]
data.table::fwrite(bins, file.path(OUT, "inflation_law_bins.csv"), sep = "\t")

byaxis <- do.call(rbind, lapply(split(M, M$axis), function(g) {
  data.frame(axis = g$axis[1], n_pairs = nrow(g),
             k_eff_over_nominal_median = median(g$k_eff_ratio),
             k_eff_max = max(g$k_eff), nominal_k_max = max(g$k))
}))
data.table::fwrite(byaxis, file.path(OUT, "inflation_law_byaxis.csv"), sep = "\t")

fit <- lm(se_ratio ~ pred, data = M)
r2 <- summary(fit)$r.squared

small <- M$k <= 9; large <- M$k >= 100
S <- data.frame(
  metric = c("n_pairs_k_ge2", "slope_obs_vs_sqrtpred", "r2_obs_vs_sqrtpred",
             "small_k_pairs", "small_k_eff_over_nom_median",
             "large_k_pairs", "large_k_eff_over_nom_median",
             "lactase_k_eff_over_nom_median", "nonlactase_k_eff_over_nom_median",
             "lactase_k_eff_max", "lactase_nominal_k_max",
             "large_k_obs_over_pred_median"),
  value = c(nrow(M), round(coef(fit)[2], 4), round(r2, 4),
            sum(small), round(median(M$k_eff_ratio[small]), 3),
            sum(large), round(median(M$k_eff_ratio[large]), 3),
            round(median(M$k_eff_ratio[M$axis == "lactase-axis projection"]), 3),
            round(median(M$k_eff_ratio[M$axis == "non-lactase taxon"]), 3),
            round(max(M$k_eff[M$axis == "lactase-axis projection"]), 1),
            max(M$k[M$axis == "lactase-axis projection"]),
            round(median(M$se_ratio[large] / M$pred[large]), 3)))
data.table::fwrite(S, file.path(OUT, "inflation_law_summary.csv"), sep = "\t")
print(bins, row.names = FALSE)
print(byaxis, row.names = FALSE)
print(S, row.names = FALSE)

mk_panel <- function(D, xlab, ylab, tag, note_at) {
  ggplot(D, aes(k, se_ratio, colour = axis)) +
    geom_function(fun = sqrt, linetype = 2, colour = "grey35", linewidth = 0.8) +
    geom_point(alpha = 0.55, size = 1.5) +
    annotate("text", x = note_at[1], y = note_at[2],
             label = "perfect redundancy: SE ratio = sqrt(k)",
             hjust = 1, size = 3.0, colour = "grey30") +
    scale_colour_manual(values = c("lactase-axis projection" = "#cb181d",
                                   "non-lactase taxon" = "#08519c"), name = NULL) +
    scale_x_log10(limits = c(1.8, 380), breaks = c(2, 3, 10, 30, 100, 300)) +
    scale_y_log10(limits = c(1.2, 26), breaks = c(1.5, 2, 3, 5, 7, 10, 15, 20)) +
    labs(x = xlab, y = ylab, tag = tag) +
    theme_bw(base_size = 11) +
    theme(legend.position = c(0.30, 0.12),
          legend.background = element_rect(fill = "white", colour = NA),
          legend.key.size = unit(0.4, "cm"),
          legend.text = element_text(size = 8),
          plot.tag = element_text(face = "bold", size = 12))
}

pA <- mk_panel(M, "Pre-clumping instruments per taxon (k)",
               "SE ratio (clumped / unclumped)", "A", c(320, 21))

EV <- file.path(OUT, "law_external_validation.csv")
if (file.exists(EV)) {
  E <- data.table::fread(EV, data.table = FALSE, showProgress = FALSE)
  E$k <- E$k_finnish_matched
  E$se_ratio <- E$se_ratio_observed
  E$axis <- ifelse(E$in_lct == 1, "lactase-axis projection", "non-lactase taxon")
  pB <- mk_panel(E, "Matched instruments per taxon (k)",
                 "SE ratio (clumped / unclumped)", "B", c(330, 21))
  p <- patchwork::wrap_plots(pA, pB, ncol = 2) +
    patchwork::plot_annotation(
      title = "Manufactured precision follows the perfect-redundancy bound, then saturates",
      theme = theme(plot.title = element_text(face = "bold", size = 11)))
  ggsave(file.path(FIG, "fig4_law.png"), p, width = 10, height = 4.6, dpi = 300)
} else {
  p <- pA + labs(title = "Manufactured precision follows the perfect-redundancy bound, then saturates") +
    theme(plot.title = element_text(face = "bold", size = 11))
  ggsave(file.path(FIG, "fig4_law.png"), p, width = 7, height = 4.6, dpi = 300)
}
cat("[04h] -> inflation_law*.csv + fig4_law.png\n")
