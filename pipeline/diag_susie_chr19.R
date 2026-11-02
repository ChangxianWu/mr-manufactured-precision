.libPaths("D:/Rlibs/4.6")
suppressPackageStartupMessages({ library(data.table); library(coloc) })
source("D:/Horzontal_project/20B3/work/course_config.R")
L <- data.table::fread("D:/Horzontal_project/20B3/work/cache/ld/ld_chr19.tsv", data.table = FALSE)
ids <- L$rsid
M <- as.matrix(L[, -(1:4)])
M[is.na(M)] <- 0
R <- suppressWarnings(cor(scale(t(M))))
R[is.na(R)] <- 0
diag(R) <- 1
rownames(R) <- colnames(R) <- ids
E <- data.table::fread("D:/Horzontal_project/20B3/work/data/coloc/GCST90017066_chr19_region.tsv", data.table = FALSE)
E <- E[E$chromosome == 19 & abs(E$base_pair_location - 49218060) <= 100000, ]
E <- E[!duplicated(E$variant_id), ]
O <- data.table::fread("D:/Horzontal_project/20B3/work/data/coloc/outcome_GASTROENT_chr19_region.tsv", data.table = FALSE)
O <- O[abs(as.integer(O$pos) - 48714803) <= 100000, ]
rows <- list()
for (i in seq_len(nrow(O))) {
  rs <- trimws(strsplit(O$rsids[i], ",")[[1]])
  rs <- rs[grepl("^rs[0-9]+$", rs)]
  for (r in rs) rows[[length(rows) + 1]] <- data.frame(
    rsid = r, beta = as.numeric(O$beta[i]), se = as.numeric(O$sebeta[i]),
    eaf = as.numeric(O$af_alt[i]))
}
D <- do.call(rbind, rows)
D <- D[!is.na(D$se) & D$se > 0 & !is.na(D$eaf) & D$eaf > 0 & D$eaf < 1, ]
D <- D[!duplicated(D$rsid), ]
M2 <- merge(data.frame(rsid = E$variant_id, be = as.numeric(E$beta),
                       se_e = as.numeric(E$standard_error), n = as.numeric(E$n)),
            D, by = "rsid")
M2$maf <- pmin(M2$eaf, 1 - M2$eaf)
M2 <- M2[M2$maf > 0.01, ]
ntot <- ENDPOINT_CASES[["GASTROENT"]] + ENDPOINT_CONTROLS[["GASTROENT"]]
keep <- intersect(M2$rsid, rownames(R))
cat("shared with LD:", length(keep), "\n")
Rk <- R[keep, keep]
Mk <- M2[match(keep, M2$rsid), ]
cat("min eig:", min(eigen(Rk, symmetric = TRUE, only.values = TRUE)$values), "\n")
s1 <- list(beta = Mk$be, varbeta = Mk$se_e^2, snp = Mk$rsid, MAF = Mk$maf,
           N = max(Mk$n), type = "quant", LD = Rk)
s2 <- list(beta = Mk$beta, varbeta = Mk$se^2, snp = Mk$rsid, MAF = Mk$maf,
           N = ntot, s = ENDPOINT_CASES[["GASTROENT"]] / ntot, type = "cc", LD = Rk)
r <- coloc.susie(s1, s2, p12 = 1e-5)
cat("--- susie summary ---\n")
print(r$summary)
cat("--- cond1 ---\n")
if (!is.null(r$results$cond1)) print(summary(r$results$cond1))
cat("--- cond2 ---\n")
if (!is.null(r$results$cond2)) print(summary(r$results$cond2))
