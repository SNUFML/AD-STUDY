# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Yumi Kim

###############################################################################
# Main Script (updated to new function names)
###############################################################################

rm(list = ls())

# ---- Packages ----
suppressPackageStartupMessages({
  library(circlize)
  library(ComplexHeatmap)
  library(grid)
  library(survival)
  library(mediation)
  library(variancePartition)
  library(ggplot2)
  library(ggpattern)
})

# ---- Load data ----
read_feature_matrix <- function(path) {
  df <- read.csv(path, check.names = FALSE)
  rownames(df) <- df[[1]]
  as.matrix(df[-1])
}
meta.d.1 <- read_feature_matrix("Example data\\Cor_C_Example.csv")
meta.d.2 <- read_feature_matrix("Example data\\Cor_S_Example.csv")

# Optional original names (unchanged behavior)
meta_realnames   <- colnames(read.csv("Example data\\Cor_C_Example.csv", check.names = FALSE))[-1]
meta.2_realnames <- colnames(read.csv("Example data\\Cor_S_Example.csv", check.names = FALSE))[-1]

# ---- Autoscale ----
autoscale_df <- function(mat) as.data.frame(scale(mat, TRUE, TRUE), check.names = FALSE)
meta.d   <- autoscale_df(meta.d.1)
meta.d.2 <- autoscale_df(meta.d.2)

# ---- Normality diagnostics ----
run_normality_tests(meta.d)
run_normality_tests(meta.d.2)

# ---- Spearman correlation across datasets ----
C.Result <- spearman_crosscorr(meta.d, meta.d.2)
SCmatrix <- C.Result[[1]]
pvalues  <- C.Result[[2]]

# ---- Save long-format coef/pvalue table ----
df_coef <- .melt_mat(SCmatrix); names(df_coef) <- c("var1", "var2", "coef")
df_pval <- .melt_mat(pvalues);  names(df_pval) <- c("Var1", "Var2", "pvalue")
melt_matrix <- merge(df_coef, df_pval, by.x = c("var1", "var2"),
                     by.y = c("Var1", "Var2"), sort = FALSE)
write.csv(melt_matrix, "Coefficient_pvalue.csv", row.names = FALSE)

# ---- Heatmap with significance overlay ----
col_rnorm <- circlize::colorRamp2(c(-1, 0, 1), c("#799fcb", "white", "#f9665e"))
star_mat <- matrix("", nrow = nrow(pvalues), ncol = ncol(pvalues),
                   dimnames = dimnames(pvalues))
star_mat[pvalues < 0.01] <- "**"
star_mat[pvalues <= 0.05 & pvalues >= 0.01] <- "*"

pdf(file = "Cor_Heatmap.pdf", width = 10, height = 10)
ComplexHeatmap::Heatmap(
  as.matrix(SCmatrix),
  column_title = NULL, row_title = NULL, col = col_rnorm,
  column_dend_side = "top",
  width  = grid::unit(5, "cm"),
  height = grid::unit(5, "cm"),
  row_names_side    = "left",
  column_names_side = "bottom",
  cluster_rows      = TRUE,
  cluster_columns   = TRUE,
  column_names_gp   = grid::gpar(fontsize = 6),
  row_names_gp      = grid::gpar(fontsize = 6),
  row_names_max_width     = ComplexHeatmap::max_text_width(rownames(SCmatrix)),
  column_names_max_height = ComplexHeatmap::max_text_width(colnames(SCmatrix)),
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid::grid.rect(x = x, y = y, width = w, height = h,
                    gp = grid::gpar(col = "grey", fill = NA))
    if (nzchar(star_mat[i, j])) {
      grid::grid.text(star_mat[i, j], x, y, vjust = 0.8,
                      gp = grid::gpar(col = "white", fontsize = 6), rot = 90)
    }
  }
)
dev.off()


# End of script
