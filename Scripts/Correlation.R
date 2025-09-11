###############################################################################
# Title: Cross-dataset correlation and heatmap with significance overlay
# Author: Yumi Kim
#
# Purpose
#   - Read two feature matrices (e.g., metabolites from two tissues/platforms)
#   - Autoscale each matrix (center/scale) for comparability
#   - Compute pairwise Spearman correlations between all features across matrices
#   - Save coefficient/p-value table and draw a heatmap with significance marks
#
# Inputs (CSV files; first column is an ID to be moved to rownames)
#   - "Example data\\Cor_C_Example.csv"
#   - "Example data\\Cor_S_Example.csv"
#
# Outputs
#   - "Coefficient_pvalue.csv" (long-format table: var1, var2, coef, pvalue)
#   - "Cor_Heatmap.pdf"        (heatmap of correlation coefficients, stars for p-values)
#
# Dependencies
#   - Custom functions expected in "R/Cor_functions.R":
#       * NT(data) : normality test per column (writes CSV/PDF)
#       * SC(d1,d2): Spearman correlation between data frames; returns list(coef, p)
#   - Packages used here: mdatools, reshape2, circlize, ComplexHeatmap, grid
#     (The script also loads survival, mediation, variancePartition, ggplot2,
#      ggpattern which may be required by other functions in Cor_functions.R.
#      They are not directly used below.)
###############################################################################

# Clean environment
rm(list = ls())

# Source user-defined functions: NT(), SC(), etc.
source("R/Cor_functions.R")

# ---- R packages (load only; install beforehand if needed) -------------------
library(mdatools)         # prep.autoscale()
library(reshape2)         # melt()
library(circlize)         # colorRamp2()
library(ComplexHeatmap)   # Heatmap(), max_text_width()
library(grid)             # grid.rect(), grid.text(), unit(), gpar()
library(survival)
library(mediation)
library(variancePartition)
library(ggplot2)
library(ggpattern)

# ---- Load data --------------------------------------------------------------
# Expectation: first column is an ID; remaining columns are numeric features.
meta.d.1 <- read.csv("Example data\\Cor_C_Example.csv")
rownames(meta.d.1) <- meta.d.1[, 1]      # move first column to rownames
meta.d.1 <- meta.d.1[, -1]
meta.d.1 <- as.matrix(meta.d.1)

meta.d.2 <- read.csv("Example data\\Cor_S_Example.csv")
rownames(meta.d.2) <- meta.d.2[, 1]
meta.d.2 <- meta.d.2[, -1]
meta.d.2 <- as.matrix(meta.d.2)

# (Optional) keep original column names for reporting
meta.d.realname    <- read.csv("Example data\\Cor_C_Example.csv")
meta_realnames     <- colnames(meta.d.realname[, 2:ncol(meta.d.realname)])
meta.d.2.realname  <- read.csv("Example data\\Cor_S_Example.csv")
meta.2_realnames   <- colnames(meta.d.2.realname[, 2:ncol(meta.d.2.realname)])

# ---- Autoscale each matrix (mean-center and unit-scale) ---------------------
# This helps make features comparable across platforms/units.
meta.d   <- prep.autoscale(meta.d.1, center = TRUE, scale = TRUE) |> as.data.frame()
meta.d.2 <- prep.autoscale(meta.d.2, center = TRUE, scale = TRUE) |> as.data.frame()

# ---- Normality diagnostics (optional; writes CSV/PDF per NT() implementation)
NT(meta.d)
NT(meta.d.2)

# ---- Spearman correlation across the two datasets ---------------------------
# SC() should return: list(correlation_matrix, pvalue_matrix)
SCresult <- SC(meta.d, meta.d.2)
SCmatrix <- SCresult[[1]]
pvalues  <- SCresult[[2]]

# ---- Save long-format table of coefficients and p-values --------------------
melt_SCmatrix <- melt(SCmatrix)
colnames(melt_SCmatrix) <- c("var1", "var2", "coef")
melt_SCpval <- melt(pvalues)
colnames(melt_SCpval) <- c("var1", "var2", "pvalue")
melt_matrix <- cbind(melt_SCmatrix, melt_SCpval)
write.csv(melt_matrix, "Coefficient_pvalue.csv", row.names = FALSE)

# Keep data.frames for plotting
SCmatrix <- as.data.frame(SCmatrix)
SCpvalue <- as.data.frame(pvalues)

# ---- Heatmap with significance overlay --------------------------------------
# Color scale: blue (−1) → white (0) → red (+1)
col_rnorm <- colorRamp2(c(-1, 0, 1), c("#799fcb", "white", "#f9665e"))

pdf(file = "Cor_Heatmap.pdf", width = 10, height = 10)
Heatmap(
  as.matrix(SCmatrix),
  column_title = NULL,
  row_title    = NULL,
  col          = col_rnorm,
  column_dend_side = "top",
  width        = unit(5, "cm"),
  height       = unit(5, "cm"),
  row_names_side    = "left",
  column_names_side = "bottom",
  cluster_rows      = TRUE,
  cluster_columns   = TRUE,
  column_names_gp   = gpar(fontsize = 6),
  row_names_gp      = gpar(fontsize = 6),
  # max_text_width() helps avoid clipping long labels
  row_names_max_width    = max_text_width(rownames(t(SCmatrix))),
  column_names_max_height = max_text_width(rownames(t(SCmatrix))),
  # cell_fun: draw grid border and significance stars based on p-values
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.rect(x = x, y = y, width = w, height = h, gp = gpar(col = "grey", fill = NA))
    if (SCpvalue[i, j] < 0.01) {
      grid.text("**", x, y, vjust = 0.8, gp = gpar(col = "white", fontsize = 6), rot = 90)
    } else if (SCpvalue[i, j] <= 0.05) {
      grid.text("*",  x, y, vjust = 0.8, gp = gpar(col = "white", fontsize = 6), rot = 90)
    }
  }
)
dev.off()

# End of script
