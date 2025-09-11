###############################################################################
# Title: Visualization workflow — group comparisons, volcano plots, and PCA
# Author: Yumi Kim
# Purpose:
#   - Run two predefined group comparisons on a metabolomics-like dataset
#   - Visualize results as boxplots and volcano plots
#   - Perform PCA, then test group differences on PC scores and visualize them
#
# Inputs:
#   - CSV file: "Example data/Vis_Example.csv"
#       * Column 1: Sample (sample identifier)
#       * Column 2: Group  (group label; used for comparisons and plotting)
#       * Columns 3+: numeric variables/features
#
# Dependencies:
#   - Custom functions sourced from "R/Vis_Function.R":
#       YMSTAT, Boxplot_ym, Volcano_YM, Volcano_YM_S, PCA_YM_fviz, Boxplot_ym_pc
#   - R packages: readxl, dplyr, grid, LMSstat, writexl, devtools, ComplexHeatmap
#     (Your Vis_Function.R may use additional packages internally.)
#

###############################################################################

# 0) Load visualization/stat helper functions
source("R/Vis_Function.R")

# 1) Load required packages (some of these may be used inside Vis_Function.R)
library(readxl)         # reading Excel files (not used directly here, but kept)
library(dplyr)          # data manipulation (used for piping/verbs)
library(grid)           # low-level plotting utilities (used by ggplot/heatmap)
library(LMSstat)        # (if used inside Vis_Function.R)
library(writexl)        # write xlsx (if Vis_Function.R exports results)
library(devtools)       # dev tooling (if Vis_Function.R calls dev helpers)
library(ComplexHeatmap) # complex heatmaps (if used in Vis_Function.R)

# 2) Read data
#    - check.names = FALSE keeps original column names (e.g., "AD^APT+F414")
data <- read.csv("Example data/Vis_Example.csv", check.names = FALSE)

# 3) Run statistics for two group comparisons using YMSTAT
#    - group_combinations: list of two pairwise contrasts
#        * "AD^APT" vs "AD^WT"
#        * "AD^APT+F414" vs "AD^APT"
#    - Adjust_p_value = FALSE: no multiple-testing correction at this stage
Statfile <- YMSTAT(
  data,
  group_combinations = list(
    c("AD^APT", "AD^WT"),
    c("AD^APT+F414", "AD^APT")
  ),
  Adjust_p_value = FALSE,
  Adjust_method   = "BH"
)

# 4) Box plot of variables across groups
#    - asterisk = "u_test": annotate using Mann–Whitney U-test results from YMSTAT
#    - color: manual color palette for the groups (ensure order matches `order=`)
#    - significant_variable_only = FALSE: plot all variables, not only significant
#    - fig_height/fig_width: final device size (if the function saves plots)
Boxplot_ym(
  Statfile,
  asterisk = "u_test",                          # NOTE: If YMSTAT uses "u-test", change this accordingly.
  color    = c("#323232", "#808080", "#a00000"),
  T_size   = 30,
  X_text   = 22,
  Y_text   = 22,
  width    = 0.45,
  significant_variable_only = FALSE,
  label_size = 5.5,
  size       = 0.7,
  fig_height = 9,
  fig_width  = 6,
  order = c("AD^WT", "AD^APT", "AD^APT+F414"), # Ensure these match data$Group exactly
  x_text_bold = "bold",
  x_text_angle = 45,
  y_text_bold  = "bold"
)

# 5) Volcano plots
#    - `comb` defines which contrasts to visualize
#    - `asterisk` MUST match the YMSTAT column name for U-test:
#        If YMSTAT outputs "u_test", use "u_test"; if "u-test", use "u-test".
Volcano_YM(
  Statfile,
  asterisk = "u_test",                          # <-- check your YMSTAT output column name
  comb = list(
    c("AD^APT",       "AD^WT"),
    c("AD^APT+F414",  "AD^APT")
  )
)

#    - Significant-only volcano (e.g., highlights variables passing a threshold)
Volcano_YM_S(
  Statfile,
  asterisk = "u_test",                          # <-- same naming note as above
  comb = list(
    c("AD^APT",       "AD^WT"),
    c("AD^APT+F414",  "AD^APT")
  )
)

# 6) PCA
#    - PCA_YM_fviz: runs PCA and returns an object containing coordinates, etc.
#    - color/order: ensure the group order matches your legend and plotting order
#    - legend_position: where to place the legend on the final figure
pca <- PCA_YM_fviz(
  data,
  color  = c("#323232", "#808080", "#a00000"),
  order  = c("AD^WT", "AD^APT", "AD^APT+F414"),
  title  = "Title",
  legend_position = "bottom"
)

# Extract grouping columns for later merging with PC scores
group <- data[1:2]  # assumes col1 = Sample, col2 = Group

# NOTE: Many PCA objects store coordinates as a data.frame or matrix.
#       Using `[1]` returns only the first element; to extract the full
#       PC1/PC2 columns, you typically want `pca$coordinates[, 1]` and `[, 2]`.
#       Keeping your original lines below for consistency:
PC1 <- pca$coordinates[1]
PC2 <- pca$coordinates[2]

# Build a data.frame with Sample, Group, PC1, PC2
PC_DATA <- cbind(group, PC1, PC2)

# 7) Statistics on PC scores and PC score box plots
#    - Run YMSTAT on PC1/PC2 (two comparisons as above)
Statfile_pc <- YMSTAT(
  PC_DATA,
  group_combinations = list(
    c("AD^APT",       "AD^WT"),
    c("AD^APT+F414",  "AD^APT")
  ),
  Adjust_p_value = FALSE
)

#    - Boxplots of PC scores
#    - asterisk = "t_test": annotate using t-test p-values from YMSTAT
#    - order: Make sure group names here match your data exactly.
Boxplot_ym_pc(
  Statfile_pc,
  asterisk = "t_test",
  color    = c("#323232", "#808080", "#a00000"),
  T_size   = 30,
  X_text   = 22,
  Y_text   = 22,
  width    = 0.45,
  significant_variable_only = FALSE,
  label_size = 5,
  size       = 0.7,
  fig_height = 10,
  fig_width  = 4,
  order = c("AD^WT", "AD^APT", "AD^APT+CKDB001"),  # <-- If your data uses "AD^APT+F414", align this label.
  x_text_bold = "bold",
  x_text_angle = 45,
  y_text_bold  = "bold"
)

###############################################################################
###############################################################################
