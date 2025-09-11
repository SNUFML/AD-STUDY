###############################################################################
# Title: Multi-organ split workflow — filter by organ, impute, and run stats
# Author: Yumi Kim
#
# What this script does
#   1) Loads helper functions from `R/MultiOrgan_Total_Split_Function.R`
#      (expected to provide `group_na_summary()` and `YMSTAT()`).
#   2) Reads a wide table where:
#        - Column 1 = Sample ID
#        - Column 2 = Group (e.g., "Feces_AP", "Serum_WP", ...)
#        - Columns 3+ = numeric features (metabolites, etc.)
#   3) Uses `group_na_summary()` to mark variables that pass an NA threshold
#      within each group, then for a set of organ/tissue keywords
#      ("Feces", "Serum", "Cecal", "Cortex", "Hippo"):
#         - Keeps variables with acceptable NA ratio for that organ/tissue
#         - Subsets rows to the organ/tissue of interest (via `Group`)
#         - Simple imputation: replace remaining NAs with (global min of
#           numeric part)/5 within that organ subset
#         - Runs `YMSTAT()` on two pre-specified contrasts and saves results
#
# Inputs
#   - R/MultiOrgan_Total_Split_Function.R
#   - Example data\\MultiOrgan_Example.csv
#
# Outputs (per keyword)
#   - Total_<keyword>.csv          : filtered & imputed data for that organ
#   - <keyword>_statresult.csv     : YMSTAT() result matrix for that organ
#
###############################################################################

# Load helper functions (defines group_na_summary(), YMSTAT(), etc.)
source("R/MultiOrgan_Total_Split_Function.R")

# --- Required packages (loaded only) -----------------------------------------
library(readxl)
library(dplyr)
library(grid)
library(LMSstat)
library(writexl)
library(devtools)
library(ComplexHeatmap)

# --- Read input data ---------------------------------------------------------
# keep original column names (special characters allowed in Group values)
data <- read.csv("Example data\\MultiOrgan_Example.csv", check.names = FALSE)

# Organ/tissue prefixes to iterate over (must match the start of Group values)
filter_keywords <- c("Feces", "Serum", "Cecal", "Cortex", "Hippo")

# Container to collect YMSTAT() results per keyword (optional)
results_list <- list()

# 1) Compute NA-allowance summary by group (provided in sourced file)
#    Expected to return a matrix/data frame allowing us to decide which
#    variables (columns) pass NA threshold within each group.
Data <- group_na_summary(data)

# 2) Loop over each organ/tissue keyword and process
for (keyword in filter_keywords) {
  
  # 2a) Pick rows in `Data` (groups) that start with the keyword (e.g., "Feces")
  filtered_data <- as.data.frame(Data[grepl(paste0("^", keyword), rownames(Data)), ])
  
  # 2b) Keep columns (variables) that are marked TRUE for at least one of the
  #     selected groups (i.e., pass NA threshold in that organ)
  filtered_idx <- which(colSums(filtered_data) != 0)
  
  # 2c) Subset the original `data` to those retained columns
  filtered_only <- data[, filtered_idx]
  
  # 2d) Keep only the rows whose Group starts with the keyword
  filtered_data_final <- filtered_only %>%
    dplyr::filter(grepl(paste0("^", keyword), Group))
  
  # 2e) Simple imputation within the numeric block (cols 3:ncol):
  #     replace remaining NAs with (global min of numeric part)/5
  filtered_data_final[is.na(filtered_data_final)] <-
    min(filtered_data_final[, 3:ncol(filtered_data_final)], na.rm = TRUE) / 5
  
  # 2f) Save the per-organ filtered/imputed table
  write.csv(filtered_data_final, paste0("Total_", keyword, ".csv"), row.names = FALSE)
  
  # 2g) If the organ subset has any rows, run YMSTAT with two pairwise contrasts:
  #     (keyword_AP vs keyword_WP) and (keyword_AS vs keyword_AP)
  if (nrow(filtered_data_final) > 0) {
    result <- YMSTAT(
      filtered_data_final,
      group_combinations = list(
        c(paste0(keyword, "_AP"), paste0(keyword, "_WP")),
        c(paste0(keyword, "_AS"), paste0(keyword, "_AP"))
      ),
      Adjust_p_value = FALSE
    )
    
    # Save the result matrix
    write.csv(result[["Result"]], paste0(keyword, "_statresult.csv"))
    
    # Also keep it in memory (optional)
    results_list[[keyword]] <- result[["Result"]]
  }
}

# End of script
