#!/usr/bin/env Rscript
library(optparse)
library(data.table)
library(TRADEtools)

option_list <- list(
  make_option("--results1", type = "character", help = "Path to CSV file with first DESeq2 results"),
  make_option("--results2", type = "character", help = "Path to CSV file with second DESeq2 results"),
  make_option("--log2FoldChange", type = "character", default = "log2FoldChange", help = "Column name for log2FoldChange [default: %default]"),
  make_option("--lfcSE", type = "character", default = "lfcSE", help = "Column name for log2FoldChange standard errors [default: %default]"),
  make_option("--pvalue", type = "character", default = "pvalue", help = "Column name for p-values [default: %default]"),
  make_option("--genes_exclude", type = "character", default = "", help = "Comma-separated list of genes to exclude"),
  make_option("--estimate_sampling_covariance", type = "logical", default = FALSE, help = "Estimate sampling covariance for shared samples [default: %default]"),
  make_option("--covariance_matrix_set", type = "character", default = "combined", help = "Covariance matrix set: mash_default, adaptive_grid, or combined [default: %default]"),
  make_option("--component_varexplained_threshold", type = "double", default = 0, help = "Variance explained threshold for adaptive grid [default: %default]"),
  make_option("--weight_nocorr", type = "double", default = 1, help = "Prior weight on 0 correlation component [default: %default]"),
  make_option("--n_sample", type = "integer", default = 0, help = "Number of samples to draw from distribution [default: 0 (no sampling)]"),
  make_option("--seed", type = "integer", default = 42, help = "Random seed for reproducibility [default: %default]"),
  make_option("--output", type = "character", help = "Output CSV file (required)")
)

args <- parse_args(OptionParser(option_list = option_list))

# Set random seed for reproducibility
set.seed(args$seed)

# Load results
results1_df <- as.data.frame(fread(args$results1))
if ("rn" %in% colnames(results1_df)) {
  rownames(results1_df) <- results1_df$rn
  results1_df$rn <- NULL
}

results2_df <- as.data.frame(fread(args$results2))
if ("rn" %in% colnames(results2_df)) {
  rownames(results2_df) <- results2_df$rn
  results2_df$rn <- NULL
}

# Parse genes to exclude
genes_exclude <- NULL
if (args$genes_exclude != "") {
  genes_exclude <- strsplit(args$genes_exclude, ",")[[1]]
}

# Parse n_sample
n_sample <- if (args$n_sample > 0) args$n_sample else NULL

# Run TRADE bivariate analysis
trade_result <- TRADE(
  mode = "bivariate",
  results1 = results1_df,
  results2 = results2_df,
  log2FoldChange = args$log2FoldChange,
  lfcSE = args$lfcSE,
  pvalue = args$pvalue,
  genes_exclude = genes_exclude,
  estimate_sampling_covariance = args$estimate_sampling_covariance,
  covariance_matrix_set = args$covariance_matrix_set,
  component_varexplained_threshold = args$component_varexplained_threshold,
  weight_nocorr = args$weight_nocorr,
  n_sample = n_sample,
  verbose = FALSE
)

# Save full result as RDS
rds_file <- sub("\\.csv$", ".rds", args$output)
saveRDS(trade_result, rds_file)

# Write summary to CSV
result_df <- data.table(
  TI_correlation = trade_result$TI_correlation,
  cor_raw = trade_result$cor_raw,
  loglik = trade_result$loglik
)
fwrite(result_df, args$output)
