#!/usr/bin/env Rscript
library(optparse)
library(data.table)
library(TRADEtools)

option_list <- list(
  make_option("--results", type = "character", help = "Path to CSV file with DESeq2 results (required columns: log2FoldChange, lfcSE, pvalue)"),
  make_option("--annot_table", type = "character", default = "", help = "Path to CSV file with gene annotations (binary matrix, genes as rows)"),
  make_option("--log2FoldChange", type = "character", default = "log2FoldChange", help = "Column name for log2FoldChange [default: %default]"),
  make_option("--lfcSE", type = "character", default = "lfcSE", help = "Column name for log2FoldChange standard errors [default: %default]"),
  make_option("--pvalue", type = "character", default = "pvalue", help = "Column name for p-values [default: %default]"),
  make_option("--model_significant", type = "logical", default = TRUE, help = "Model significant genes separately [default: %default]"),
  make_option("--genes_exclude", type = "character", default = "", help = "Comma-separated list of genes to exclude"),
  make_option("--n_sample", type = "integer", default = 0, help = "Number of samples to draw from distribution [default: 0 (no sampling)]"),
  make_option("--seed", type = "integer", default = 42, help = "Random seed for reproducibility [default: %default]"),
  make_option("--output", type = "character", help = "Output CSV file (required)")
)

args <- parse_args(OptionParser(option_list = option_list))

# Set random seed for reproducibility
set.seed(args$seed)

# Load results
results_df <- as.data.frame(fread(args$results))
if ("rn" %in% colnames(results_df)) {
  rownames(results_df) <- results_df$rn
  results_df$rn <- NULL
}

# Load annotation table if provided
annot_table <- NULL
if (args$annot_table != "") {
  annot_table <- as.data.frame(fread(args$annot_table))
  if ("rn" %in% colnames(annot_table)) {
    rownames(annot_table) <- annot_table$rn
    annot_table$rn <- NULL
  }
}

# Parse genes to exclude
genes_exclude <- NULL
if (args$genes_exclude != "") {
  genes_exclude <- strsplit(args$genes_exclude, ",")[[1]]
}

# Parse n_sample
n_sample <- if (args$n_sample > 0) args$n_sample else NULL

# Run TRADE univariate analysis
trade_result <- TRADE(
  mode = "univariate",
  results1 = results_df,
  annot_table = annot_table,
  log2FoldChange = args$log2FoldChange,
  lfcSE = args$lfcSE,
  pvalue = args$pvalue,
  model_significant = args$model_significant,
  genes_exclude = genes_exclude,
  n_sample = n_sample,
  verbose = FALSE
)

# Extract key results for output
output_data <- list()

# Distribution summary (main results)
if (!is.null(trade_result$distribution_summary)) {
  output_data$transcriptome_wide_impact <- trade_result$distribution_summary$transcriptome_wide_impact
  output_data$Me <- trade_result$distribution_summary$Me
  output_data$mean <- trade_result$distribution_summary$mean
}

# Save full result as RDS
rds_file <- sub("\\.csv$", ".rds", args$output)
saveRDS(trade_result, rds_file)

# Write summary to CSV
result_df <- data.table(
  transcriptome_wide_impact = output_data$transcriptome_wide_impact,
  Me = output_data$Me,
  mean = output_data$mean
)
fwrite(result_df, args$output)
