"""MCP tools for TRADEtools - Transcriptome-wide Analysis of Differential Expression"""

import subprocess
import tempfile
from pathlib import Path
from typing import Annotated, Optional
import pandas as pd
from mcp.server.fastmcp import FastMCP

mcp = FastMCP("tradetools-intro")

# Point to the R scripts directory for this tutorial
R_SCRIPT_DIR = Path(__file__).parent.parent / "r_scripts" / "tradetools_intro"


@mcp.tool()
def trade_univariate(
    results_csv: Annotated[str,
        "Path to CSV file with DESeq2 differential expression results. "
        "Required columns: log2FoldChange (effect sizes), lfcSE (standard errors), pvalue (unadjusted p-values). "
        "Rownames should contain gene identifiers (e.g., ENSG IDs). "
        "Example: DESeq2::results() output saved as CSV with rownames."],
    annot_table_csv: Annotated[Optional[str],
        "Path to CSV file with binary gene annotations for enrichment analysis. "
        "Format: genes as rows (matching results rownames), annotations as columns (0/1 values). "
        "Example: gene sets where 1 indicates membership in that set."] = None,
    log2FoldChange_col: Annotated[str,
        "Column name for log2FoldChange in results CSV."] = "log2FoldChange",
    lfcSE_col: Annotated[str,
        "Column name for log2FoldChange standard errors in results CSV."] = "lfcSE",
    pvalue_col: Annotated[str,
        "Column name for unadjusted p-values in results CSV."] = "pvalue",
    model_significant: Annotated[bool,
        "Whether to model significant genes separately and compute fraction of signal in significant genes."] = True,
    genes_exclude: Annotated[Optional[str],
        "Comma-separated list of gene IDs to exclude from analysis (e.g., perturbed genes themselves)."] = None,
    n_sample: Annotated[int,
        "Number of samples to draw from the inferred effect size distribution (0 = no sampling)."] = 0,
    seed: Annotated[int,
        "Random seed for reproducibility."] = 42,
) -> dict:
    """
    Run univariate TRADE analysis to estimate transcriptome-wide impact of a perturbation.

    TRADE (Transcriptome-wide Analysis of Differential Expression) infers the distribution of
    differential expression effects and estimates:
    - Transcriptome-wide impact: variance of effect sizes (log2FC^2)
    - Me: effective number of DEGs (based on kurtosis)
    - Enrichments: differential expression signal enrichment in gene sets

    The method uses adaptive shrinkage (ashr) to model the effect size distribution as a
    mixture of half-uniform components, accounting for measurement uncertainty.
    """
    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as f:
        output_csv = f.name

    cmd = [
        "Rscript", str(R_SCRIPT_DIR / "trade_univariate.R"),
        "--results", results_csv,
        "--log2FoldChange", log2FoldChange_col,
        "--lfcSE", lfcSE_col,
        "--pvalue", pvalue_col,
        "--model_significant", str(model_significant).upper(),
        "--n_sample", str(n_sample),
        "--seed", str(seed),
        "--output", output_csv
    ]

    if annot_table_csv:
        cmd.extend(["--annot_table", annot_table_csv])

    if genes_exclude:
        cmd.extend(["--genes_exclude", genes_exclude])

    subprocess.run(cmd, check=True)

    result_df = pd.read_csv(output_csv)
    result_rds = output_csv.replace(".csv", ".rds")
    Path(output_csv).unlink()

    result = result_df.to_dict(orient="records")[0]

    return {
        "message": "TRADE univariate analysis completed",
        "reference": "https://github.com/SONGDONGYUAN1994/TRADEtools/blob/main/vignettes/TRADEtools-intro.Rmd",
        "transcriptome_wide_impact": float(result["transcriptome_wide_impact"]),
        "Me": float(result["Me"]),
        "mean": float(result["mean"]),
        "result_rds": result_rds,
    }


@mcp.tool()
def trade_bivariate(
    results1_csv: Annotated[str,
        "Path to CSV file with first DESeq2 differential expression results. "
        "Same format requirements as univariate mode."],
    results2_csv: Annotated[str,
        "Path to CSV file with second DESeq2 differential expression results. "
        "Gene naming must match results1_csv exactly."],
    log2FoldChange_col: Annotated[str,
        "Column name for log2FoldChange in results CSVs."] = "log2FoldChange",
    lfcSE_col: Annotated[str,
        "Column name for log2FoldChange standard errors in results CSVs."] = "lfcSE",
    pvalue_col: Annotated[str,
        "Column name for unadjusted p-values in results CSVs."] = "pvalue",
    genes_exclude: Annotated[Optional[str],
        "Comma-separated list of gene IDs to exclude from analysis."] = None,
    estimate_sampling_covariance: Annotated[bool,
        "Whether to estimate sampling covariance for shared samples/controls using mashr."] = False,
    covariance_matrix_set: Annotated[str,
        "Basis set of covariance matrices: 'mash_default', 'adaptive_grid', or 'combined'."] = "combined",
    component_varexplained_threshold: Annotated[float,
        "Variance explained threshold for retaining components in adaptive grid (0-1)."] = 0.0,
    weight_nocorr: Annotated[float,
        "Prior weight on 0-correlation component (1 = no penalty, >1 = penalty on correlation)."] = 1.0,
    n_sample: Annotated[int,
        "Number of samples to draw from the inferred effect size distribution (0 = no sampling)."] = 0,
    seed: Annotated[int,
        "Random seed for reproducibility."] = 42,
) -> dict:
    """
    Run bivariate TRADE analysis to estimate correlation of differential expression effects
    between two perturbations.

    TRADE bivariate mode estimates the joint distribution of effects and computes:
    - TI correlation: transcriptome-wide impact correlation (correlation of true effects)
    - Raw correlation: Pearson correlation of observed log2FoldChanges
    - Covariance/correlation matrices: inferred effect size relationships

    This uses mashr (multivariate adaptive shrinkage) to jointly model the two sets of
    summary statistics, accounting for measurement uncertainty and potential sampling covariance.
    """
    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as f:
        output_csv = f.name

    cmd = [
        "Rscript", str(R_SCRIPT_DIR / "trade_bivariate.R"),
        "--results1", results1_csv,
        "--results2", results2_csv,
        "--log2FoldChange", log2FoldChange_col,
        "--lfcSE", lfcSE_col,
        "--pvalue", pvalue_col,
        "--estimate_sampling_covariance", str(estimate_sampling_covariance).upper(),
        "--covariance_matrix_set", covariance_matrix_set,
        "--component_varexplained_threshold", str(component_varexplained_threshold),
        "--weight_nocorr", str(weight_nocorr),
        "--n_sample", str(n_sample),
        "--seed", str(seed),
        "--output", output_csv
    ]

    if genes_exclude:
        cmd.extend(["--genes_exclude", genes_exclude])

    subprocess.run(cmd, check=True)

    result_df = pd.read_csv(output_csv)
    result_rds = output_csv.replace(".csv", ".rds")
    Path(output_csv).unlink()

    result = result_df.to_dict(orient="records")[0]

    return {
        "message": "TRADE bivariate analysis completed",
        "reference": "https://github.com/SONGDONGYUAN1994/TRADEtools/blob/main/vignettes/TRADEtools-intro.Rmd",
        "TI_correlation": float(result["TI_correlation"]),
        "cor_raw": float(result["cor_raw"]),
        "loglik": float(result["loglik"]),
        "result_rds": result_rds,
    }


if __name__ == "__main__":
    mcp.run()
