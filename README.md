# [Paper2Agent](https://github.com/jmiao24/Paper2Agent): TRADEtools Demo

A demonstration of turning the [TRADE paper](https://www.biorxiv.org/content/10.1101/2024.04.26.591331v1) into an interactive AI agent. This project transforms Song et al.'s TRADEtools (Transcriptome-wide Analysis of Differential Expression) method for estimating transcriptome-wide impact of genetic perturbations into a conversational agent that can run univariate and bivariate TRADE analyses through natural language.

## Folder Structure

```
trade_agent/
├── src/
│   ├── TRADEtools_mcp.py         # MCP server entry point
│   ├── requirements.txt          # Python dependencies
│   ├── tools/                    # Python wrappers for R scripts
│   │   └── tradetools_intro.py   # TRADE tool implementations
│   └── r_scripts/                # R scripts for each tool
│       └── tradetools_intro/
│           ├── trade_univariate.R   # Univariate TRADE analysis
│           └── trade_bivariate.R    # Bivariate TRADE analysis
└── README.md
```

## Quick Start

### 1. Clone the Repository

```bash
git clone https://github.com/jmiao24/trade_agent.git
cd trade_agent
```

### 2. Install Gemini CLI

Install the [Google Gemini CLI](https://github.com/google-gemini/gemini-cli):

```bash
brew install gemini-cli
```

### 3. Install R Dependencies

The agent requires R and the `TRADEtools` package. Install R from [CRAN](https://cran.r-project.org/), then install the required packages:

```r
# Install dependencies
install.packages(c("optparse", "data.table", "ashr"))

# Install mashr from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("mashr")

# Install TRADEtools from GitHub
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("SONGDONGYUAN1994/TRADEtools")
```

### 4. Install FastMCP

```bash
pip install fastmcp
```

### 5. Install MCP Server

Install the TRADEtools MCP server using fastmcp:

```bash
fastmcp install gemini-cli ./src/TRADEtools_mcp.py --with-requirements ./src/requirements.txt
```

### 6. Start the Agent

Start Gemini CLI in the repository folder:

```bash
gemini
```

You will now have access to the TRADEtools agent with all available tools.

## Example Query

```
Run a univariate TRADE analysis on my DESeq2 results to estimate the transcriptome-wide
impact of my perturbation. Exclude the perturbed gene from the analysis.
```

## Available Agent Tools

The agent provides the following capabilities through natural language:

### Univariate Analysis
- `trade_univariate`: Estimate transcriptome-wide impact of a single perturbation
  - Uses adaptive shrinkage (ashr) to model effect size distribution
  - Accounts for measurement uncertainty in effect estimates
  - Computes transcriptome-wide impact (variance of true effects)
  - Estimates effective number of differentially expressed genes (Me)
  - Supports gene set enrichment analysis with annotation tables

### Bivariate Analysis
- `trade_bivariate`: Estimate correlation of effects between two perturbations
  - Uses multivariate adaptive shrinkage (mashr) for joint modeling
  - Computes correlation between true (unobserved) effects
  - Accounts for sampling covariance when samples are shared
  - Supports flexible covariance matrix specifications
  - Returns both raw and corrected correlation estimates

## Input Data Format

Both tools expect DESeq2-style results files in CSV format with the following columns:
- `log2FoldChange`: Effect size estimates
- `lfcSE`: Standard errors of effect sizes
- `pvalue`: P-values for differential expression

Gene names should be in the row names or first column.

## Analysis Workflow

### Univariate Mode
1. Prepare DESeq2 results from a perturbation experiment
2. Call `trade_univariate` with the results file
3. Optionally exclude perturbed genes and provide annotation tables
4. Get transcriptome-wide impact and effective number of DEGs

### Bivariate Mode
1. Prepare DESeq2 results from two perturbation experiments
2. Ensure gene naming is consistent between files
3. Call `trade_bivariate` with both results files
4. Get correlation between perturbation effects across the transcriptome

## About TRADE

TRADE (Transcriptome-wide Analysis of Differential Expression) is a statistical framework for quantifying the global impact of genetic perturbations on gene expression. Key features include:

- Adaptive shrinkage modeling of effect size distributions
- Proper accounting of measurement uncertainty
- Estimation of transcriptome-wide impact (total variance of effects)
- Correlation estimation between multiple perturbations
- Based on ashr and mashr statistical frameworks

For more details, see the [TRADEtools documentation](https://github.com/SONGDONGYUAN1994/TRADEtools).
