# Fix multiprocessing deadlock in asyncio context
# Must be called before any other imports
import multiprocessing
multiprocessing.set_start_method("spawn", force=True)

"""
Model Context Protocol (MCP) for TRADEtools

TRADEtools (Transcriptome-wide Analysis of Differential Expression) provides statistical
methods for analyzing differential expression data from RNA-seq experiments. It estimates
transcriptome-wide impact of perturbations by inferring effect size distributions using
adaptive shrinkage, and computes correlation between perturbation effects in bivariate mode.

This MCP Server provides Python interfaces to R tools extracted from the following tutorial files:
1. tradetools_intro
    - trade_univariate: Run univariate TRADE analysis to estimate transcriptome-wide impact (calls R via Rscript)
    - trade_bivariate: Run bivariate TRADE analysis to estimate correlation between perturbations (calls R via Rscript)

Note: All tools execute R code via Rscript subprocess calls. Ensure R is installed
and the TRADEtools package dependencies are available in the renv environment at repo/TRADEtools/.
"""

from mcp.server.fastmcp import FastMCP

# Import tool functions from modules
from tools.tradetools_intro import trade_univariate, trade_bivariate

# Server definition
mcp = FastMCP(name="TRADEtools")

# Register tools on the unified server
mcp.tool()(trade_univariate)
mcp.tool()(trade_bivariate)

if __name__ == "__main__":
    mcp.run()
