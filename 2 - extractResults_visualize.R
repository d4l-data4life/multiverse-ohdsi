# =============================================================================
# Multiverse analysis with OHDSI tools: extraction and visualisation template
# -----------------------------------------------------------------------------
# Run after 1 - runMultiverseAnalysis.R has finished. Works on any Strategus
# results folder containing cm_result.csv.
#
# Steps
#   1  prepare dependencies
#   2  read and inspect the multiverse   
#   3  visualise
# =============================================================================


# =============================================================================
# 1  Prepare dependencies
# =============================================================================

library(ggplot2)
library(jsonlite)
library(patchwork)

source("FUNS_multiverseOHDSI.R")


# =============================================================================
# 2  Read and inspect the multiverse
# =============================================================================

# Point at the folder passed as `resultsFolder` to
# Strategus::createCdmExecutionSettings(). Searched recursively, so the parent
# output folder also works.

# resultsFolder <- file.path(getwd(), "resultsFolder", "results_folder")
resultsFolder <- file.path(getwd(), "runs", "multiverse_20260901_132313", "results_folder")

mv <- readMultiverse(resultsFolder)

nrow(mv$estimates)   # analyses returned
mv$estimates         # one row per analysis: effect estimate, CI, p, plus every argument joined on analysis_id
mv$grid              # the specification grid alone, full dotted argument paths
mv$varying           # arguments taking more than one value --> candidate rows for the specification dashboard
mv$derived           # columns describing which functions were used, not their arguments (psAdjustment, uses.*)
mv$gridSource        # which file the grid was read from

# Print every varying argument, its distinct values, and the distinct combinations realised.
spec <- inspectMultiverseSpec(mv, all = FALSE) # all = TRUE --> also lists arguments held constant

spec$values                             # one row per argument, with its distinct values
spec$combinations                       # distinct combinations, and how many analyses realise each


# =============================================================================
# 3  Visualise
# =============================================================================

# --- Specification curve -----------------------------------------------------
# Estimates ranked smallest to largest, over a dashboard of the decisions you
# name. Names become the row labels; row order follows the order given.

spec_curve <- createSpecificationCurve(
  mv,
  decisions = c(
    "PS adjustment"      = "psAdjustment",
    "Risk window (days)" = "createStudyPopulationArgs.riskWindowEnd",
    "Caliper"            = "matchOnPsArgs.caliper",
    "Match ratio"        = "matchOnPsArgs.maxRatio",
    "PS strata"          = "stratifyByPsArgs.numberOfStrata",
    "Trim fraction"      = "trimByPsArgs.trimFraction",
    "Max weight"         = "truncateIptwArgs.maxWeight"
  ),
  baseSize = 14
)
spec_curve

# --- Volcano plot ------------------------------------------------------------
# -log10(p) against the estimate, one point per specification. The shaded band
# is p > 0.05; dotted lines mark the 1st, 50th and 99th percentiles of the
# estimate. The spread between the outer two is the relative effect size ratio.
# colourBy takes one decision only.

createVolcanoPlot(mv, colourBy = c("PS adjustment" = "psAdjustment"))

volcano
