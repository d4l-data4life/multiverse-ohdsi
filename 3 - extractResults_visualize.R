source("FUNS_multiverseOHDSI.R")

mv <- readMultiverse("../resultsFolder/results_folder")

inspectMultiverseSpec(mv)

createSpecificationCurve(mv,
                         decisions = c("PS adjustment"  = "psAdjustment",
                                       "Washout (days)" = "createStudyPopArgs.washoutPeriod",
                                       "Caliper"        = "matchOnPsArgs.caliper",
                                       "Match ratio"    = "matchOnPsArgs.maxRatio",
                                       "PS strata"      = "stratifyByPsArgs.numberOfStrata",
                                       "Trim fraction"  = "trimByPsArgs.trimFraction",
                                       "Max weight"     = "truncateIptwArgs.maxWeight"))

createVolcanoPlot(mv, colourBy = c("PS adjustment" = "psAdjustment"))
