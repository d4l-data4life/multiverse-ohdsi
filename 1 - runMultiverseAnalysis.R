# =============================================================================
# Multiverse analysis template using the OHDSI R package Strategus
# -----------------------------------------------------------------------------
# Purpose:   Illustrative multiverse analysis example using the Eunomia GiBleed
#            demo CDM.
#
# Estimand:  Effect of new use of CELECOXIB (target, cohort 1) relative to new
#            use of DICLOFENAC (comparator, cohort 2) on incident GI BLEED
#            (outcome, cohort 3), fixed risk window, ATT scale throughout.
#
# Multiverse declaration (42 specifications):
#     matching        caliper {0.2, 0.0001} x maxRatio {1, 10}        ->  4
#     stratification  numberOfStrata {4, 5, 6, 7, 8, 9}               ->  6
#     IPTW            trim {0, 1%} x truncation {none, maxWeight 10}  ->  4
#   Risk window end   {30, 60, 90} days                               ->  x3
#
# HELD FIXED (each a researcher degree of freedom NOT explored here):
#   washoutPeriod = 180, riskWindowStart = 1, both anchors at cohort start
#   removeDuplicateSubjects = "remove all"
#   removeSubjectsWithPriorOutcome = TRUE
#   covariates = FeatureExtraction defaults, excluding 1118084 (diclofenac)
#     and 1124300 (celecoxib) plus descendants
#   PS estimator = "att" for all 42 specifications
#   outcome model = Cox
#   No positive/negative control outcomes, no empirical calibration
# -----------------------------------------------------------------------------
#
# Output:  runs/multiverse_<YYYYMMDD_HHMMSS>/
#            analysisSpecification.json   (machine-readable declaration of the
#              full multiverse - the pre-registerable artifact)
#            results_folder/CohortMethodModule/*.csv
#              (consumed by 2_-_extractResults_visualize.R)
# =============================================================================


# =============================================================================
# Environment
# =============================================================================
# Restores the exact package versions recorded in renv.lock. Run once per
# machine; re-run after any renv::snapshot(). HADES packages are version-
# sensitive, so results are only reproducible against this lockfile.

# download.file(
#  "https://raw.githubusercontent.com/ohdsi-studies/StrategusStudyRepoTemplate/main/renv.lock",
#  destfile = "renv.lock")
# renv::snapshot()
# renv::status()
# renv::restore()

# =============================================================================
# Output Folder
# =============================================================================
# Each execution gets its own timestamped folder, so runs accumulate rather
# than overwrite. Defined up front so the multiverse declaration can be
# written before execution begins.

runId        <- format(Sys.time(), "%Y%m%d_%H%M%S")
outputFolder <- file.path(getwd(), "runs", paste0("multiverse_", runId))
dir.create(outputFolder, recursive = TRUE)

# =============================================================================
# Define Database and Cohorts
# =============================================================================
# Loads the target, comparator, and outcome cohort definitions shipped with
# Strategus (Eunomia GiBleed demo) and registers them as a shared resource,
# so every one of the 42 specifications is built on identical cohorts.

cohortDefinitionSet <- CohortGenerator::getCohortDefinitionSet(
  settingsFileName = "testdata/Cohorts.csv",
  jsonFolder       = "testdata/cohorts",
  sqlFolder        = "testdata/sql",
  packageName      = "Strategus"
)

cgModule <- Strategus::CohortGeneratorModule$new()

cohortDefinitionSharedResource <- cgModule$createCohortSharedResourceSpecifications(
  cohortDefinitionSet = cohortDefinitionSet
)

# =============================================================================
# Define Study
# =============================================================================
# Everything the multiverse holds constant: covariate construction, the
# target/comparator/outcome triplet, and the three study populations that
# differ only in risk window length.

covs2exclude <- c(1118084, 1124300)

covSettings <- FeatureExtraction::createDefaultCovariateSettings(
  excludedCovariateConceptIds = covs2exclude,
  addDescendantsToExclude     = TRUE
)


tco <- CohortMethod::createTargetComparatorOutcomes(
  targetId     = 1,   # celecoxib
  comparatorId = 2,   # diclofenac
  outcomes     = list(CohortMethod::createOutcome(outcomeId = 3))
)

getDbCmDataArgs <- CohortMethod::createGetDbCohortMethodDataArgs(
  covariateSettings       = covSettings,
  washoutPeriod           = 180,
  removeDuplicateSubjects = "remove all",
  firstExposureOnly       = TRUE,
  restrictToCommonPeriod  = TRUE
)

makeStudyPopArgs <- function(riskWindowEnd) {
  CohortMethod::createCreateStudyPopulationArgs(
    removeSubjectsWithPriorOutcome = TRUE,
    startAnchor                    = "cohort start",
    riskWindowStart                = 1,
    endAnchor                      = "cohort start",
    riskWindowEnd                  = riskWindowEnd
  )
}

studyPop_30 <- makeStudyPopArgs(30)
studyPop_60 <- makeStudyPopArgs(60)
studyPop_90 <- makeStudyPopArgs(90)

# =============================================================================
# Helper: Matched Cox Analysis
# =============================================================================
# Factory for one PS-matched specification. Varies caliper (0.2 = Austin's
# recommendation, 0.0001 ~ exact matching) and maxRatio (1 = pair matching,
# 10 = variable ratio). Estimator is ATT by default.

createMatchedCoxAnalysis <- function(analysisId, caliper, maxRatio, studyPopArgs,
                                     createPsArgs = CohortMethod::createCreatePsArgs(estimator = "att")) {
  description <- glue::glue("1:{maxRatio} matched Cox, caliper = {caliper}")
  
  CohortMethod::createCmAnalysis(
    analysisId                        = analysisId,
    description                       = description,
    getDbCohortMethodDataArgs         = getDbCmDataArgs,
    createStudyPopulationArgs         = studyPopArgs,
    createPsArgs                      = createPsArgs,
    matchOnPsArgs                     = CohortMethod::createMatchOnPsArgs(caliper = caliper, maxRatio = maxRatio),
    #computeSharedCovariateBalanceArgs = CohortMethod::createComputeCovariateBalanceArgs(),
    #computeCovariateBalanceArgs       = CohortMethod::createComputeCovariateBalanceArgs(
    #  covariateFilter = FeatureExtraction::getDefaultTable1Specifications()
    #),
    fitOutcomeModelArgs               = CohortMethod::createFitOutcomeModelArgs(modelType = "cox")
  )
}

# =============================================================================
# Matching Analyses (1-12)
# =============================================================================
# The matching specs: 2 calipers x 2 ratios x 3 risk windows. = 12.

cmAnalysis1  <- createMatchedCoxAnalysis(analysisId = 1,  caliper = 0.2,    maxRatio = 1,  studyPopArgs = studyPop_30)
cmAnalysis2  <- createMatchedCoxAnalysis(analysisId = 2,  caliper = 0.2,    maxRatio = 10, studyPopArgs = studyPop_30)
cmAnalysis3  <- createMatchedCoxAnalysis(analysisId = 3,  caliper = 0.0001, maxRatio = 1,  studyPopArgs = studyPop_30)
cmAnalysis4  <- createMatchedCoxAnalysis(analysisId = 4,  caliper = 0.0001, maxRatio = 10, studyPopArgs = studyPop_30)

cmAnalysis5  <- createMatchedCoxAnalysis(analysisId = 5,  caliper = 0.2,    maxRatio = 1,  studyPopArgs = studyPop_60)
cmAnalysis6  <- createMatchedCoxAnalysis(analysisId = 6,  caliper = 0.2,    maxRatio = 10, studyPopArgs = studyPop_60)
cmAnalysis7  <- createMatchedCoxAnalysis(analysisId = 7,  caliper = 0.0001, maxRatio = 1,  studyPopArgs = studyPop_60)
cmAnalysis8  <- createMatchedCoxAnalysis(analysisId = 8,  caliper = 0.0001, maxRatio = 10, studyPopArgs = studyPop_60)

cmAnalysis9  <- createMatchedCoxAnalysis(analysisId = 9,  caliper = 0.2,    maxRatio = 1,  studyPopArgs = studyPop_90)
cmAnalysis10 <- createMatchedCoxAnalysis(analysisId = 10, caliper = 0.2,    maxRatio = 10, studyPopArgs = studyPop_90)
cmAnalysis11 <- createMatchedCoxAnalysis(analysisId = 11, caliper = 0.0001, maxRatio = 1,  studyPopArgs = studyPop_90)
cmAnalysis12 <- createMatchedCoxAnalysis(analysisId = 12, caliper = 0.0001, maxRatio = 10, studyPopArgs = studyPop_90)

# =============================================================================
# Helper: Stratified Cox Analysis
# =============================================================================
# Factory for one PS-stratified specification. Varies only the number of
# strata; the outcome model is stratified to match.

createStratifiedCoxAnalysis <- function(analysisId, numberOfStrata, studyPopArgs,
                                        createPsArgs = CohortMethod::createCreatePsArgs(estimator = "att")) {
  description <- glue::glue("Stratified Cox, strata = {numberOfStrata}")
  
  CohortMethod::createCmAnalysis(
    analysisId                        = analysisId,
    description                       = description,
    getDbCohortMethodDataArgs         = getDbCmDataArgs,
    createStudyPopulationArgs         = studyPopArgs,
    createPsArgs                      = createPsArgs,
    stratifyByPsArgs                  = CohortMethod::createStratifyByPsArgs(numberOfStrata = numberOfStrata),
    #computeSharedCovariateBalanceArgs = CohortMethod::createComputeCovariateBalanceArgs(),
    #computeCovariateBalanceArgs       = CohortMethod::createComputeCovariateBalanceArgs(
    #  covariateFilter = FeatureExtraction::getDefaultTable1Specifications()
    #),
    fitOutcomeModelArgs               = CohortMethod::createFitOutcomeModelArgs(modelType = "cox", stratified = TRUE)
  )
}

# =============================================================================
# Stratification Analyses (13-30)
# =============================================================================
# The stratification specs: 6 strata counts (4-9) x 3 risk windows. = 18.

cmAnalysis13 <- createStratifiedCoxAnalysis(analysisId = 13, numberOfStrata = 4, studyPopArgs = studyPop_30)
cmAnalysis14 <- createStratifiedCoxAnalysis(analysisId = 14, numberOfStrata = 5, studyPopArgs = studyPop_30)
cmAnalysis15 <- createStratifiedCoxAnalysis(analysisId = 15, numberOfStrata = 6, studyPopArgs = studyPop_30)
cmAnalysis16 <- createStratifiedCoxAnalysis(analysisId = 16, numberOfStrata = 7, studyPopArgs = studyPop_30)
cmAnalysis17 <- createStratifiedCoxAnalysis(analysisId = 17, numberOfStrata = 8, studyPopArgs = studyPop_30)
cmAnalysis18 <- createStratifiedCoxAnalysis(analysisId = 18, numberOfStrata = 9, studyPopArgs = studyPop_30)

cmAnalysis19 <- createStratifiedCoxAnalysis(analysisId = 19, numberOfStrata = 4, studyPopArgs = studyPop_60)
cmAnalysis20 <- createStratifiedCoxAnalysis(analysisId = 20, numberOfStrata = 5, studyPopArgs = studyPop_60)
cmAnalysis21 <- createStratifiedCoxAnalysis(analysisId = 21, numberOfStrata = 6, studyPopArgs = studyPop_60)
cmAnalysis22 <- createStratifiedCoxAnalysis(analysisId = 22, numberOfStrata = 7, studyPopArgs = studyPop_60)
cmAnalysis23 <- createStratifiedCoxAnalysis(analysisId = 23, numberOfStrata = 8, studyPopArgs = studyPop_60)
cmAnalysis24 <- createStratifiedCoxAnalysis(analysisId = 24, numberOfStrata = 9, studyPopArgs = studyPop_60)

cmAnalysis25 <- createStratifiedCoxAnalysis(analysisId = 25, numberOfStrata = 4, studyPopArgs = studyPop_90)
cmAnalysis26 <- createStratifiedCoxAnalysis(analysisId = 26, numberOfStrata = 5, studyPopArgs = studyPop_90)
cmAnalysis27 <- createStratifiedCoxAnalysis(analysisId = 27, numberOfStrata = 6, studyPopArgs = studyPop_90)
cmAnalysis28 <- createStratifiedCoxAnalysis(analysisId = 28, numberOfStrata = 7, studyPopArgs = studyPop_90)
cmAnalysis29 <- createStratifiedCoxAnalysis(analysisId = 29, numberOfStrata = 8, studyPopArgs = studyPop_90)
cmAnalysis30 <- createStratifiedCoxAnalysis(analysisId = 30, numberOfStrata = 9, studyPopArgs = studyPop_90)

# =============================================================================
# Helper: Weighted (IPTW) Cox Analysis
# =============================================================================
# Factory for one IPTW specification. Trimming and truncation are optional:
# passing trimPercentile = 0 or maxWeight = NULL omits the corresponding
# step from the cmAnalysis object entirely rather than passing a null-effect
# argument.

createWeightedCoxAnalysis <- function(analysisId, estimator,
                                      maxWeight      = NULL,
                                      trimPercentile = 0,
                                      studyPopArgs) {
  description <- glue::glue(
    "IPTW Cox, estimator = {estimator}, trim = {trimPercentile}%, max weight = {ifelse(is.null(maxWeight), 'none', maxWeight)}"
  )
  
  trimByPsArgs <- if (trimPercentile > 0) {
    CohortMethod::createTrimByPsArgs(trimFraction = trimPercentile / 100)
  } else {
    NULL
  }
  
  truncateIptwArgs <- if (!is.null(maxWeight)) {
    CohortMethod::createTruncateIptwArgs(maxWeight = maxWeight)
  } else {
    NULL
  }
  
  CohortMethod::createCmAnalysis(
    analysisId                        = analysisId,
    description                       = description,
    getDbCohortMethodDataArgs         = getDbCmDataArgs,
    createStudyPopulationArgs         = studyPopArgs,
    createPsArgs                      = CohortMethod::createCreatePsArgs(estimator = estimator),
    trimByPsArgs                      = trimByPsArgs,
    truncateIptwArgs                  = truncateIptwArgs,
    #computeSharedCovariateBalanceArgs = CohortMethod::createComputeCovariateBalanceArgs(),
    #computeCovariateBalanceArgs       = CohortMethod::createComputeCovariateBalanceArgs(
    #  covariateFilter = FeatureExtraction::getDefaultTable1Specifications()
    #),
    fitOutcomeModelArgs               = CohortMethod::createFitOutcomeModelArgs(
      modelType          = "cox",
      inversePtWeighting = TRUE
    )
  )
}

# =============================================================================
# Weighting Analyses (31-42)
# =============================================================================
# The IPTW specs: {trim, no trim} x {truncate, no truncate} x 3 risk windows
# periods = 12.

cmAnalysis31 <- createWeightedCoxAnalysis(31, estimator = "att", trimPercentile = 0, maxWeight = NULL, studyPopArgs = studyPop_30)  # no trim, no truncation
cmAnalysis32 <- createWeightedCoxAnalysis(32, estimator = "att", trimPercentile = 1, maxWeight = NULL, studyPopArgs = studyPop_30)  # trim only
cmAnalysis33 <- createWeightedCoxAnalysis(33, estimator = "att", trimPercentile = 0, maxWeight = 10,   studyPopArgs = studyPop_30)  # truncate only
cmAnalysis34 <- createWeightedCoxAnalysis(34, estimator = "att", trimPercentile = 1, maxWeight = 10,   studyPopArgs = studyPop_30)  # both

cmAnalysis35 <- createWeightedCoxAnalysis(35, estimator = "att", trimPercentile = 0, maxWeight = NULL, studyPopArgs = studyPop_60)   # no trim, no truncation
cmAnalysis36 <- createWeightedCoxAnalysis(36, estimator = "att", trimPercentile = 1, maxWeight = NULL, studyPopArgs = studyPop_60)   # trim only
cmAnalysis37 <- createWeightedCoxAnalysis(37, estimator = "att", trimPercentile = 0, maxWeight = 10,   studyPopArgs = studyPop_60)   # truncate only
cmAnalysis38 <- createWeightedCoxAnalysis(38, estimator = "att", trimPercentile = 1, maxWeight = 10,   studyPopArgs = studyPop_60)   # both

cmAnalysis39 <- createWeightedCoxAnalysis(39, estimator = "att", trimPercentile = 0, maxWeight = NULL, studyPopArgs = studyPop_90)    # no trim, no truncation
cmAnalysis40 <- createWeightedCoxAnalysis(40, estimator = "att", trimPercentile = 1, maxWeight = NULL, studyPopArgs = studyPop_90)    # trim only
cmAnalysis41 <- createWeightedCoxAnalysis(41, estimator = "att", trimPercentile = 0, maxWeight = 10,   studyPopArgs = studyPop_90)    # truncate only
cmAnalysis42 <- createWeightedCoxAnalysis(42, estimator = "att", trimPercentile = 1, maxWeight = 10,   studyPopArgs = studyPop_90)    # both

# =============================================================================
# Setup CohortMethod Module
# =============================================================================
# Collects all 42 cmAnalysis objects into a single module specification.
# mget() assumes the objects are named cmAnalysis1 - cmAnalysis42 with no gaps.

cmModule <- Strategus::CohortMethodModule$new()

cmAnalysesSpec <- CohortMethod::createCmAnalysesSpecifications(
  cmAnalysisList               = mget(paste0("cmAnalysis", 1:42)),
  targetComparatorOutcomesList = list(tco)
)

cohortMethodModuleSpecifications <- cmModule$createModuleSpecifications(
  cmAnalysesSpecifications = cmAnalysesSpec$toList()
)

# =============================================================================
# Setup CohortGenerator Module
# =============================================================================
# Instructs Strategus to instantiate the cohorts and emit inclusion-rule
# attrition statistics before the estimation module runs. Reuses the cgModule
# instantiated above under "Define Database and Cohorts".

cohortGeneratorModuleSpecifications <- cgModule$createModuleSpecifications(
  generateStats = TRUE
)

# =============================================================================
# Define Analysis Specification
# =============================================================================
# Assembles shared resources and both modules into the complete, executable
# declaration of the multiverse. This object *is* the multiverse; everything
# after this point is execution.
#
# saveSettingsToJson() writes that declaration to disk before execution. The
# JSON is the shareable, pre-registerable artifact: it fully determines the 42
# specifications independently of this script and can be re-executed via
# Strategus::loadAnalysisSpecifications(). Writing it first means a failed run
# still leaves the declaration behind.

analysisSpecifications <- Strategus::createEmptyAnalysisSpecifications() |>
  Strategus::addSharedResources(cohortDefinitionSharedResource) |>
  Strategus::addModuleSpecifications(cohortGeneratorModuleSpecifications) |>
  Strategus::addModuleSpecifications(cohortMethodModuleSpecifications)

# save the spec as a machine-readable JSON file
ParallelLogger::saveSettingsToJson(
  analysisSpecifications,
  file.path(outputFolder, "analysisSpecification.json")
)

# =============================================================================
# Connection Details
# =============================================================================
# Points at the Eunomia demo CDM (a temporary local SQLite copy) and creates
# its cohort tables. Replace this block alone to run against a real OMOP CDM;
# nothing above it needs to change.

connectionDetails <- Eunomia::getEunomiaConnectionDetails()

# Eunomia::createCohorts(
#   connectionDetails,
#   cdmDatabaseSchema    = "main",
#   cohortDatabaseSchema = "main",
#   cohortTable          = "cohort"
# )

# =============================================================================
# Execution Settings
# =============================================================================
# Working and results directories inside the timestamped run folder, plus
# minimum cell count for disclosure control.

executionSettings <- Strategus::createCdmExecutionSettings(
  workDatabaseSchema = "main",
  cdmDatabaseSchema  = "main",
  cohortTableNames   = CohortGenerator::getCohortTableNames(),
  workFolder         = file.path(outputFolder, "work_folder"),
  resultsFolder      = file.path(outputFolder, "results_folder"),
  minCellCount       = 5
)

# =============================================================================
# Execute
# =============================================================================
# Runs all 42 specifications sequentially and writes per-analysis results to
# <run folder>/results_folder/CohortMethodModule/.

Strategus::execute(
  analysisSpecifications,
  executionSettings,
  connectionDetails
)
