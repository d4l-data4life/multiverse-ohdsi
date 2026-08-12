# =============================================================================
# Multiverse Analysis with Strategus
# =============================================================================
# renv::snapshot()
 renv::restore()

# =============================================================================
# Define Database and Cohorts
# =============================================================================

cohortDefinitionSet <- CohortGenerator::getCohortDefinitionSet(
  settingsFileName = "testdata/Cohorts.csv",
  jsonFolder = "testdata/cohorts",
  sqlFolder = "testdata/sql",
  packageName = "Strategus"
)

cgModule <- Strategus::CohortGeneratorModule$new()

cohortDefinitionSharedResource <- cgModule$createCohortSharedResourceSpecifications(
  cohortDefinitionSet = cohortDefinitionSet
)

# =============================================================================
# Define Study (Target Estimand)
# =============================================================================

covs2exclude <- c(1118084, 1124300)

covSettings <- FeatureExtraction::createDefaultCovariateSettings(
  excludedCovariateConceptIds = covs2exclude,
  addDescendantsToExclude = TRUE
)

getDbCmDataArgs <- CohortMethod::createGetDbCohortMethodDataArgs(
  covariateSettings = covSettings
)

tco <- CohortMethod::createTargetComparatorOutcomes(
  targetId     = 2,
  comparatorId = 1,
  outcomes     = list(CohortMethod::createOutcome(outcomeId = 3))
)

studyPopArgs_wash180 <- CohortMethod::createCreateStudyPopulationArgs(
  washoutPeriod                  = 180,
  removeDuplicateSubjects        = "remove all",
  removeSubjectsWithPriorOutcome = TRUE,
  startAnchor                    = "cohort start",
  riskWindowStart                = 1,
  endAnchor                      = "cohort start",
  riskWindowEnd                  = 280
)

studyPopArgs_wash90 <- CohortMethod::createCreateStudyPopulationArgs(
  washoutPeriod                  = 90,
  removeDuplicateSubjects        = "remove all",
  removeSubjectsWithPriorOutcome = TRUE,
  startAnchor                    = "cohort start",
  riskWindowStart                = 1,
  endAnchor                      = "cohort start",
  riskWindowEnd                  = 280
)

studyPopArgs_wash0 <- CohortMethod::createCreateStudyPopulationArgs(
  washoutPeriod                  = 0,
  removeDuplicateSubjects        = "remove all",
  removeSubjectsWithPriorOutcome = TRUE,
  startAnchor                    = "cohort start",
  riskWindowStart                = 1,
  endAnchor                      = "cohort start",
  riskWindowEnd                  = 280
)

# =============================================================================
# Helper: Matched Cox Analysis
# =============================================================================

createMatchedCoxAnalysis <- function(analysisId, caliper, maxRatio, createStudyPopArgs,
                                     createPsArgs = CohortMethod::createCreatePsArgs(estimator = "att")) {
  description <- glue::glue("1:{maxRatio} matched Cox, caliper = {caliper}")

  CohortMethod::createCmAnalysis(
    analysisId                        = analysisId,
    description                       = description,
    getDbCohortMethodDataArgs         = getDbCmDataArgs,
    createStudyPopArgs                = createStudyPopArgs,
    createPsArgs                      = createPsArgs,
    matchOnPsArgs                     = CohortMethod::createMatchOnPsArgs(caliper = caliper, maxRatio = maxRatio),
    computeSharedCovariateBalanceArgs = CohortMethod::createComputeCovariateBalanceArgs(),
    computeCovariateBalanceArgs       = CohortMethod::createComputeCovariateBalanceArgs(
      covariateFilter = FeatureExtraction::getDefaultTable1Specifications()
    ),
    fitOutcomeModelArgs               = CohortMethod::createFitOutcomeModelArgs(modelType = "cox")
  )
}

# =============================================================================
# Matching Analyses (1-12)
# =============================================================================

cmAnalysis1  <- createMatchedCoxAnalysis(analysisId = 1,  caliper = 0.2,    maxRatio = 1,  studyPopArgs_wash180)
cmAnalysis2  <- createMatchedCoxAnalysis(analysisId = 2,  caliper = 0.2,    maxRatio = 10, studyPopArgs_wash180)
cmAnalysis3  <- createMatchedCoxAnalysis(analysisId = 3,  caliper = 0.0001, maxRatio = 1,  studyPopArgs_wash180)
cmAnalysis4  <- createMatchedCoxAnalysis(analysisId = 4,  caliper = 0.0001, maxRatio = 10, studyPopArgs_wash180)

cmAnalysis5  <- createMatchedCoxAnalysis(analysisId = 5,  caliper = 0.2,    maxRatio = 1,  studyPopArgs_wash90)
cmAnalysis6  <- createMatchedCoxAnalysis(analysisId = 6,  caliper = 0.2,    maxRatio = 10, studyPopArgs_wash90)
cmAnalysis7  <- createMatchedCoxAnalysis(analysisId = 7,  caliper = 0.0001, maxRatio = 1,  studyPopArgs_wash90)
cmAnalysis8  <- createMatchedCoxAnalysis(analysisId = 8,  caliper = 0.0001, maxRatio = 10, studyPopArgs_wash90)

cmAnalysis9   <- createMatchedCoxAnalysis(analysisId = 9,   caliper = 0.2,    maxRatio = 1,  studyPopArgs_wash0)
cmAnalysis10  <- createMatchedCoxAnalysis(analysisId = 10,  caliper = 0.2,    maxRatio = 10, studyPopArgs_wash0)
cmAnalysis11  <- createMatchedCoxAnalysis(analysisId = 11,  caliper = 0.0001, maxRatio = 1,  studyPopArgs_wash0)
cmAnalysis12  <- createMatchedCoxAnalysis(analysisId = 12,  caliper = 0.0001, maxRatio = 10, studyPopArgs_wash0)

# =============================================================================
# Helper: Stratified Cox Analysis
# =============================================================================

createStratifiedCoxAnalysis <- function(analysisId, numberOfStrata, createStudyPopArgs,
                                        createPsArgs = CohortMethod::createCreatePsArgs(estimator = "att")) {
  description <- glue::glue("Stratified Cox, strata = {numberOfStrata}")

  CohortMethod::createCmAnalysis(
    analysisId                        = analysisId,
    description                       = description,
    getDbCohortMethodDataArgs         = getDbCmDataArgs,
    createStudyPopArgs                = createStudyPopArgs,
    createPsArgs                      = createPsArgs,
    stratifyByPsArgs                  = CohortMethod::createStratifyByPsArgs(numberOfStrata = numberOfStrata),
    computeSharedCovariateBalanceArgs = CohortMethod::createComputeCovariateBalanceArgs(),
    computeCovariateBalanceArgs       = CohortMethod::createComputeCovariateBalanceArgs(
      covariateFilter = FeatureExtraction::getDefaultTable1Specifications()
    ),
    fitOutcomeModelArgs               = CohortMethod::createFitOutcomeModelArgs(modelType = "cox", stratified = TRUE)
  )
}

# =============================================================================
# Stratification Analyses (13-30)
# =============================================================================

cmAnalysis13  <- createStratifiedCoxAnalysis(analysisId = 13,  numberOfStrata = 4, studyPopArgs_wash180)
cmAnalysis14  <- createStratifiedCoxAnalysis(analysisId = 14,  numberOfStrata = 5, studyPopArgs_wash180)
cmAnalysis15  <- createStratifiedCoxAnalysis(analysisId = 15,  numberOfStrata = 6, studyPopArgs_wash180)
cmAnalysis16  <- createStratifiedCoxAnalysis(analysisId = 16,  numberOfStrata = 7, studyPopArgs_wash180)
cmAnalysis17  <- createStratifiedCoxAnalysis(analysisId = 17,  numberOfStrata = 8, studyPopArgs_wash180)
cmAnalysis18  <- createStratifiedCoxAnalysis(analysisId = 18,  numberOfStrata = 9, studyPopArgs_wash180)

cmAnalysis19  <- createStratifiedCoxAnalysis(analysisId = 19,  numberOfStrata = 4, studyPopArgs_wash90)
cmAnalysis20  <- createStratifiedCoxAnalysis(analysisId = 20,  numberOfStrata = 5, studyPopArgs_wash90)
cmAnalysis21  <- createStratifiedCoxAnalysis(analysisId = 21,  numberOfStrata = 6, studyPopArgs_wash90)
cmAnalysis22  <- createStratifiedCoxAnalysis(analysisId = 22,  numberOfStrata = 7, studyPopArgs_wash90)
cmAnalysis23  <- createStratifiedCoxAnalysis(analysisId = 23,  numberOfStrata = 8, studyPopArgs_wash90)
cmAnalysis24  <- createStratifiedCoxAnalysis(analysisId = 24,  numberOfStrata = 9, studyPopArgs_wash90)

cmAnalysis25  <- createStratifiedCoxAnalysis(analysisId = 25,  numberOfStrata = 4, studyPopArgs_wash0)
cmAnalysis26  <- createStratifiedCoxAnalysis(analysisId = 26,  numberOfStrata = 5, studyPopArgs_wash0)
cmAnalysis27  <- createStratifiedCoxAnalysis(analysisId = 27,  numberOfStrata = 6, studyPopArgs_wash0)
cmAnalysis28  <- createStratifiedCoxAnalysis(analysisId = 28,  numberOfStrata = 7, studyPopArgs_wash0)
cmAnalysis29  <- createStratifiedCoxAnalysis(analysisId = 29,  numberOfStrata = 8, studyPopArgs_wash0)
cmAnalysis30  <- createStratifiedCoxAnalysis(analysisId = 30,  numberOfStrata = 9, studyPopArgs_wash0)

# =============================================================================
# Helper: Weighted (IPTW) Cox Analysis
# =============================================================================

createWeightedCoxAnalysis <- function(analysisId, estimator,
                                      maxWeight = NULL,
                                      trimPercentile = 0,
                                      createStudyPopArgs,
                                      createPsArgs) {
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
    createStudyPopArgs                = createStudyPopArgs,
    createPsArgs                      = CohortMethod::createCreatePsArgs(estimator = estimator),
    trimByPsArgs                      = trimByPsArgs,
    truncateIptwArgs                  = truncateIptwArgs,
    computeSharedCovariateBalanceArgs = CohortMethod::createComputeCovariateBalanceArgs(),
    computeCovariateBalanceArgs       = CohortMethod::createComputeCovariateBalanceArgs(
      covariateFilter = FeatureExtraction::getDefaultTable1Specifications()
    ),
    fitOutcomeModelArgs               = CohortMethod::createFitOutcomeModelArgs(
      modelType        = "cox",
      inversePtWeighting = TRUE
    )
  )
}

# =============================================================================
# Weighting Analyses (31-22)
# =============================================================================

# ATT
cmAnalysis31 <- createWeightedCoxAnalysis(31, estimator = "att", trimPercentile = 0, maxWeight = NULL, createStudyPopArgs = studyPopArgs_wash180)  # no trim, no truncation
cmAnalysis32 <- createWeightedCoxAnalysis(32, estimator = "att", trimPercentile = 1, maxWeight = NULL, createStudyPopArgs = studyPopArgs_wash180)  # trim only
cmAnalysis33 <- createWeightedCoxAnalysis(33, estimator = "att", trimPercentile = 0, maxWeight = 10,   createStudyPopArgs = studyPopArgs_wash180)  # truncate only
cmAnalysis34 <- createWeightedCoxAnalysis(34, estimator = "att", trimPercentile = 1, maxWeight = 10,   createStudyPopArgs = studyPopArgs_wash180)  # both

cmAnalysis35 <- createWeightedCoxAnalysis(35, estimator = "att", trimPercentile = 0, maxWeight = NULL, createStudyPopArgs = studyPopArgs_wash90)  # no trim, no truncation
cmAnalysis36 <- createWeightedCoxAnalysis(36, estimator = "att", trimPercentile = 1, maxWeight = NULL, createStudyPopArgs = studyPopArgs_wash90)  # trim only
cmAnalysis37 <- createWeightedCoxAnalysis(37, estimator = "att", trimPercentile = 0, maxWeight = 10,   createStudyPopArgs = studyPopArgs_wash90)  # truncate only
cmAnalysis38 <- createWeightedCoxAnalysis(38, estimator = "att", trimPercentile = 1, maxWeight = 10,   createStudyPopArgs = studyPopArgs_wash90)  # both

cmAnalysis39 <- createWeightedCoxAnalysis(39, estimator = "att", trimPercentile = 0, maxWeight = NULL, createStudyPopArgs = studyPopArgs_wash0)  # no trim, no truncation
cmAnalysis40 <- createWeightedCoxAnalysis(40, estimator = "att", trimPercentile = 1, maxWeight = NULL, createStudyPopArgs = studyPopArgs_wash0)  # trim only
cmAnalysis41 <- createWeightedCoxAnalysis(41, estimator = "att", trimPercentile = 0, maxWeight = 10,   createStudyPopArgs = studyPopArgs_wash0)  # truncate only
cmAnalysis42 <- createWeightedCoxAnalysis(42, estimator = "att", trimPercentile = 1, maxWeight = 10,   createStudyPopArgs = studyPopArgs_wash0)  # both


# # ATE
# cmAnalysis15 <- createWeightedCoxAnalysis(15, estimator = "ate", trimPercentile = 0, maxWeight = NULL, createStudyPopArgs = studyPopArgs)  # no trim, no truncation
# cmAnalysis16 <- createWeightedCoxAnalysis(16, estimator = "ate", trimPercentile = 1, maxWeight = NULL, createStudyPopArgs = studyPopArgs)  # trim only
# cmAnalysis17 <- createWeightedCoxAnalysis(17, estimator = "ate", trimPercentile = 0, maxWeight = 10,   createStudyPopArgs = studyPopArgs)  # truncate only
# cmAnalysis18 <- createWeightedCoxAnalysis(18, estimator = "ate", trimPercentile = 1, maxWeight = 10,   createStudyPopArgs = studyPopArgs)  # both
# 
# # ATO
# cmAnalysis19 <- createWeightedCoxAnalysis(19, estimator = "ato", trimPercentile = 0, maxWeight = NULL, createStudyPopArgs = studyPopArgs)  # no trim, no truncation
# cmAnalysis20 <- createWeightedCoxAnalysis(20, estimator = "ato", trimPercentile = 1, maxWeight = NULL, createStudyPopArgs = studyPopArgs)  # trim only
# cmAnalysis21 <- createWeightedCoxAnalysis(21, estimator = "ato", trimPercentile = 0, maxWeight = 10,   createStudyPopArgs = studyPopArgs)  # truncate only
# cmAnalysis22 <- createWeightedCoxAnalysis(22, estimator = "ato", trimPercentile = 1, maxWeight = 10,   createStudyPopArgs = studyPopArgs)  # both

# =============================================================================
# Setup CohortMethod Module
# =============================================================================

cmModule <- Strategus::CohortMethodModule$new()

cohortMethodModuleSpecifications <- cmModule$createModuleSpecifications(
  cmAnalysisList               = mget(paste0("cmAnalysis", 1:42)),
  targetComparatorOutcomesList = list(tco)
)

# =============================================================================
# Setup CohortGenerator Module
# =============================================================================

cgModule <- Strategus::CohortGeneratorModule$new()

cohortGeneratorModuleSpecifications <- cgModule$createModuleSpecifications(
  generateStats = TRUE
)

# =============================================================================
# Define Analysis Specification
# =============================================================================

analysisSpecifications <- Strategus::createEmptyAnalysisSpecificiations() |>
  Strategus::addSharedResources(cohortDefinitionSharedResource) |>
  Strategus::addModuleSpecifications(cohortGeneratorModuleSpecifications) |>
  Strategus::addModuleSpecifications(cohortMethodModuleSpecifications)

# =============================================================================
# Connection Details
# =============================================================================

connectionDetails <- Eunomia::getEunomiaConnectionDetails()
# create cohorts
Eunomia::createCohorts(connectionDetails,
                       cdmDatabaseSchema = "main",
                       cohortDatabaseSchema = "main", 
                       cohortTable = "cohort")


# =============================================================================
# Execution Settings
# =============================================================================
# setwd("C:/Users/MarioLawes/Desktop")
outputFolder <- file.path(getwd(), "resultsFolder")
# Delete existing folder and recreate fresh
unlink(outputFolder, recursive = TRUE)
dir.create(outputFolder, showWarnings = FALSE)

executionSettings <- Strategus::createCdmExecutionSettings(
  workDatabaseSchema = "main",
  cdmDatabaseSchema  = "main",
  cohortTableNames   = CohortGenerator::getCohortTableNames(),
  workFolder         = file.path(outputFolder, "work_folder"),
  resultsFolder      = file.path(outputFolder, "results_folder"),
  minCellCount       = 5
)

ParallelLogger::saveSettingsToJson(
  analysisSpecifications,
  file.path(outputFolder, "analysisSpecification.json")
)

# =============================================================================
# Execute
# =============================================================================

Strategus::execute(
  analysisSpecifications,
  executionSettings,
  connectionDetails
)

