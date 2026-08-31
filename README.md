# Multiverse analysis for OHDSI

This repository contains template code for conducting a multiverse analysis with
[Strategus](https://github.com/OHDSI/Strategus) and the
[HADES](https://ohdsi.github.io/Hades/) toolchain.

## Motivation

The OHDSI community has built substantial infrastructure for reproducible observational
research. The OMOP Common Data Model standardises the structure and semantics of health
data so that one analysis can run across heterogeneous sources (Hripcsak et al., 2015),
and an open-source toolchain supplies shared implementations of cohort definition,
covariate construction, and population-level effect estimation (HADES). Methodological
best practice accompanies these tools: new-user active-comparator designs, large-scale
propensity score adjustment (Suchard et al., 2013), negative control outcomes with
empirical calibration (Schuemie et al., 2018), and pre-specified diagnostics that gate
whether an estimate may be unblinded (Conover et al., 2025).

Harmonising the data and the tooling removes variability in implementation of real world evidence studies, and best
practice settles many design questions. However, it does not settle all of them. How long a washout period should be, whether to match,
stratify, or weight on the propensity score, how wide a caliper to set, how many strata
to use, whether to trim or truncate extreme weights — guidance favours no single option,
yet each can move the estimate.

Multiverse analysis addresses precisely this (Steegen et al., 2016). Instead of reporting
one pipeline, the analyst defines all defensible alternatives, executes them,
and reports the distribution of findings. Related approaches include specification curve
analysis (Simonsohn et al., 2020) and, from epidemiology, vibration of effects
(Patel et al., 2015; Vinatier et al., 2025).

Building on OHDSI methodology, we propose a structured framework for multiverse analysis in large-scale RWE and 
give guidance on how to interpret the resulting distribution of estimates.
## Illustrative example

For demonstration purposes the template runs against the synthetic
[Eunomia](https://github.com/OHDSI/Eunomia) GiBleed CDM.

| | |
|---|---|
| **Target** | Diclofenac (cohort 2) |
| **Comparator** | Celecoxib (cohort 1) |
| **Outcome** | Gastrointestinal bleed (cohort 3) |
| **Data** | Eunomia GiBleed synthetic CDM |
| **Estimand** | ATT |
| **Outcome model** | Cox proportional hazards |
| **Risk window** | Days 1–280, anchored on cohort start |

Covariates for the two exposure ingredients, and their descendants, are excluded from the
propensity model.

### Specification space

The 42 specifications are the crossing of three adjustment families with three washout
periods:

| Family | Varied arguments | Levels | n |
|---|---|---|---|
| PS matching | `caliper` ∈ {0.2, 0.0001} × `maxRatio` ∈ {1, 10} | 4 | 4 |
| PS stratification | `numberOfStrata` ∈ {4, 5, 6, 7, 8, 9} | 6 | 6 |
| IPTW | `trimFraction` ∈ {0, 0.01} × `maxWeight` ∈ {none, 10} | 4 | 4 |
| | **Washout period** ∈ {0, 90, 180} days | 3 | ×3 |
| | | | **42** |

Analysis IDs map as follows: 1–12 matching, 13–30 stratification, 31–42 weighting.
Commented-out blocks at the end of `1_-_runMultiverseAnalysis.R` sketch ATE and ATO
extensions that are not part of the current specification space.

## Repository contents

```
1_-_runMultiverseAnalysis.R    Build the 42-specification Strategus analysis
                               specification and execute it against Eunomia
2_-_extractResults_visualize.R Parse the executed specifications back out of
                               cm_analysis.csv, join to estimates, and draw
                               the specification curve
```

Execution writes to `resultsFolder/`, which is deleted and recreated on each run:

```
resultsFolder/
├── work_folder/
└── results_folder/
    └── CohortMethodModule/
        ├── cm_analysis.csv   One row per specification; `definition` holds
        │                     the analysis as JSON
        └── cm_result.csv     Effect estimates, CIs, p-values per analysis
```

## Requirements

`renv` is used for dependency pinning; restore the recorded library with:

```r
renv::restore()
```

The cohort definitions are read from `testdata/` **inside the installed Strategus
package** (`getCohortDefinitionSet(..., packageName = "Strategus")`), not from this
repository, so no cohort JSON is vendored here.

## Running the analysis

```r
source("1_-_runMultiverseAnalysis.R")    # builds and executes all 42 specifications
source("2_-_extractResults_visualize.R") # extracts results and plots
```

The first script downloads and instantiates the Eunomia CDM, creates the demo cohorts,
and calls `Strategus::execute()`. Runtime is a few minutes on the demo data.

The second script recovers the varied arguments by parsing the `definition` JSON column
of `cm_analysis.csv`. The specifications are read back out of the executed artefacts
rather than re-declared, so the plot cannot drift from what was actually run. Output is a
two-panel specification curve: ranked estimates with 95% CIs on top, and the corresponding
specification grid below, faceted by argument family.

## Citation

Lawes M, Sampri A, Illigens B. *Multiverse analysis for large-scale real-world evidence
using OHDSI tools.* Manuscript in preparation.

## References

Conover MM, Ryan PB, Chen Y, Suchard MA, Hripcsak G, Schuemie MJ. Objective study
validity diagnostics: a framework requiring pre-specified, empirical verification to
increase trust in the reliability of real-world evidence. *Journal of the American
Medical Informatics Association*. 2025;32(3):518–525.
doi:[10.1093/jamia/ocae317](https://doi.org/10.1093/jamia/ocae317)

Hripcsak G, Duke JD, Shah NH, et al. Observational Health Data Sciences and Informatics
(OHDSI): opportunities for observational researchers. *Studies in Health Technology and
Informatics*. 2015;216:574–578.

Observational Health Data Sciences and Informatics. *HADES: Health Analytics
Data-to-Evidence Suite.* <https://ohdsi.github.io/Hades/>

Patel CJ, Burford B, Ioannidis JPA. Assessment of vibration of effects due to model
specification can demonstrate the instability of observational associations. *Journal of
Clinical Epidemiology*. 2015;68(9):1046–1058.
doi:[10.1016/j.jclinepi.2015.05.029](https://doi.org/10.1016/j.jclinepi.2015.05.029)

Schuemie MJ, Hripcsak G, Ryan PB, Madigan D, Suchard MA. Empirical confidence interval
calibration for population-level effect estimation studies in observational healthcare
data. *Proceedings of the National Academy of Sciences*. 2018;115(11):2571–2577.
doi:[10.1073/pnas.1708282114](https://doi.org/10.1073/pnas.1708282114)

Simonsohn U, Simmons JP, Nelson LD. Specification curve analysis. *Nature Human
Behaviour*. 2020;4(11):1208–1214.
doi:[10.1038/s41562-020-0912-z](https://doi.org/10.1038/s41562-020-0912-z)

Steegen S, Tuerlinckx F, Gelman A, Vanpaemel W. Increasing transparency through a
multiverse analysis. *Perspectives on Psychological Science*. 2016;11(5):702–712.
doi:[10.1177/1745691616658637](https://doi.org/10.1177/1745691616658637)

Suchard MA, Simpson SE, Zorych I, Ryan P, Madigan D. Massive parallelization of serial
inference algorithms for a complex generalized linear model. *ACM Transactions on
Modeling and Computer Simulation*. 2013;23(1):10.
doi:[10.1145/2414416.2414791](https://doi.org/10.1145/2414416.2414791)

Vinatier C, Hoffmann S, Patel C, et al. What is the vibration of effects? *BMJ
Evidence-Based Medicine*. 2025;30(1):e112747.
doi:[10.1136/bmjebm-2023-112747](https://doi.org/10.1136/bmjebm-2023-112747)
