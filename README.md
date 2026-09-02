# Multiverse analysis for OHDSI

This repo contains template code to conduct a multiverse analysis using
[Strategus](https://github.com/OHDSI/Strategus) and the wider
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

Harmonising the data and the tooling reduces variability in the implementation of
real-world evidence studies, and best practice settles many design questions. However, it
does not settle all of them. How long the risk window should be, whether to match,
stratify, or weight on the propensity score, how wide a caliper to set, how many strata to
use, whether to trim or truncate extreme weights — guidance favours no single option, yet
each can move the estimate. Multiverse analysis addresses precisely this (Steegen et al.,
2016): instead of reporting one pipeline, the analyst defines all defensible alternatives,
executes them, and reports the distribution of findings. Related approaches include
specification curve analysis (Simonsohn et al., 2020) and, from epidemiology, vibration of
effects (Patel et al., 2015; Vinatier et al., 2025).

Importantly, multiverse analysis is not a single method. Recent work distinguishes three
variants with different purposes (Lemster et al., 2026): a confirmatory multiverse
restricts itself to a small set of specifications targeting a precisely defined estimand
and draws an inferential conclusion with adjustment for multiplicity; an exploratory
multiverse admits a broader but still defensible set, foregoes formal inference, and serves
robustness assessment and hypothesis generation; a meta-research multiverse takes the
research process rather than the clinical question as its object, and is the only variant
in which specifications that are not equally defensible may legitimately be included,
precisely because researchers demonstrably use them. So framed, multiverse analysis
complements rather than duplicates existing large-scale evidence generation in OHDSI:
programmes such as LEGEND vary the question across many exposure–comparator–outcome
triplets and databases while holding the design close to fixed (Suchard et al., 2019),
whereas a multiverse varies the design while holding the question fixed.

Building on OHDSI methodology, we propose a structured framework for multiverse analysis in
large-scale RWE and give guidance on how to interpret the resulting distribution of
estimates.

## Illustrative example

The framework is illustrated on the synthetic
[Eunomia](https://github.com/OHDSI/Eunomia) GiBleed demo CDM.

| | |
|---|---|
| **Target** | Celecoxib (cohort 1) |
| **Comparator** | Diclofenac (cohort 2) |
| **Outcome** | Gastrointestinal bleed (cohort 3) |
| **Data** | Eunomia GiBleed synthetic CDM |
| **Estimand** | ATT |
| **Outcome model** | Cox proportional hazards |

Covariates for the two exposure ingredients (and their descendants) are excluded from
the propensity model.

This multiverse is exploratory in the sense of Lemster et al. (2026): the specification set
is broad rather than minimal, and no inferential conclusion is drawn from it.

### Specification space

The 42 specifications are the crossing of three adjustment families with three risk windows:

| Family | Varied arguments | n |
|---|---|---|
| PS matching | `caliper` ∈ {0.2, 0.0001} × `maxRatio` ∈ {1, 10} | 4 |
| PS stratification | `numberOfStrata` ∈ {4, 5, 6, 7, 8, 9} | 6 |
| IPTW | `trimFraction` ∈ {0, 0.01} × `maxWeight` ∈ {none, 10} | 4 |
| | **Subtotal** | **14** |
| | **Risk window** ∈ {30, 60, 90} days | ×3 |
| | **Total** | **42** |

Analysis IDs map as follows: 1–12 matching, 13–30 stratification, 31–42 weighting.

## Repository contents

```
1 - runMultiverseAnalysis.R     Build the 42-specification Strategus analysis
                                specification and execute it against Eunomia
2 - extractResults_visualize.R  Parse the executed specifications back out,
                                join them to the estimates, and draw the
                                specification curve and volcano plot
FUNS_multiverseOHDSI.R          Sourceable helpers behind script 2:
                                readMultiverse(), inspectMultiverseSpec(),
                                createSpecificationCurve(), createVolcanoPlot()
renv.lock                       Pinned package versions (see Requirements)
```

Each execution writes to its own timestamped folder, so runs accumulate rather
than overwrite:

```
runs/multiverse_<YYYYMMDD_HHMMSS>/
├── analysisSpecification.json  Machine-readable declaration of the full
│                               multiverse, written before execution begins;
│                               the pre-registerable artifact
├── work_folder/
└── results_folder/
    └── CohortMethodModule/
        ├── cm_analysis.csv     One row per specification; `definition` holds
        │                       the analysis as JSON
        └── cm_result.csv       Effect estimates, CIs, p-values per analysis
```

## Requirements

R (4.4.1) with Strategus, CohortGenerator, CohortMethod, FeatureExtraction and Eunomia for
execution, and ggplot2, jsonlite and patchwork for extraction and plotting.

`renv` is used for dependency pinning; restore the recorded library with:

```r
renv::restore()
```

## Running the analysis

```r
source("1 - runMultiverseAnalysis.R")     # builds and executes all 42 specifications
# then set `resultsFolder` in the script below to the run folder just created
source("2 - extractResults_visualize.R")  # extracts results and plots
```

The first script downloads and instantiates the Eunomia CDM, creates the demo cohorts,
writes `analysisSpecification.json`, and calls `Strategus::execute()`.

The second script recovers the varied arguments from the serialised `cmAnalysisList` —
preferably `analysisSpecification.json`, falling back to the `definition` JSON column of
`cm_analysis.csv` — rather than re-declaring them in the plotting code, so the dashboard
cannot drift from the analysis list that was submitted. Because each run gets its own
timestamped folder, the second script has to be pointed at the run you want:

```r
resultsFolder <- file.path(getwd(), "runs", "multiverse_<YYYYMMDD_HHMMSS>", "results_folder")
```

`mv$gridSource` records which file was used. `inspectMultiverseSpec()` then prints every
argument that varies, its distinct values, and the distinct combinations realised, which is
also the check that the executed multiverse is the size that was declared.

Two visualisations are available: (1) a specification curve, with ranked estimates and 95%
CIs on top and the corresponding specification grid below, and (2) a volcano plot showing
−log10(p) against the estimate, coloured by PS adjustment family, with the 1st, 50th and
99th percentiles of the estimate marked.

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

Lemster S, Bonneville EF, Castro CAP, Columbus A, Dunkler D, Schmidt CO, Hothorn A,
Hoffmann S, Boulesteix A-L. Multiverse-style analyses for confirmatory, exploratory, and
meta-research purposes: design, execution, and interpretation. 2026. Preprint.
<https://www.ibe.med.uni-muenchen.de/mitarbeiter/mitarbeiter/simon-lemster/preprint_multiverse.pdf>

Levitt M, Zonta F, Ioannidis JPA. Excess death estimates from multiverse analysis in
2009–2021. *European Journal of Epidemiology*. 2023;38(11):1129–1139.
doi:[10.1007/s10654-023-00998-2](https://doi.org/10.1007/s10654-023-00998-2)

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

Suchard MA, Schuemie MJ, Krumholz HM, et al. Comprehensive comparative effectiveness and
safety of first-line antihypertensive drug classes: a systematic, multinational,
large-scale analysis. *The Lancet*. 2019;394(10211):1816–1827.
doi:[10.1016/S0140-6736(19)32317-7](https://doi.org/10.1016/S0140-6736(19)32317-7)

Vinatier C, Hoffmann S, Patel C, et al. What is the vibration of effects? *BMJ
Evidence-Based Medicine*. 2025;30(1):e112747.
doi:[10.1136/bmjebm-2023-112747](https://doi.org/10.1136/bmjebm-2023-112747)