# =============================================================================
# multiverseOHDSI (minimal)
# -----------------------------------------------------------------------------
#   source("FUNS_multiverseOHDSI.R")
#
#   mv <- readMultiverse("runs/multiverse_20260828_141530/results_folder")
#
#   # 1. See what actually varies across the analyses, and in which combinations
#   inspectMultiverseSpec(mv)
#
#   # 2. Choose the arguments to show under the curve, and how to label them
#   createSpecificationCurve(
#     mv,
#     decisions = c("PS adjustment"      = "psAdjustment",
#                   "Outcome model"      = "outcomeModel",
#                   "Risk window (days)" = "createStudyPopulationArgs.riskWindowEnd",
#                   "Caliper"            = "matchOnPsArgs.caliper",
#                   "Match ratio"        = "matchOnPsArgs.maxRatio",
#                   "PS strata"          = "stratifyByPsArgs.numberOfStrata",
#                   "Trim fraction"      = "trimByPsArgs.trimFraction",
#                   "Max weight"         = "truncateIptwArgs.maxWeight")
#   )
#
#   createVolcanoPlot(mv, colourBy = c("PS adjustment" = "psAdjustment"))
#
# Both plot functions take `baseSize` to scale every font in the figure at
# once; raise it for posters and slides, lower it for print.
#
# Requires ggplot2 and jsonlite; patchwork for the stacked specification curve.
#
# -----------------------------------------------------------------------------
# readMultiverse(resultsFolder)
#   Reads cm_result.csv and joins it to a grid of every argument found in the
#   serialised cmAnalysisList. Columns keep their full dotted path
#   ("matchOnPsArgs.caliper"), so nothing is renamed or dropped behind your
#   back. Arguments absent from an analysis read "not applicable"; arguments
#   explicitly passed as NULL read "none".
#
#   Three kinds of derived column are added, describing how each specification
#   was built rather than what its arguments were:
#     psAdjustment          "Matching", "Stratification", "IPTW" or "None"
#     outcomeModel          "Cox", "Logistic", "Poisson", ...
#     uses.<block>          "yes"/"no" for each createXxxArgs() block, e.g.
#                           uses.matchOnPsArgs, uses.trimByPsArgs
#   All are usable in `decisions` and `colourBy` like any other column.
#
# inspectMultiverseSpec(x, all = FALSE)
#   Reports which arguments vary across analyses, their distinct values, and
#   the distinct combinations of those arguments. Use this to decide what
#   belongs under the plot, and to check the multiverse is the size you
#   declared. `all = TRUE` also lists the arguments held constant.
#
# createSpecificationCurve(x, decisions = NULL, baseSize = 14)
#   Estimates ranked smallest to largest, over a dashboard of the arguments
#   named in `decisions`. Pass a named vector to relabel the dashboard rows:
#     c("Caliper" = "matchOnPsArgs.caliper", "PS strata" = "stratifyByPsArgs.numberOfStrata")
#   Row order follows the order given. Defaults to every argument that varies.
#
# createVolcanoPlot(x, colourBy = NULL, baseSize = 14)
#   -log10(p) against the estimate. The region below p = 0.05 is shaded and
#   the 1st/50th/99th percentiles of the estimate are marked; the spread
#   between the outer two is the relative effect size ratio. Reference lines
#   and their labels are black, so colour carries only the `colourBy`
#   grouping. Specifications landing on the same point are collapsed into one
#   marker sized by how many they represent, so exact ties stay visible.
#   `colourBy` takes the same input as `decisions`, but one entry only:
#     colourBy = c("PS adjustment" = "psAdjustment")
#   The name becomes the legend title.
#
# Both plots put the estimate on a log-spaced axis, so 0.5 and 2 sit equally
# far from the null. The axis is labelled "Estimate" rather than with a named
# effect measure: the OHDSI results model stores every effect estimate in a
# column called `rr` whatever outcome model was fitted, so what the number
# means depends on the specification. Add `outcomeModel` to `decisions` to
# show which model produced each estimate.
# =============================================================================


# =============================================================================
# READ
# =============================================================================

readMultiverse <- function(resultsFolder) {
  resultFile <- .findOne(resultsFolder, "cm_result.csv")
  if (is.null(resultFile)) {
    stop("No 'cm_result.csv' found beneath '", resultsFolder, "'.")
  }
  
  results <- .readCsv(resultFile)
  
  grid <- .specificationGrid(
    jsonFile = .findOne(resultsFolder, "analysisSpecification.json"),
    csvFile  = .findOne(resultsFolder, "cm_analysis.csv")
  )
  
  estimates <- merge(results, grid$grid, by = "analysis_id", all.x = TRUE)
  estimates$significant <- !is.na(estimates$p) & estimates$p < 0.05
  estimates <- estimates[order(estimates$rr, na.last = TRUE), ]
  
  list(
    estimates  = estimates,
    grid       = grid$grid,        # every argument, full dotted paths
    varying    = grid$varying,     # those taking more than one value
    derived    = grid$derived,     # columns describing how each spec was built
    gridSource = grid$source
  )
}

.findOne <- function(root, pattern) {
  hits <- list.files(root, pattern = paste0("^", pattern, "$"),
                     recursive = TRUE, full.names = TRUE)
  if (length(hits) == 0) NULL else hits[1]
}

.readCsv <- function(path) {
  df <- utils::read.csv(path, stringsAsFactors = FALSE)
  names(df) <- tolower(names(df))
  df
}


# =============================================================================
# SPECIFICATION GRID
# =============================================================================
# Every argument in the serialised cmAnalysisList, one column per argument,
# one row per analysis. Preferred source is analysisSpecification.json, which
# holds the complete declared list; fallback is the `definition` column of
# cm_analysis.csv.

.specificationGrid <- function(jsonFile = NULL, csvFile = NULL) {
  if (!is.null(jsonFile) && file.exists(jsonFile)) {
    spec <- jsonlite::fromJSON(jsonFile, simplifyVector = FALSE)
    analyses <- .findElement(spec, "cmAnalysisList")
    if (!is.null(analyses)) {
      return(.assembleGrid(lapply(analyses, .flatten),
                           "analysisSpecification.json"))
    }
  }
  
  df <- .readCsv(csvFile)
  rows <- lapply(seq_len(nrow(df)), function(i) {
    flat <- .flatten(jsonlite::fromJSON(df$definition[i], simplifyVector = FALSE))
    flat[["analysisId"]] <- df$analysis_id[i]
    flat
  })
  .assembleGrid(rows, "cm_analysis.csv")
}

# Flatten a nested list to a named character vector. Names are dot-separated
# paths; NULL becomes "none", so "argument absent" stays a visible choice.
.flatten <- function(x, prefix = "") {
  if (is.null(x)) {
    return(stats::setNames("none", prefix))
  }
  if (!is.list(x)) {
    if (length(x) == 0) x <- "none"
    return(stats::setNames(paste(as.character(x), collapse = "|"), prefix))
  }
  nms <- names(x)
  if (is.null(nms)) nms <- as.character(seq_along(x))
  parts <- mapply(
    function(el, nm) {
      .flatten(el, if (nzchar(prefix)) paste(prefix, nm, sep = ".") else nm)
    },
    x, nms, SIMPLIFY = FALSE
  )
  # unname() before unlist(): the inner vectors already carry the full dotted
  # path, and unlist() would otherwise prepend the outer name again.
  unlist(unname(parts))
}

.assembleGrid <- function(rows, source) {
  allNames <- unique(unlist(lapply(rows, names)))
  mat <- do.call(rbind, lapply(rows, function(r) {
    v <- r[allNames]
    v[is.na(v)] <- "not applicable"
    stats::setNames(v, allNames)
  }))
  grid <- as.data.frame(mat, stringsAsFactors = FALSE)
  
  # ParallelLogger writes R class metadata into the JSON (e.g.
  # "matchOnPsArgs.attr_class"); these are serialisation artifacts, not
  # arguments, and are the one thing dropped without asking.
  grid <- grid[, !grepl("attr_|attributes|\\.class$|^class$", names(grid)),
               drop = FALSE]
  
  idCol <- grep("analysisId$", names(grid), value = TRUE)[1]
  grid$analysis_id <- as.integer(grid[[idCol]])
  
  fields <- setdiff(names(grid),
                    c(idCol, "analysis_id",
                      grep("description|analysisId", names(grid), value = TRUE)))
  varying <- fields[vapply(fields,
                           function(f) length(unique(grid[[f]])) > 1,
                           logical(1))]
  
  grid <- grid[, c("analysis_id", fields), drop = FALSE]
  grid <- grid[order(grid$analysis_id), ]
  
  # How each specification was built, not just what its arguments were.
  grid <- .derivedColumns(grid)
  derived <- setdiff(names(grid), c("analysis_id", fields))
  
  fields <- c(fields, derived)
  varying <- fields[vapply(fields,
                           function(f) length(unique(grid[[f]])) > 1,
                           logical(1))]
  
  list(grid = grid, varying = varying, derived = derived, source = source)
}

# Add columns describing how each specification was built, alongside the raw
# arguments. These are derived rather than read from the JSON:
#
#   uses.matchOnPsArgs, uses.stratifyByPsArgs, ...   yes / no
#   psAdjustment                                    Matching / Stratification /
#                                                   IPTW / None
#   outcomeModel                                    Cox / Logistic / Poisson
#
# psAdjustment rules: matchOnPsArgs present -> Matching; stratifyByPsArgs
# present -> Stratification; fitOutcomeModelArgs.inversePtWeighting TRUE ->
# IPTW. IPTW is detected from the weighting flag rather than from trimming or
# truncation, because an untrimmed, untruncated IPTW analysis carries neither
# of those blocks. Families are combined with " + " if more than one applies.
#
# outcomeModel reads fitOutcomeModelArgs.modelType and reports the model
# family only. Stratification is deliberately not folded in: a stratified Cox
# model is still a Cox model producing a hazard ratio, and the stratification
# is already visible through psAdjustment. outcomeModel matters because the
# results model stores every effect estimate in a column named `rr` whatever
# model was fitted, so the estimate's meaning can only be read from the
# specification.
#
# A block counts as absent when every one of its columns reads either
# "not applicable" (the argument does not appear for this analysis) or "none"
# (the block was serialised as an explicit JSON null). Both mean the block was
# not used; treating only the former as absent marks every block present on
# every analysis.
.derivedColumns <- function(grid) {
  fields <- setdiff(names(grid), "analysis_id")
  blocks <- unique(sub("\\..*$", "", fields))
  blocks <- blocks[grepl("Args$", blocks)]
  
  absent <- c("not applicable", "none")
  
  used <- list()
  for (b in blocks) {
    cols <- fields[sub("\\..*$", "", fields) == b]
    present <- apply(grid[, cols, drop = FALSE], 1,
                     function(r) any(!r %in% absent))
    used[[b]] <- present
    grid[[paste0("uses.", b)]] <- ifelse(present, "yes", "no")
  }
  
  isIptw <- if ("fitOutcomeModelArgs.inversePtWeighting" %in% fields) {
    grid[["fitOutcomeModelArgs.inversePtWeighting"]] %in% c("TRUE", "true", "1")
  } else {
    rep(FALSE, nrow(grid))
  }
  
  fam <- rep("", nrow(grid))
  addFamily <- function(fam, flag, label) {
    ifelse(flag, ifelse(nzchar(fam), paste(fam, label, sep = " + "), label), fam)
  }
  if (!is.null(used[["matchOnPsArgs"]])) {
    fam <- addFamily(fam, used[["matchOnPsArgs"]], "Matching")
  }
  if (!is.null(used[["stratifyByPsArgs"]])) {
    fam <- addFamily(fam, used[["stratifyByPsArgs"]], "Stratification")
  }
  fam <- addFamily(fam, isIptw, "IPTW")
  fam[!nzchar(fam)] <- "None"
  grid$psAdjustment <- fam
  
  if ("fitOutcomeModelArgs.modelType" %in% fields) {
    modelType <- grid[["fitOutcomeModelArgs.modelType"]]
    pretty <- c(cox = "Cox", logistic = "Logistic", poisson = "Poisson")
    label <- unname(pretty[modelType])
    label[is.na(label)] <- modelType[is.na(label)]
    grid$outcomeModel <- label
  }
  
  grid
}

# First element with the given name, depth first.
.findElement <- function(x, name) {
  if (!is.list(x)) return(NULL)
  if (!is.null(names(x)) && name %in% names(x)) return(x[[name]])
  for (el in x) {
    hit <- .findElement(el, name)
    if (!is.null(hit)) return(hit)
  }
  NULL
}


# =============================================================================
# INSPECT
# =============================================================================
# What varies, with which values, in which combinations. Read this before
# choosing what to put under the specification curve.

inspectMultiverseSpec <- function(x, all = FALSE) {
  if (is.character(x)) x <- readMultiverse(x)
  grid <- x$grid
  varying <- x$varying
  
  fields <- setdiff(names(grid), "analysis_id")
  show <- if (all) fields else varying
  
  values <- data.frame(
    argument = show,
    source   = ifelse(show %in% x$derived, "derived", "argument"),
    nValues  = vapply(show, function(f) length(unique(grid[[f]])), integer(1)),
    varies   = show %in% varying,
    values   = vapply(show, function(f) {
      paste(sort(unique(grid[[f]])), collapse = ", ")
    }, character(1)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  values <- values[order(-values$nValues, values$argument), ]
  
  # Distinct combinations of the varying arguments, with how many analyses
  # realise each. The row count is the size of the declared multiverse.
  combos <- NULL
  if (length(varying)) {
    cmb <- grid[, varying, drop = FALSE]
    key <- do.call(paste, c(cmb, sep = "\r"))
    counts <- table(key)
    keep <- !duplicated(key)
    combos <- cmb[keep, , drop = FALSE]
    combos$n <- as.integer(counts[key[keep]])
    combos <- combos[order(-combos$n), , drop = FALSE]
    row.names(combos) <- NULL
  }
  
  cat("Multiverse specification\n")
  cat("------------------------\n")
  cat("Source       :", x$gridSource, "\n")
  cat("Analyses     :", nrow(grid), "\n")
  cat("Arguments    :", length(fields) - length(x$derived), "read,",
      length(x$derived), "derived,", length(varying), "varying\n")
  if (!is.null(combos)) {
    cat("Combinations :", nrow(combos), "distinct\n")
  }
  cat("\nArguments",
      if (all) "(all)" else "(varying only; use all = TRUE for the rest)", "\n\n")
  print(values, right = FALSE)
  
  invisible(list(values = values, combinations = combos, grid = grid))
}


# =============================================================================
# VISUALISATION
# =============================================================================

# Shared preparation: drop non-estimable rows, order, rank.
.prepEstimates <- function(x) {
  if (is.character(x)) x <- readMultiverse(x)
  d <- x$estimates
  d <- d[!is.na(d$rr) & is.finite(d$rr) & d$rr > 0, ]
  d$.est <- d$rr
  d <- d[order(d$.est), ]
  d$.rank <- seq_len(nrow(d))
  list(data = d, spec = x)
}

# Resolve a decisions/colourBy vector into column names plus display labels.
# Names supply the labels; unnamed entries label themselves.
.resolveDecisions <- function(decisions, d) {
  cols <- unname(decisions)
  absent <- setdiff(cols, names(d))
  if (length(absent)) {
    stop("Not in the grid: ", paste(absent, collapse = ", "),
         "\nRun inspectMultiverseSpec() to list the available arguments.")
  }
  labels <- names(decisions)
  if (is.null(labels)) labels <- cols
  labels[!nzchar(labels)] <- cols[!nzchar(labels)]
  list(cols = cols, labels = labels)
}

# Readable breaks for a ratio axis: symmetric about 1, log-spaced, trimmed to
# the range actually present.
.ratioBreaks <- function(v) {
  v <- v[is.finite(v) & v > 0]
  if (!length(v)) return(1)
  candidates <- c(0.1, 0.125, 0.2, 0.25, 0.33, 0.5, 0.67, 0.8,
                  1, 1.25, 1.5, 2, 3, 4, 5, 8, 10)
  keep <- candidates >= min(v) * 0.95 & candidates <= max(v) * 1.05
  b <- candidates[keep]
  if (length(b) < 3) b <- signif(exp(seq(log(min(v)), log(max(v)),
                                         length.out = 5)), 2)
  sort(unique(c(b, 1)))
}

# Shared look for both figures: boxed panels, black axis text, horizontal
# gridlines only, and every font scaled from one number so the figure stays
# legible when it is shrunk into a poster or a slide.
.multiverseTheme <- function(baseSize = 14) {
  ggplot2::theme_bw(base_size = baseSize) +
    ggplot2::theme(
      panel.border       = ggplot2::element_rect(colour = "black", fill = NA,
                                                 linewidth = 0.8),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = "grey90",
                                                 linewidth = 0.4),
      axis.text          = ggplot2::element_text(colour = "black",
                                                 size = baseSize),
      axis.title         = ggplot2::element_text(size = baseSize + 2),
      axis.ticks         = ggplot2::element_line(colour = "black"),
      strip.background   = ggplot2::element_blank(),
      strip.placement    = "outside",
      strip.text         = ggplot2::element_text(size = baseSize,
                                                 colour = "black"),
      legend.position    = "top",
      legend.title       = ggplot2::element_text(size = baseSize),
      legend.text        = ggplot2::element_text(size = baseSize),
      legend.key         = ggplot2::element_blank()
    )
}


# ----------------------------------------------------- specification curve

createSpecificationCurve <- function(x, decisions = NULL, baseSize = 14) {
  prepped <- .prepEstimates(x)
  d <- prepped$data
  
  if (is.null(decisions)) decisions <- prepped$spec$varying
  dec <- .resolveDecisions(decisions, d)
  
  top <- ggplot2::ggplot(d, ggplot2::aes(x = .rank, y = .est)) +
    ggplot2::geom_linerange(
      ggplot2::aes(ymin = ci_95_lb, ymax = ci_95_ub,
                   colour = significant), linewidth = 0.7, alpha = 0.7) +
    ggplot2::geom_point(ggplot2::aes(colour = significant), size = 2.6) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed",
                        colour = "grey40", linewidth = 0.6) +
    ggplot2::scale_colour_manual(
      values = c("FALSE" = "#2E75B6", "TRUE" = "#E03531"),
      labels = c("p >= 0.05", "p < 0.05"), name = NULL) +
    ggplot2::scale_y_continuous(
      trans = "log", breaks = .ratioBreaks(c(d$ci_95_lb, d$ci_95_ub, d$.est))) +
    ggplot2::labs(x = NULL, y = "Estimate") +
    .multiverseTheme(baseSize) +
    ggplot2::theme(axis.text.x  = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.text.y = element_text(size = baseSize - 4))
  
  long <- do.call(rbind, lapply(seq_along(dec$cols), function(i) {
    data.frame(rank = d$.rank, decision = dec$labels[i],
               value = as.character(d[[dec$cols[i]]]),
               stringsAsFactors = FALSE)
  }))
  long <- long[long$value != "not applicable", ]
  long$decision <- factor(long$decision, levels = dec$labels)
  
  bottom <- ggplot2::ggplot(long, ggplot2::aes(x = rank, y = value)) +
    ggplot2::geom_point(shape = 15, size = 2.6, colour = "#2E75B6") +
    ggplot2::facet_grid(decision ~ ., scales = "free_y", space = "free_y",
                        switch = "y") +
    ggplot2::labs(x = "Specification (ranked by estimate)", y = NULL) +
    .multiverseTheme(baseSize) +
    ggplot2::theme(
      axis.text.x       = ggplot2::element_blank(),
      axis.ticks.x      = ggplot2::element_blank(),
      strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1,
                                                size = baseSize)
    )
  
  patchwork::wrap_plots(top, bottom, ncol = 1, heights = c(1.5, 3.5))
}


# ------------------------------------------------------------------ volcano

createVolcanoPlot <- function(x, colourBy = NULL, baseSize = 14) {
  prepped <- .prepEstimates(x)
  d <- prepped$data
  
  d$.negLogP <- -log10(pmax(d$p, .Machine$double.xmin))
  sigLine <- -log10(0.05)
  
  legendTitle <- NULL
  if (!is.null(colourBy)) {
    if (length(colourBy) > 1) {
      stop("colourBy takes one argument only; got ", length(colourBy), ".")
    }
    cb <- .resolveDecisions(colourBy, d)
    d$.colour <- as.character(d[[cb$cols]])
    legendTitle <- cb$labels
  }
  
  # Specifications with identical estimates overplot exactly, so a cluster of
  # a dozen is indistinguishable from a single point. Count them and keep one
  # row per distinct position, sizing the marker by the count. Points are
  # separated by colour group as well as position, so a group never absorbs
  # another group's ties.
  key <- paste(signif(d$.est, 8), signif(d$.negLogP, 8),
               if (is.null(colourBy)) "" else d$.colour)
  d$.n <- as.integer(table(key)[key])
  d <- d[!duplicated(key), ]
  
  q <- stats::quantile(d$.est, probs = c(0.01, 0.5, 0.99), na.rm = TRUE)
  qDf <- data.frame(x = as.numeric(q), label = c("1st", "50th", "99th"))
  
  mapping <- if (is.null(colourBy)) {
    ggplot2::aes(x = .est, y = .negLogP, size = .n)
  } else {
    ggplot2::aes(x = .est, y = .negLogP, colour = .colour, size = .n)
  }
  
  ggplot2::ggplot(d, mapping) +
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf,
                      ymin = -Inf, ymax = sigLine,
                      fill = "grey50", alpha = 0.10) +
    ggplot2::geom_vline(data = qDf, ggplot2::aes(xintercept = x),
                        colour = "black", linetype = "dotted",
                        linewidth = 0.6, inherit.aes = FALSE) +
    ggplot2::geom_text(data = qDf, ggplot2::aes(x = x, y = Inf, label = label),
                       inherit.aes = FALSE, vjust = 1.4, hjust = -0.15,
                       size = baseSize / 3.2, colour = "black") +
    ggplot2::geom_hline(yintercept = sigLine, colour = "black",
                        linetype = "dashed", linewidth = 0.6) +
    ggplot2::geom_vline(xintercept = 1, colour = "black",
                        linetype = "longdash", linewidth = 0.6) +
    ggplot2::geom_point(alpha = 0.75) +
    ggplot2::annotate("text", x = -Inf, y = sigLine, label = "p = 0.05",
                      hjust = -0.15, vjust = -0.6, size = baseSize / 3.2,
                      colour = "black") +
    ggplot2::scale_x_continuous(trans = "log", breaks = .ratioBreaks(d$.est)) +
    ggplot2::scale_size_area(max_size = 10,
                             breaks = function(v) unique(round(pretty(v)))) +
    ggplot2::labs(x = "Estimate", y = expression(-log[10](p)),
                  colour = legendTitle, size = "Specifications") +
    .multiverseTheme(baseSize)
}