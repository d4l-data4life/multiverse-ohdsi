library(dplyr)
library(ggplot2)

# null coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x

# define paths and create objects
setwd("C:/Users/MarioLawes/Desktop")
cm_path <- file.path(getwd(), "resultsFolder", "results_folder", "CohortMethodModule")
estimates <- read.csv(file.path(cm_path, "cm_result.csv"))
analysis <- read.csv(file.path(cm_path, "cm_analysis.csv"))

analysis
specDetails <- purrr::map_dfr(analysis$definition, \(json_str) {
  x <- jsonlite::fromJSON(json_str)
  
  tibble::tibble(
    analysis_id        = x$analysisId,
    description        = x$description,
    estimand          = x$createPsArgs$estimator              %||% NA,
    caliper            = x$matchOnPsArgs$caliper               %||% NA,
    maxRatio           = x$matchOnPsArgs$maxRatio              %||% NA,
    numberOfStrata     = x$stratifyByPsArgs$numberOfStrata     %||% NA,
    trimFraction       = x$trimByPsArgs$trimFraction           %||% NA,
    maxWeight          = x$truncateIptwArgs$maxWeight          %||% NA,
    inversePtWeighting = x$fitOutcomeModelArgs$inversePtWeighting %||% FALSE,
    stratified         = x$fitOutcomeModelArgs$stratified      %||% FALSE,
    washoutPeriod = x$createStudyPopArgs$washoutPeriod %||% NA
  )
})

table(specDetails$analysis_id)
table(estimates$analysis_id)

# join data
plot_data <- estimates %>%
  left_join(specDetails, by = "analysis_id") %>%
  arrange(rr) %>%                             
  mutate(rank = row_number(),
         significant = p < 0.05)

# --- Top panel: effect estimates ---
p_top <- ggplot(plot_data, aes(x = rank, y = rr, color = significant)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = ci_95_lb, ymax = ci_95_ub), 
                width = 0, alpha = 0.6) +
  geom_point(size = 3) +
  scale_color_manual(values = c("TRUE" = "#d62728", "FALSE" = "#1f77b4"),
                     labels = c("Not significant", "Significant"),
                     name = NULL) +
  scale_x_continuous(breaks = NULL,
                     expand = expansion(mult = 0.02)) +
  labs(y = "RR", x = NULL, title = "Specification Curve") +
  theme_bw() +
  theme(legend.position = "top",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# --- Bottom panel: specification grid ---
spec_long <- plot_data |>
  select(rank, washoutPeriod, estimand, caliper, maxRatio, numberOfStrata, trimFraction, maxWeight) |>
  mutate(
    washoutPeriod  = as.character(washoutPeriod),
    caliper        = as.character(caliper),
    maxRatio       = as.character(maxRatio),
    numberOfStrata = as.character(numberOfStrata),
    trimFraction   = as.character(trimFraction),
    maxWeight      = as.character(maxWeight)
    ) |>
  tidyr::pivot_longer(-rank, names_to = "spec_type", values_to = "spec_value") |>
  filter(!is.na(spec_value) & spec_value != "FALSE" & spec_value != "NA") |>
  mutate(group = case_when(
    spec_type %in% c("caliper", "maxRatio")          ~ "Matching",
    spec_type %in% c("maxWeight", "trimFraction")    ~ "Weighting",
    spec_type %in% c("numberOfStrata")               ~ "Stratification",
    spec_type %in% c("washoutPeriod")                ~ "General",
    spec_type %in% c("estimand")                    ~ "General",
    .default = "Other"
  ))

x_limits <- c(1, max(plot_data$rank))

panel_heights <- spec_long %>%
  group_by(group, spec_type) %>%
  summarise(n_unique = n_distinct(spec_value), .groups = "drop") %>%
  arrange(group, spec_type)  # must match facet order

panel_heights

p_bottom <- ggplot(spec_long, aes(x = rank, y = spec_value)) +
  geom_point(shape = 15, size = 3, color = "#2c7bb6") +
  ggh4x::facet_nested(group + spec_type ~ ., 
                      scales = "free_y",
                      space  = "free",        # <-- this does the proportional sizing
                      switch = "y",
                      nest_line = element_line(color = "grey70")) +
  scale_x_continuous(breaks = NULL, limits = x_limits,
                     expand = expansion(mult = 0.02)) +
  labs(x = "Specification (ranked by RR)", y = NULL) +
  theme_bw() +
  theme(
    strip.placement   = "outside",
    strip.background  = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    axis.text.x       = element_blank(),
    axis.ticks.x      = element_blank()
  )

plot_finished <- p_top / p_bottom + 
  plot_layout(heights = c(1, 2))

plot_finished
# 
# library(specr)
# 
# results <- run_specs(df = example_data, 
#                      y = c("y1", "y2"), 
#                      x = c("x1", "x2"), 
#                      model = c("lm"), 
#                      controls = c("c1", "c2"), 
#                      subsets = list(group1 = unique(example_data$group1)))
# 
# # Check
# results
