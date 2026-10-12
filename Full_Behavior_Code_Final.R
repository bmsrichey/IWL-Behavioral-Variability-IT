###############################################################################
# Invasion-related behavior analyses and manuscript figures
#
# Cleaned, frozen data. Runs statistical analyses and generates figures.
###############################################################################

# =============================================================================
# Install and load packages
# =============================================================================

# List the packages required for the final statistical analyses and figures
required_packages <- c("tidyverse", "lme4", "car", "emmeans", "rptR", "gridExtra")

# Install packages that are not currently available
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

# Load required packages
invisible(lapply(required_packages, library, character.only = TRUE))

# Use sum contrasts for Type-III tests of repeated behaviors
options(scipen = 999, contrasts = c("contr.sum", "contr.poly"))

# Random-number generation for Monte Carlo Fisher and repeatability tests
set.seed(123)

# Prevent a blank Rplots.pdf from being created by noninteractive graphics calls.
if (!interactive()) grDevices::pdf(NULL)

# Create the two directories used for manuscript-ready outputs
dir.create("results", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)

# Register frozen data files
analysis_files <- c(individual = "frozen_data/final_data_consensus_individual.csv", boldness = "frozen_data/bold_summary_recorder_level.csv", 
  exploration = "frozen_data/exp_summary_recorder_level.csv", repeated = "frozen_data/long_flex_data_consensus.csv")

# Fix population order
model_population_order <- c("cal", "cap", "elb", "gor", "liv")

# Fix assay order
assay_order <- c("Aggression", "Boldness", "Neophilia", "Exploration")

# =============================================================================
# Import analysis-ready data
# =============================================================================

# Convert identifiers and design variables to the factor types used by models
prepare_analysis_data <- function(data, factor_columns = character()) {
  data %>% dplyr::mutate(id = as.character(id), dplyr::across(dplyr::any_of(c("sex_binary", "sex_ternary", factor_columns)), 
    as.factor), pop_id = factor(pop_id, levels = model_population_order))
}

# Read one consensus record per individual for correlations and aggression
individual_data <- readr::read_csv(analysis_files[["individual"]], show_col_types = FALSE) %>% prepare_analysis_data(c("nba_day_agg", 
  "nba_order_agg", "group_number_agg", "nba_day_bold", "nba_order_bold", "group_number_bold", "nba_day_neo", 
  "nba_order_neo", "group_number_neo", "nba_day_exp", "nba_order_exp", "group_number_exp"))

# Read recorder-level boldness scores so technical replication remains clustered
boldness_data <- readr::read_csv(analysis_files[["boldness"]], show_col_types = FALSE) %>% prepare_analysis_data(c("nba_order_bold", 
  "nba_day_bold", "group_number_bold", "recorder_bold"))

# Read recorder-level exploration scores so technical replication remains clustered
exploration_data <- readr::read_csv(analysis_files[["exploration"]], show_col_types = FALSE) %>% prepare_analysis_data(c("nba_order_exp", 
  "nba_day_exp", "group_number_exp", "recorder_exp"))

# Read consensus climbing and sheltering proportions for all assays
repeated_data <- readr::read_csv(analysis_files[["repeated"]], show_col_types = FALSE) %>% prepare_analysis_data(c("behavioral_assay", 
  "nba_order", "nba_day", "group_number")) %>% dplyr::mutate(behavioral_assay = factor(behavioral_assay, levels = assay_order))

# =============================================================================
# Prepare analysis variables
# =============================================================================

# Standardize SVL once so coefficients are comparable across models
svl_center <- mean(individual_data$svl, na.rm = TRUE)
svl_scale <- stats::sd(individual_data$svl, na.rm = TRUE)

# Describe the one-standard-deviation SVL scale used in every retained model
svl_standardization_note <- paste0("SVL standardized at mean ", format(signif(svl_center, 5), trim = TRUE), " and SD ", 
  format(signif(svl_scale, 5), trim = TRUE), "; this test represents a one-SD increase in SVL.")

# Apply the SVL center and scale to every analysis table
add_standardized_svl <- function(data) {
  dplyr::mutate(data, svl_z = (svl - svl_center)/svl_scale)
}

# Create the boldness score used only for directional correlations
individual_data <- individual_data %>% add_standardized_svl() %>% dplyr::mutate(boldness_score = -log1p(bold_latency))

# Log-transform emergence latency for the Gaussian boldness model
boldness_data <- boldness_data %>% add_standardized_svl() %>% dplyr::mutate(log_bold_latency = log1p(bold_latency))

# Count corners visited and unvisited among three available corners
exploration_data <- exploration_data %>% add_standardized_svl() %>% dplyr::mutate(exp_unvisited_corners = 3L - 
  exp_unique_corners)

# Bound proportions at 0.001 and 0.999 so logit values remain finite
epsilon <- 0.001
repeated_data <- repeated_data %>% add_standardized_svl() %>% dplyr::mutate(climbing_prop_adjusted = pmin(pmax(climbing_prop, 
  epsilon), 1 - epsilon), shelter_prop_adjusted = pmin(pmax(shelter_prop, epsilon), 1 - epsilon), logit_climbing = stats::qlogis(climbing_prop_adjusted), 
  logit_sheltering = stats::qlogis(shelter_prop_adjusted), id = as.factor(id))

# =============================================================================
# Statistical helper functions
# =============================================================================

# Return the first available column name
pick_column <- function(data, candidates, label) {
  found <- intersect(candidates, names(data))
  if (length(found) == 0) {
    stop("Could not identify ", label, " column.", call. = FALSE)
  }
  found[[1]]
}

# Store a model formula as one readable string for the consolidated table
formula_text <- function(model) {
  paste(deparse(stats::formula(model)), collapse = " ")
}

# Extract the variance attributed to one random grouping factor
random_variance <- function(model, group) {
  variance_table <- as.data.frame(lme4::VarCorr(model))
  values <- variance_table %>% dplyr::filter(.data$grp == .env$group) %>% dplyr::pull(vcov)
  values <- values[is.finite(values)]
  if (length(values) == 0) 
    NA_real_
  else values[[1]]
}

# Convert Type-II focal-model tests to a common manuscript-facing structure
focal_anova_rows <- function(model, response, outcome) {
  quasi_model <- inherits(model, "glm") && startsWith(stats::family(model)$family, "quasi")
  test_table <- if (quasi_model) {
    car::Anova(model, type = "II", test.statistic = "F")
  } else {
    car::Anova(model, type = "II")
  }
  test_table <- test_table %>% as.data.frame() %>% tibble::rownames_to_column("term") %>% tibble::as_tibble()
  statistic_column <- pick_column(test_table, c("Chisq", "F value", "F"), "test statistic")
  df_column <- pick_column(test_table, c("Df", "num Df"), "degrees of freedom")
  p_column <- pick_column(test_table, grep("^Pr\\(", names(test_table), value = TRUE), "P value")
  tibble::tibble(analysis_family = "Focal behavior", response, outcome, term = test_table$term, test_method = if (quasi_model) {
    "Type-II quasi-likelihood F test"
  }
  else if (statistic_column == "Chisq") {
    "Type-II Wald chi-square"
  } else {
    "Type-II F test"
  }, statistic = as.numeric(test_table[[statistic_column]]), df = as.numeric(test_table[[df_column]]), p_value = as.numeric(test_table[[p_column]]), 
    model = formula_text(model), note = NA_character_)
}

# Run a Monte Carlo Fisher exact test for a sparse categorical comparison
safe_fisher_test <- function(x, y, replicates = 100000L) {
  complete <- stats::complete.cases(x, y)
  contingency <- table(droplevels(as.factor(x[complete])), y[complete])
  if (nrow(contingency) < 2 || ncol(contingency) < 2) {
    return(list(p.value = NA_real_, method = "Not estimable because fewer than two levels were observed"))
  }
  stats::fisher.test(contingency, simulate.p.value = TRUE, B = replicates)
}

# Convert estimated marginal means and confidence limits to one table
estimated_mean_table <- function(emm, response, scale, transform = identity, multiplier = 1) {
  table <- as.data.frame(stats::confint(emm, level = 0.95))
  estimate_column <- pick_column(table, c("emmean", "response", "prob", "rate"), "estimated mean")
  lower_column <- pick_column(table, c("lower.CL", "asymp.LCL"), "lower confidence limit")
  upper_column <- pick_column(table, c("upper.CL", "asymp.UCL"), "upper confidence limit")
  tibble::tibble(analysis_family = "Focal behavior", response, pop_id = as.character(table$pop_id), behavioral_assay = NA_character_, 
    estimate = transform(table[[estimate_column]]) * multiplier, lower_95 = transform(table[[lower_column]]) * 
      multiplier, upper_95 = transform(table[[upper_column]]) * multiplier, scale, source = "Model-adjusted estimated marginal mean", 
    n = NA_integer_)
}

# Create an empty contrast table with the final required column types
empty_contrast_table <- function() {
  tibble::tibble(analysis_family = character(), response = character(), effect = character(), conditioning = character(), 
    contrast = character(), estimate = double(), lower_95 = double(), upper_95 = double(), scale = character(), 
    statistic = double(), df = double(), p_value = double(), p_holm = double(), note = character())
}

# Extract pairwise contrasts and apply Holm correction within each contrast family
pairwise_table <- function(emm, analysis_family, response, effect, scale, conditioning_column = NULL, transform = identity, 
  note = NA_character_) {
  contrasts <- emmeans::contrast(emm, method = "pairwise", adjust = "none")
  table <- as.data.frame(summary(contrasts, infer = c(TRUE, TRUE), adjust = "none"))
  lower_column <- pick_column(table, c("lower.CL", "asymp.LCL"), "contrast lower confidence limit")
  upper_column <- pick_column(table, c("upper.CL", "asymp.UCL"), "contrast upper confidence limit")
  statistic_column <- pick_column(table, c("t.ratio", "z.ratio"), "contrast statistic")
  conditioning <- if (is.null(conditioning_column) || !conditioning_column %in% names(table)) {
    rep(NA_character_, nrow(table))
  } else {
    paste0(conditioning_column, "=", table[[conditioning_column]])
  }
  tibble::tibble(analysis_family, response, effect, conditioning, contrast = as.character(table$contrast), estimate = transform(table$estimate), 
    lower_95 = transform(table[[lower_column]]), upper_95 = transform(table[[upper_column]]), scale, statistic = as.numeric(table[[statistic_column]]), 
    df = if ("df" %in% names(table)) 
      as.numeric(table$df)
    else NA_real_, p_value = as.numeric(table$p.value), p_holm = stats::p.adjust(as.numeric(table$p.value), 
      method = "holm"), note)
}

# Retrieve one  term P value from a model-test table
term_p_value <- function(table, response_name, term_name) {
  values <- table %>% dplyr::filter(response == response_name, term == term_name) %>% dplyr::pull(p_value)
  if (length(values) == 1) 
    values[[1]]
  else NA_real_
}

# Calculate six pairwise Spearman correlations among four behaviors and Holm-adjust within groups
correlation_table <- function(data, variables, group_label) {
  purrr::map_dfr(utils::combn(variables, 2, simplify = FALSE), function(pair) {
    pair_data <- data %>% dplyr::select(dplyr::all_of(pair)) %>% dplyr::filter(stats::complete.cases(.))
    estimable <- nrow(pair_data) >= 3 && stats::sd(pair_data[[pair[[1]]]]) > 0 && stats::sd(pair_data[[pair[[2]]]]) > 
      0
    if (!estimable) {
      return(tibble::tibble(Group = group_label, Var1 = pair[[1]], Var2 = pair[[2]], Correlation = NA_real_, 
        P_Value = NA_real_, N = nrow(pair_data)))
    }
    test <- stats::cor.test(pair_data[[pair[[1]]]], pair_data[[pair[[2]]]], method = "spearman", exact = FALSE)
    tibble::tibble(Group = group_label, Var1 = pair[[1]], Var2 = pair[[2]], Correlation = unname(test$estimate), 
      P_Value = test$p.value, N = nrow(pair_data))
  }) %>% dplyr::mutate(P_Holm = stats::p.adjust(P_Value, method = "holm"))
}

# Test whether population-specific correlations differ in Fisher z space
correlation_homogeneity <- function(correlation_data) {
  correlation_data %>% dplyr::distinct(Var1, Var2) %>% purrr::pmap_dfr(function(Var1, Var2) {
    pair_data <- correlation_data %>% dplyr::filter(.data$Var1 == .env$Var1, .data$Var2 == .env$Var2, !is.na(Correlation), 
      N > 3, abs(Correlation) < 1)
    if (nrow(pair_data) < 2) {
      return(tibble::tibble(Var1, Var2, k_populations_used = nrow(pair_data), populations_used = paste(pair_data$Group, 
        collapse = "; "), Q = NA_real_, df = NA_integer_, P_Value = NA_real_, note = "Fewer than two estimable population correlations"))
    }
    fisher_z <- atanh(pair_data$Correlation)
    weights <- pair_data$N - 3
    weighted_mean <- sum(weights * fisher_z)/sum(weights)
    q_statistic <- sum(weights * (fisher_z - weighted_mean)^2)
    degrees_freedom <- nrow(pair_data) - 1L
    tibble::tibble(Var1, Var2, k_populations_used = nrow(pair_data), populations_used = paste(pair_data$Group, 
      collapse = "; "), Q = q_statistic, df = degrees_freedom, P_Value = stats::pchisq(q_statistic, degrees_freedom, 
      lower.tail = FALSE), note = "Approximate Fisher r-to-z homogeneity test")
  })
}

# Use rank-based interactions as sensitivity tests for population differences
rank_interaction_tests <- function(data, variables) {
  purrr::map_dfr(utils::combn(variables, 2, simplify = FALSE), function(pair) {
    analysis_data <- data %>% dplyr::select(pop_id, dplyr::all_of(pair)) %>% dplyr::filter(stats::complete.cases(.)) %>% 
      droplevels()
    estimable <- nrow(analysis_data) >= 10 && nlevels(as.factor(analysis_data$pop_id)) >= 2 && stats::sd(analysis_data[[pair[[1]]]]) > 
      0 && stats::sd(analysis_data[[pair[[2]]]]) > 0
    if (!estimable) {
      return(tibble::tibble(Var1 = pair[[1]], Var2 = pair[[2]], N = nrow(analysis_data), F = NA_real_, df = NA_real_, 
        P_Value = NA_real_, note = "Insufficient data, population levels, or response variation"))
    }
    analysis_data <- analysis_data %>% dplyr::mutate(pop_id = droplevels(as.factor(pop_id)), x_rank = rank(.data[[pair[[1]]]], 
      ties.method = "average"), y_rank = rank(.data[[pair[[2]]]], ties.method = "average"))
    interaction <- car::Anova(stats::lm(y_rank ~ x_rank * pop_id, data = analysis_data), type = "II") %>% as.data.frame() %>% 
      tibble::rownames_to_column("term") %>% dplyr::filter(term == "x_rank:pop_id")
    tibble::tibble(Var1 = pair[[1]], Var2 = pair[[2]], N = nrow(analysis_data), F = if (nrow(interaction) == 
      1) 
      interaction$`F value`
    else NA_real_, df = if (nrow(interaction) == 1) 
      interaction$Df
    else NA_real_, P_Value = if (nrow(interaction) == 1) {
      interaction$`Pr(>F)`
    } else {
      NA_real_
    }, note = if (nrow(interaction) == 1) {
      "Rank-transformed behavior-by-population interaction"
    } else {
      "Interaction term was not estimable"
    })
  })
}

# =============================================================================
# Question 1: Focal-behavior correlations
# =============================================================================

# Define four focal variables 
correlation_variables <- c("agg_bite_freq", "boldness_score", "neo_freq", "exp_unique_corners")

# Estimate correlations across all individuals
correlations_overall <- correlation_table(individual_data, correlation_variables, "All individuals")

# Estimate correlations within each population
correlations_by_population <- purrr::map_dfr(model_population_order, function(population) {
  correlation_table(dplyr::filter(individual_data, pop_id == population), correlation_variables, population)
})

# Combine overall and population estimates for tables and heatmap
correlation_results <- dplyr::bind_rows(correlations_overall, correlations_by_population) %>% dplyr::mutate(scope = ifelse(Group == 
  "All individuals", "Overall", "Within population"), direction_note = ifelse(Var1 == "boldness_score" | Var2 == 
  "boldness_score", "Higher boldness_score means shorter emergence latency", NA_character_), .before = 1)

# Compare population-specific correlations using Fisher z tests
correlation_fisher_tests <- correlation_homogeneity(correlations_by_population) %>% dplyr::mutate(P_Holm = stats::p.adjust(P_Value, 
  method = "holm"))

# Run complementary rank-interaction tests for population comparisons
correlation_rank_tests <- rank_interaction_tests(individual_data, correlation_variables) %>% dplyr::mutate(P_Holm = stats::p.adjust(P_Value, 
  method = "holm"))

# Combine both population-comparison approaches
correlation_population_tests <- dplyr::full_join(correlation_fisher_tests %>% dplyr::transmute(Var1, Var2, k_populations_used, 
  populations_used, fisher_Q = Q, fisher_df = df, fisher_p = P_Value, fisher_p_holm = P_Holm, fisher_note = note), 
  correlation_rank_tests %>% dplyr::transmute(Var1, Var2, rank_model_N = N, rank_interaction_F = F, rank_interaction_df = df, 
    rank_interaction_p = P_Value, rank_interaction_p_holm = P_Holm, rank_interaction_note = note), by = c("Var1", 
    "Var2"))

# =============================================================================
# Question 2: Focal behaviors among populations
# =============================================================================

# Set stable optimizer and iteration limit for mixed models
lmer_control <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+05))
glmer_control <- lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+05))

# Define tolerance for effectively zero recorder variance
recorder_variance_tolerance <- 1e-08

# Fit boldness with individual and recorder intercepts
boldness_candidate <- lme4::lmer(log_bold_latency ~ pop_id + sex_ternary + svl_z + nba_order_bold + (1 | id) + 
  (1 | recorder_bold), data = boldness_data, REML = TRUE, control = lmer_control)

# Measure recorder variance before selecting final boldness model
boldness_recorder_variance <- random_variance(boldness_candidate, "recorder_bold")

# Remove recorder intercept only when variance is effectively zero
boldness_model <- if (!is.finite(boldness_recorder_variance) || boldness_recorder_variance < recorder_variance_tolerance) {
  lme4::lmer(log_bold_latency ~ pop_id + sex_ternary + svl_z + nba_order_bold + (1 | id), data = boldness_data, 
    REML = TRUE, control = lmer_control)
} else {
  boldness_candidate
}

# Use one equal-weight recorder-consensus neophilia response per individual
neophilia_individual <- individual_data %>% dplyr::filter(stats::complete.cases(id, neo_freq, pop_id, sex_ternary, 
  svl_z, nba_order_neo)) %>% droplevels()

# Fit quasi-Poisson neophilia because consensus means are fractional and overdispersed
neophilia_model <- stats::glm(neo_freq ~ pop_id + sex_ternary + svl_z + nba_order_neo, data = neophilia_individual, 
  family = stats::quasipoisson(link = "log"))

# Model corners visited and unvisited among three available corners
exploration_candidate <- lme4::glmer(cbind(exp_unique_corners, exp_unvisited_corners) ~ pop_id + sex_ternary + 
  svl_z + nba_day_exp + (1 | id) + (1 | recorder_exp), data = exploration_data, family = binomial, control = glmer_control)

# Measure recorder variance before selecting final exploration model
exploration_recorder_variance <- random_variance(exploration_candidate, "recorder_exp")

# Remove recorder intercept only when variance is effectively zero
exploration_model <- if (!is.finite(exploration_recorder_variance) || exploration_recorder_variance < recorder_variance_tolerance) {
  lme4::glmer(cbind(exp_unique_corners, exp_unvisited_corners) ~ pop_id + sex_ternary + svl_z + nba_day_exp + 
    (1 | id), data = exploration_data, family = binomial, control = glmer_control)
} else {
  exploration_candidate
}

# Record recorder treatment for each focal model
boldness_recorder_note <- if (is.finite(boldness_recorder_variance) && boldness_recorder_variance >= recorder_variance_tolerance) {
  paste0("Recorder random intercept retained; estimated variance = ", format(signif(boldness_recorder_variance, 
    4), scientific = TRUE))
} else {
  paste0("Recorder random intercept omitted because its variance was not estimable ", "or was below ", recorder_variance_tolerance)
}
exploration_recorder_note <- if (is.finite(exploration_recorder_variance) && exploration_recorder_variance >= recorder_variance_tolerance) {
  paste0("Recorder random intercept retained; estimated variance = ", format(signif(exploration_recorder_variance, 
    4), scientific = TRUE))
} else {
  paste0("Recorder random intercept omitted because its variance was not estimable ", "or was below ", recorder_variance_tolerance)
}

# Extract population, sex, SVL, and scheduling tests from focal models
focal_model_tests <- dplyr::bind_rows(focal_anova_rows(boldness_model, "Boldness", "log1p emergence latency; lower values indicate greater boldness"), 
  focal_anova_rows(neophilia_model, "Neophilia", "mean recorder-consensus novel-object contacts"), focal_anova_rows(exploration_model, 
    "Exploration", "unique corners visited out of three")) %>% dplyr::mutate(note = dplyr::case_when(response == 
  "Boldness" ~ boldness_recorder_note, response == "Neophilia" ~ paste("One equal-weight recorder-consensus response per individual;", 
  "quasi-Poisson variance accounts for overdispersion"), response == "Exploration" ~ exploration_recorder_note, 
  TRUE ~ note))

# Convert aggression to whether each individual bit at least once
aggression_individual <- individual_data %>% dplyr::filter(!is.na(agg_bite_freq)) %>% dplyr::mutate(agg_bite_present = as.integer(agg_bite_freq > 
  0))

# Define prespecified aggression comparisons
aggression_test_variables <- list(pop_id = aggression_individual$pop_id, sex_ternary = aggression_individual$sex_ternary, 
  nba_order_agg = aggression_individual$nba_order_agg)

# Use Monte Carlo Fisher tests because bites are too sparse for GLMM
aggression_tests <- purrr::imap_dfr(aggression_test_variables, function(group, term) {
  test <- safe_fisher_test(group, aggression_individual$agg_bite_present)
  tibble::tibble(analysis_family = "Focal behavior", response = "Aggression", outcome = "consensus bite presence", 
    term, test_method = test$method, statistic = NA_real_, df = NA_real_, p_value = test$p.value, model = "Separate Monte Carlo Fisher test; one row per individual", 
    note = paste("Sparse outcome; no multivariable aggression GLMM or", "model-derived marginal means"))
})

# Estimate population means from each focal model
boldness_emm <- emmeans::emmeans(boldness_model, ~pop_id, weights = "proportional")
neophilia_emm <- emmeans::emmeans(neophilia_model, ~pop_id, weights = "proportional")
exploration_emm <- emmeans::emmeans(exploration_model, ~pop_id, weights = "proportional")

# Return focal estimates on original response scales
focal_estimates <- dplyr::bind_rows(estimated_mean_table(boldness_emm, "Boldness", "estimated emergence latency in seconds", 
  transform = expm1), estimated_mean_table(emmeans::emmeans(neophilia_model, ~pop_id, type = "response", weights = "proportional"), 
  "Neophilia", "expected recorder-consensus contacts"), estimated_mean_table(emmeans::emmeans(exploration_model, 
  ~pop_id, type = "response", weights = "proportional"), "Exploration", "expected unique corners out of three", 
  multiplier = 3))

# Calculate observed aggression probabilities and exact binomial intervals
aggression_estimates <- aggression_individual %>% dplyr::group_by(pop_id) %>% dplyr::summarise(n = dplyr::n(), 
  events = sum(agg_bite_present), .groups = "drop") %>% dplyr::mutate(confidence_interval = purrr::map2(events, 
  n, ~stats::binom.test(.x, .y)$conf.int), analysis_family = "Focal behavior", response = "Aggression", pop_id = as.character(pop_id), 
  behavioral_assay = NA_character_, estimate = events/n, lower_95 = purrr::map_dbl(confidence_interval, 1), upper_95 = purrr::map_dbl(confidence_interval, 
    2), scale = "observed probability of any bite", source = "Consensus proportion with exact binomial confidence interval") %>% 
  dplyr::select(analysis_family, response, pop_id, behavioral_assay, estimate, lower_95, upper_95, scale, source, 
    n)

# Add observed aggression probabilities to focal estimate table
focal_estimates <- dplyr::bind_rows(focal_estimates, aggression_estimates)

# Start empty post hoc table and add contrasts only after supported tests
posthoc_contrasts <- empty_contrast_table()

# Compare boldness populations only after supported population test
if (isTRUE(term_p_value(focal_model_tests, "Boldness", "pop_id") < 0.05)) {
  posthoc_contrasts <- dplyr::bind_rows(posthoc_contrasts, pairwise_table(boldness_emm, "Focal behavior", "Boldness", 
    "pop_id", "ratio of latency plus one", transform = exp))
}

# Compare neophilia populations only after supported population test
if (isTRUE(term_p_value(focal_model_tests, "Neophilia", "pop_id") < 0.05)) {
  posthoc_contrasts <- dplyr::bind_rows(posthoc_contrasts, pairwise_table(neophilia_emm, "Focal behavior", "Neophilia", 
    "pop_id", "rate ratio", transform = exp))
}

# Compare exploration populations only after supported population test
if (isTRUE(term_p_value(focal_model_tests, "Exploration", "pop_id") < 0.05)) {
  posthoc_contrasts <- dplyr::bind_rows(posthoc_contrasts, pairwise_table(exploration_emm, "Focal behavior", 
    "Exploration", "pop_id", "odds ratio", transform = exp))
}

# Compare aggression populations only after supported global Fisher test
if (isTRUE(term_p_value(aggression_tests, "Aggression", "pop_id") < 0.05)) {
  aggression_pairs <- utils::combn(levels(droplevels(aggression_individual$pop_id)), 2, simplify = FALSE)
  aggression_pairwise <- purrr::map_dfr(aggression_pairs, function(pair) {
    pair_data <- aggression_individual %>% dplyr::filter(pop_id %in% pair) %>% dplyr::mutate(pop_id = droplevels(pop_id))
    test <- stats::fisher.test(table(pair_data$pop_id, factor(pair_data$agg_bite_present, levels = c(0, 1))))
    estimate <- if (length(test$estimate) == 0) {
      NA_real_
    } else {
      unname(test$estimate)
    }
    interval <- if (is.null(test$conf.int)) {
      c(NA_real_, NA_real_)
    } else {
      test$conf.int
    }
    tibble::tibble(analysis_family = "Focal behavior", response = "Aggression", effect = "pop_id", conditioning = NA_character_, 
      contrast = paste(pair, collapse = " / "), estimate, lower_95 = interval[[1]], upper_95 = interval[[2]], 
      scale = "conditional odds ratio", statistic = NA_real_, df = NA_real_, p_value = test$p.value, note = "Exact pairwise Fisher test after the supported global test")
  }) %>% dplyr::mutate(p_holm = stats::p.adjust(p_value, method = "holm")) %>% dplyr::select(dplyr::all_of(names(empty_contrast_table())))
  posthoc_contrasts <- dplyr::bind_rows(posthoc_contrasts, aggression_pairwise)
}

# =============================================================================
# Question 3: Repeated climbing and sheltering
# =============================================================================

# Convert Type-III repeated-behavior tests to common result structure
repeated_anova_rows <- function(model, response) {
  test_table <- car::Anova(model, type = "III") %>% as.data.frame() %>% tibble::rownames_to_column("term") %>% 
    tibble::as_tibble() %>% dplyr::filter(term != "(Intercept)")
  statistic_column <- pick_column(test_table, "Chisq", "test statistic")
  df_column <- pick_column(test_table, "Df", "degrees of freedom")
  p_column <- pick_column(test_table, "Pr(>Chisq)", "P value")
  tibble::tibble(analysis_family = "Repeated behavior", response, outcome = "logit-transformed proportion of observation time", 
    term = test_table$term, test_method = "Type-III Wald chi-square with sum contrasts", statistic = as.numeric(test_table[[statistic_column]]), 
    df = as.numeric(test_table[[df_column]]), p_value = as.numeric(test_table[[p_column]]), model = formula_text(model), 
    note = ifelse(test_table$term == "behavioral_assay:pop_id", "Prespecified test of whether assay patterns differ among populations", 
      NA_character_))
}

# Fit assay-by-population mixed model for each repeated behavior
fit_repeated_behavior <- function(response_variable, response_label) {
  analysis_data <- repeated_data %>% dplyr::filter(dplyr::if_all(dplyr::all_of(c(response_variable, "behavioral_assay", 
    "pop_id", "sex_ternary", "svl_z", "id")), ~!is.na(.x))) %>% droplevels()
  model_formula <- stats::as.formula(paste(response_variable, "~ behavioral_assay * pop_id + sex_ternary + svl_z + (1 | id)"))
  model <- lme4::lmer(model_formula, data = analysis_data, REML = TRUE, control = lmer_control)
  tests <- repeated_anova_rows(model, response_label)
  interaction_p <- term_p_value(tests, response_label, "behavioral_assay:pop_id")
  list(response = response_label, response_variable = response_variable, data = analysis_data, model = model, 
    tests = tests, interaction_supported = isTRUE(interaction_p < 0.05))
}

# Fit prespecified models for climbing and sheltering
climbing_fit <- fit_repeated_behavior("logit_climbing", "Climbing")
sheltering_fit <- fit_repeated_behavior("logit_sheltering", "Sheltering")

# Combine climbing and sheltering fixed-effect tests
repeated_model_tests <- dplyr::bind_rows(climbing_fit$tests, sheltering_fit$tests)

# Back-transform adjusted assay-by-population means to proportions
repeated_estimates_for <- function(fit) {
  marginal_means <- emmeans::emmeans(fit$model, ~behavioral_assay * pop_id, weights = "proportional")
  table <- as.data.frame(stats::confint(marginal_means, level = 0.95))
  tibble::tibble(analysis_family = "Repeated behavior", response = fit$response, pop_id = as.character(table$pop_id), 
    behavioral_assay = as.character(table$behavioral_assay), estimate = stats::plogis(table$emmean), lower_95 = stats::plogis(table$lower.CL), 
    upper_95 = stats::plogis(table$upper.CL), scale = "adjusted proportion on the inverse-logit scale", source = "Prespecified assay-by-population mixed model", 
    n = NA_integer_)
}

# Calculate response-scale estimates for both repeated behaviors
repeated_estimates <- dplyr::bind_rows(repeated_estimates_for(climbing_fit), repeated_estimates_for(sheltering_fit))

# Select post hoc contrasts based on supported fixed-effect structure
repeated_contrasts_for <- function(fit) {
  if (fit$interaction_supported) {
    marginal_means <- emmeans::emmeans(fit$model, ~behavioral_assay | pop_id, weights = "proportional")
    return(pairwise_table(marginal_means, "Repeated behavior", fit$response, "behavioral_assay within population", 
      "odds ratio of observed proportion", conditioning_column = "pop_id", transform = exp, note = paste("Generated because the assay-by-population interaction", 
        "was supported")))
  }
  output <- empty_contrast_table()
  assay_p <- term_p_value(fit$tests, fit$response, "behavioral_assay")
  population_p <- term_p_value(fit$tests, fit$response, "pop_id")
  if (isTRUE(assay_p < 0.05)) {
    output <- dplyr::bind_rows(output, pairwise_table(emmeans::emmeans(fit$model, ~behavioral_assay, weights = "proportional"), 
      "Repeated behavior", fit$response, "behavioral_assay", "odds ratio of observed proportion", transform = exp, 
      note = "Generated after a supported overall assay effect"))
  }
  if (isTRUE(population_p < 0.05)) {
    output <- dplyr::bind_rows(output, pairwise_table(emmeans::emmeans(fit$model, ~pop_id, weights = "proportional"), 
      "Repeated behavior", fit$response, "pop_id", "odds ratio of observed proportion", transform = exp, 
      note = "Generated after a supported overall population effect"))
  }
  output
}

# Add supported repeated-behavior contrasts
posthoc_contrasts <- dplyr::bind_rows(posthoc_contrasts, repeated_contrasts_for(climbing_fit), repeated_contrasts_for(sheltering_fit))

# Extract first finite scalar from version-dependent rptR outputs
first_finite <- function(x) {
  values <- suppressWarnings(as.numeric(unlist(x)))
  values <- values[is.finite(values)]
  if (length(values) == 0) 
    NA_real_
  else values[[1]]
}

# Calculate individual and residual variance from Gaussian mixed model
repeatability_variances <- function(model) {
  variance <- as.data.frame(lme4::VarCorr(model))
  id_variance <- variance %>% dplyr::filter(grp == "id") %>% dplyr::pull(vcov)
  residual_variance <- variance %>% dplyr::filter(grp == "Residual") %>% dplyr::pull(vcov)
  tibble::tibble(var_id = id_variance[[1]], var_residual = residual_variance[[1]], R_from_model = id_variance[[1]]/(id_variance[[1]] + 
    residual_variance[[1]]))
}

# Estimate global adjusted repeatability with 1,000 bootstrap and permutation draws
repeatability_row <- function(fit, resamples = 1000L) {
  if (lme4::isSingular(fit$model, tol = 1e-05)) {
    stop(fit$response, " primary repeatability model is singular and should not be interpreted.", call. = FALSE)
  }
  result <- suppressWarnings(suppressMessages(rptR::rpt(stats::formula(fit$model), grname = "id", data = fit$data, 
    datatype = "Gaussian", nboot = resamples, npermut = resamples)))
  confidence_interval <- suppressWarnings(as.numeric(unlist(result$CI_emp)))
  confidence_interval <- confidence_interval[is.finite(confidence_interval)]
  repeatability_estimate <- first_finite(result$R)
  if (!is.finite(repeatability_estimate)) {
    stop(fit$response, " repeatability estimate was not returned by rptR.", call. = FALSE)
  }
  if (length(confidence_interval) < 2) {
    stop(fit$response, " repeatability bootstrap did not return a two-sided interval.", call. = FALSE)
  }
  p_values_named <- suppressWarnings(unlist(result$P))
  p_values <- suppressWarnings(as.numeric(p_values_named))
  finite_p <- is.finite(p_values)
  p_values <- p_values[finite_p]
  p_names <- names(p_values_named)[finite_p]
  permutation_index <- grep("permut", p_names, ignore.case = TRUE)
  variance <- repeatability_variances(fit$model)
  id_counts <- table(fit$data$id)
  tibble::tibble(response = fit$response, adjusted_for = paste("assay * population + sex/reproductive status + standardized SVL"), 
    n_observations = nrow(fit$data), n_individuals = dplyr::n_distinct(fit$data$id), n_individuals_repeated = sum(id_counts >= 
      2), R = repeatability_estimate, SE = first_finite(result$se), CI_lower = confidence_interval[[1]], 
    CI_upper = confidence_interval[[2]], P_permutation = if (length(permutation_index) > 0) {
      p_values[[permutation_index[[1]]]]
    }
    else if (length(p_values) >= 2) {
      p_values[[2]]
    } else {
      first_finite(result$P)
    }, var_id = variance$var_id, var_residual = variance$var_residual, R_from_model = variance$R_from_model, 
    n_bootstrap = resamples, n_permutations = resamples)
}

# Estimate global repeatability for climbing and sheltering
repeatability_results <- dplyr::bind_rows(repeatability_row(climbing_fit), repeatability_row(sheltering_fit))

# =============================================================================
# Consolidated manuscript results
# =============================================================================

# Combine focal and repeated-behavior fixed-effect tests
model_tests <- dplyr::bind_rows(focal_model_tests, aggression_tests, repeated_model_tests)

# Combine focal and repeated-behavior estimates
adjusted_estimates <- dplyr::bind_rows(focal_estimates, repeated_estimates)

# Format model tests for primary results table
primary_model_rows <- model_tests %>% dplyr::transmute(result_group = "Model test", analysis = analysis_family, 
  response, comparison = term, subgroup = NA_character_, estimate = NA_real_, standard_error = NA_real_, lower_95 = NA_real_, 
  upper_95 = NA_real_, statistic, df, p_value, p_adjusted = NA_real_, n = NA_real_, scale = outcome, method = test_method, 
  model, note = ifelse(term == "svl_z", ifelse(is.na(note), svl_standardization_note, paste(note, svl_standardization_note)), 
    note))

# Format overall and population correlations for primary results table
primary_correlation_rows <- correlation_results %>% dplyr::transmute(result_group = "Behavioral correlation", analysis = scope, 
  response = paste(Var1, "vs", Var2), comparison = "Spearman rank correlation", subgroup = Group, estimate = Correlation, 
  standard_error = NA_real_, lower_95 = NA_real_, upper_95 = NA_real_, statistic = NA_real_, df = NA_real_, p_value = P_Value, 
  p_adjusted = P_Holm, n = as.numeric(N), scale = "Spearman rho", method = "Pairwise-complete Spearman correlation", 
  model = NA_character_, note = direction_note)

# Format Fisher z population-comparison tests for primary results table
primary_fisher_rows <- correlation_population_tests %>% dplyr::transmute(result_group = "Correlation population comparison", 
  analysis = "Fisher r-to-z homogeneity", response = paste(Var1, "vs", Var2), comparison = "Population difference in correlation", 
  subgroup = populations_used, estimate = NA_real_, standard_error = NA_real_, lower_95 = NA_real_, upper_95 = NA_real_, 
  statistic = fisher_Q, df = as.numeric(fisher_df), p_value = fisher_p, p_adjusted = fisher_p_holm, n = as.numeric(k_populations_used), 
  scale = "Fisher z heterogeneity", method = "Approximate Fisher r-to-z homogeneity test", model = NA_character_, 
  note = fisher_note)

# Format rank-interaction sensitivity tests for primary results table
primary_rank_rows <- correlation_population_tests %>% dplyr::transmute(result_group = "Correlation population comparison", 
  analysis = "Rank interaction sensitivity test", response = paste(Var1, "vs", Var2), comparison = "Behavior-by-population interaction", 
  subgroup = NA_character_, estimate = NA_real_, standard_error = NA_real_, lower_95 = NA_real_, upper_95 = NA_real_, 
  statistic = rank_interaction_F, df = as.numeric(rank_interaction_df), p_value = rank_interaction_p, p_adjusted = rank_interaction_p_holm, 
  n = as.numeric(rank_model_N), scale = "Rank-transformed F test", method = "Rank-transformed linear interaction model", 
  model = "rank(y) ~ rank(x) * population", note = rank_interaction_note)

# Format global climbing and sheltering repeatability for primary results table
primary_repeatability_rows <- repeatability_results %>% dplyr::transmute(result_group = "Repeatability", analysis = "Global adjusted repeatability", 
  response, comparison = "Individual random-intercept repeatability", subgroup = "All populations and sexes", 
  estimate = R, standard_error = SE, lower_95 = CI_lower, upper_95 = CI_upper, statistic = NA_real_, df = NA_real_, 
  p_value = P_permutation, p_adjusted = NA_real_, n = as.numeric(n_observations), scale = "Adjusted repeatability R", 
  method = "rptR bootstrap interval and permutation test", model = adjusted_for, note = paste0(n_individuals, 
    " individuals; ", n_individuals_repeated, " measured in at least two assays; ", n_bootstrap, " bootstrap and permutation draws. Singular resample fits are ", 
    "boundary estimates; the primary fitted model must be nonsingular."))

# Combine analysis results in primary results table
primary_results <- dplyr::bind_rows(primary_model_rows, primary_correlation_rows, primary_fisher_rows, primary_rank_rows, 
  primary_repeatability_rows)

# Format adjusted estimates for estimates-and-contrasts table
estimate_rows <- adjusted_estimates %>% dplyr::transmute(result_group = "Adjusted estimate", analysis = analysis_family, 
  response, comparison = ifelse(is.na(behavioral_assay), paste0("pop_id=", pop_id), paste0("pop_id=", pop_id, 
    "; assay=", behavioral_assay)), subgroup = as.character(pop_id), estimate, standard_error = NA_real_, lower_95, 
  upper_95, statistic = NA_real_, df = NA_real_, p_value = NA_real_, p_adjusted = NA_real_, n = as.numeric(n), 
  scale, method = source, model = NA_character_, note = NA_character_)

# Format supported post hoc contrasts for estimates-and-contrasts table
contrast_rows <- posthoc_contrasts %>% dplyr::transmute(result_group = "Post hoc contrast", analysis = analysis_family, 
  response, comparison = contrast, subgroup = conditioning, estimate, standard_error = NA_real_, lower_95, upper_95, 
  statistic, df, p_value, p_adjusted = p_holm, n = NA_real_, scale, method = paste("Pairwise contrast for", effect), 
  model = NA_character_, note)

# Combine estimates and contrasts in second results table
estimates_and_contrasts <- dplyr::bind_rows(estimate_rows, contrast_rows)

# The broader audit tables are retained above, but they are not manuscript tables.
# readr::write_csv(primary_results, "results/Manuscript_primary_results.csv")
# readr::write_csv(estimates_and_contrasts, "results/Manuscript_estimates_and_contrasts.csv")
# print(primary_results, n = Inf, width = Inf)
# print(estimates_and_contrasts, n = Inf, width = Inf)

# =============================================================================
# Tables 1-3 in the manuscript
# =============================================================================

# Table 1 reports the complete capture sample, including animals without a usable
# behavioral record; its population-history fields therefore come from the study
# sampling metadata rather than being reconstructed from the analysis datasets.
manuscript_table_1 <- tibble::tribble(
  ~Population, ~Geography, ~Status, ~`Introduction or colonization history`, ~`Sample size total (male, female)`,
  "Calci", "Mainland Tuscany", "Native", "Native-range population", "48 (28, 20)",
  "Livorno", "Mainland Tuscany", "Native", "Native-range population", "31 (18, 13)",
  "Capraia", "Island", "Biological dispersal", "Established through biological dispersal", "52 (29, 23)",
  "Elba", "Island", "Transient, complex history", "Land-bridge history with local extirpation or reintroduction", "26 (13, 13)",
  "Gorgona", "Island", "Recent introduction", "Recently introduced, reported within last 20 years", "26 (17, 9)"
)

# Format displayed P values exactly as in the manuscript tables.
format_table_p <- function(value) {
  ifelse(is.na(value), "—", ifelse(value < 0.001, "<0.001", sprintf("%.3f", value)))
}

# Assemble Table 2 from the overall tests and the ranges across five populations.
behavior_pair_labels <- tibble::tribble(
  ~Var1, ~Var2, ~pair_order, ~`Behavioral comparison`,
  "agg_bite_freq", "boldness_score", 1, "Aggression × boldness",
  "agg_bite_freq", "neo_freq", 2, "Aggression × neophilia",
  "agg_bite_freq", "exp_unique_corners", 3, "Aggression × exploration",
  "boldness_score", "neo_freq", 4, "Boldness × neophilia",
  "boldness_score", "exp_unique_corners", 5, "Boldness × exploration",
  "neo_freq", "exp_unique_corners", 6, "Neophilia × exploration"
)

within_population_ranges <- correlations_by_population %>% dplyr::group_by(Var1, Var2) %>% dplyr::summarise(
  rho_min = min(Correlation, na.rm = TRUE), rho_max = max(Correlation, na.rm = TRUE),
  p_min = min(P_Value, na.rm = TRUE), p_max = max(P_Value, na.rm = TRUE),
  p_holm_min = min(P_Holm, na.rm = TRUE), p_holm_max = max(P_Holm, na.rm = TRUE), .groups = "drop"
)

manuscript_table_2 <- correlations_overall %>% dplyr::left_join(within_population_ranges, by = c("Var1", "Var2")) %>%
  dplyr::left_join(behavior_pair_labels, by = c("Var1", "Var2")) %>% dplyr::arrange(pair_order) %>%
  dplyr::transmute(
    `Behavioral comparison`, `Overall ρ` = sprintf("%.3f", Correlation),
    `Overall P (Holm P)` = paste0(format_table_p(P_Value), " (", format_table_p(P_Holm), ")"),
    `Overall N` = N, `Within-population ρ range` = paste0(sprintf("%.3f", rho_min), " to ", sprintf("%.3f", rho_max)),
    `Within-population raw P range (Holm P range)` = paste0(format_table_p(p_min), " - ", format_table_p(p_max),
      " (", format_table_p(p_holm_min), " - ", format_table_p(p_holm_max), ")")
  )

# Define the displayed row order and labels for Table 3, then insert the fitted results.
table_3_layout <- tibble::tribble(
  ~response, ~term, ~row_order, ~Response, ~Model, ~Predictor,
  "Aggression", "pop_id", 1, "Aggression: bite presence", "Monte Carlo Fisher exact test", "Population",
  "Aggression", "sex_ternary", 2, "Aggression: bite presence", "Monte Carlo Fisher exact test", "Sex/reproductive status",
  "Aggression", "nba_order_agg", 3, "Aggression: bite presence", "Monte Carlo Fisher exact test", "NBA assay order",
  "Boldness", "pop_id", 4, "Boldness: log (1 + latency)", "LMM", "Population",
  "Boldness", "sex_ternary", 5, "Boldness: log (1 + latency)", "LMM", "Sex/reproductive status",
  "Boldness", "svl_z", 6, "Boldness: log (1 + latency)", "LMM", "Standardized SVL",
  "Boldness", "nba_order_bold", 7, "Boldness: log (1 + latency)", "LMM", "NBA assay order",
  "Neophilia", "pop_id", 8, "Neophilia: contact (consensus)", "Quasi-Poisson GLM", "Population",
  "Neophilia", "sex_ternary", 9, "Neophilia: contact (consensus)", "Quasi-Poisson GLM", "Sex/reproductive status",
  "Neophilia", "svl_z", 10, "Neophilia: contact (consensus)", "Quasi-Poisson GLM", "Standardized SVL",
  "Neophilia", "nba_order_neo", 11, "Neophilia: contact (consensus)", "Quasi-Poisson GLM", "NBA assay order",
  "Exploration", "pop_id", 12, "Exploration: unique corners", "Binomial GLMM", "Population",
  "Exploration", "sex_ternary", 13, "Exploration: unique corners", "Binomial GLMM", "Sex/reproductive status",
  "Exploration", "svl_z", 14, "Exploration: unique corners", "Binomial GLMM", "Standardized SVL",
  "Exploration", "nba_day_exp", 15, "Exploration: unique corners", "Binomial GLMM", "Exploration day",
  "Climbing", "behavioral_assay", 16, "Climbing proportion", "LMM, logit proportion", "Behavioral assay",
  "Climbing", "pop_id", 17, "Climbing proportion", "LMM, logit proportion", "Population",
  "Climbing", "sex_ternary", 18, "Climbing proportion", "LMM, logit proportion", "Sex/reproductive status",
  "Climbing", "svl_z", 19, "Climbing proportion", "LMM, logit proportion", "Standardized SVL",
  "Climbing", "behavioral_assay:pop_id", 20, "Climbing proportion", "LMM, logit proportion", "Assay × population",
  "Sheltering", "behavioral_assay", 21, "Sheltering proportion", "LMM, logit proportion", "Behavioral assay",
  "Sheltering", "pop_id", 22, "Sheltering proportion", "LMM, logit proportion", "Population",
  "Sheltering", "sex_ternary", 23, "Sheltering proportion", "LMM, logit proportion", "Sex/reproductive status",
  "Sheltering", "svl_z", 24, "Sheltering proportion", "LMM, logit proportion", "Standardized SVL",
  "Sheltering", "behavioral_assay:pop_id", 25, "Sheltering proportion", "LMM, logit proportion", "Assay × population"
)

neophilia_residual_df <- focal_model_tests %>% dplyr::filter(response == "Neophilia", term == "Residuals") %>%
  dplyr::pull(df) %>% dplyr::first()

manuscript_table_3 <- table_3_layout %>% dplyr::left_join(
  model_tests %>% dplyr::select(response, term, statistic, df, p_value), by = c("response", "term")
) %>% dplyr::arrange(row_order) %>% dplyr::transmute(
  Response, Model, Predictor,
  `Test statistic` = dplyr::case_when(
    is.na(statistic) ~ "—", response == "Neophilia" ~ paste0("F = ", sprintf("%.2f", statistic)),
    TRUE ~ paste0("χ² = ", sprintf("%.2f", statistic))
  ),
  df = dplyr::case_when(
    is.na(df) ~ "—", response == "Neophilia" ~ paste0(as.integer(df), ", ", as.integer(neophilia_residual_df)),
    TRUE ~ as.character(as.integer(df))
  ),
  p = format_table_p(p_value)
)

# Save and display only Tables 1-3 from the manuscript.
readr::write_csv(manuscript_table_1, "results/Table1_population_history_and_sample_composition.csv")
readr::write_csv(manuscript_table_2, "results/Table2_behavioral_correlations.csv")
readr::write_csv(manuscript_table_3, "results/Table3_fixed_effect_tests.csv")
print(manuscript_table_1, n = Inf, width = Inf)
print(manuscript_table_2, n = Inf, width = Inf)
print(manuscript_table_3, n = Inf, width = Inf)

# =============================================================================
# Descriptive summaries used only for figures
# =============================================================================

# Return finite mean while preserving missing groups
mean_or_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) 
    NA_real_
  else mean(x)
}

# Return SD when at least two values are available
sd_or_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2) 
    NA_real_
  else stats::sd(x)
}

# Summarize focal responses by population for bar and violin figures
summarize_population_behavior <- function(data, response, behavior) {
  data %>% dplyr::group_by(pop_id) %>% dplyr::summarise(n = sum(!is.na(.data[[response]])), mean = mean_or_na(.data[[response]]), 
    sd = sd_or_na(.data[[response]]), se = sd/sqrt(n), .groups = "drop") %>% dplyr::mutate(behavior, response, 
    .before = 1)
}

# Combine focal population summaries without extra files
population_behavior_summary <- dplyr::bind_rows(summarize_population_behavior(individual_data, "agg_bite_freq", 
  "Aggression"), summarize_population_behavior(individual_data, "bold_latency", "Boldness"), summarize_population_behavior(individual_data, 
  "neo_freq", "Neophilia"), summarize_population_behavior(individual_data, "exp_unique_corners", "Exploration"))

# Summarize repeated behaviors by population and assay
summarize_repeated_behavior <- function(data, response, behavior) {
  data %>% dplyr::group_by(pop_id, behavioral_assay) %>% dplyr::summarise(n = sum(!is.na(.data[[response]])), 
    mean_prop = mean_or_na(.data[[response]]), sd_prop = sd_or_na(.data[[response]]), se_prop = sd_prop/sqrt(n), 
    .groups = "drop") %>% dplyr::mutate(behavior, response, summary_scope = "population_by_assay", .before = 1)
}

# Combine climbing and sheltering summaries for assay figure
repeated_behavior_summary <- dplyr::bind_rows(summarize_repeated_behavior(repeated_data, "climbing_prop", "Climbing"), 
  summarize_repeated_behavior(repeated_data, "shelter_prop", "Sheltering"))

# Calculate focal means and SEs by population and sex
behavior_sex_plot_data <- individual_data %>% dplyr::mutate(agg_bite_present = dplyr::if_else(is.na(agg_bite_freq), 
  NA_real_, as.numeric(agg_bite_freq > 0))) %>% dplyr::select(pop_id, sex_binary, agg_bite_present, bold_latency, 
  neo_freq, exp_unique_corners) %>% tidyr::pivot_longer(cols = c(agg_bite_present, bold_latency, neo_freq, exp_unique_corners), 
  names_to = "response", values_to = "value") %>% dplyr::filter(!is.na(sex_binary), !is.na(value)) %>% dplyr::group_by(pop_id, 
  sex_binary, response) %>% dplyr::summarise(n = dplyr::n(), estimate = mean(value), se = if (dplyr::n() > 1) {
  stats::sd(value)/sqrt(dplyr::n())
} else {
  NA_real_
}, .groups = "drop") %>% dplyr::mutate(lower = pmax(estimate - se, 0), upper = dplyr::case_when(response == "agg_bite_present" ~ 
  pmin(estimate + se, 1), response == "exp_unique_corners" ~ pmin(estimate + se, 3), TRUE ~ estimate + se), behavior = dplyr::recode(response, 
  agg_bite_present = "Aggression\n(probability of biting)", bold_latency = "Boldness\n(emergence latency, s)", 
  neo_freq = "Neophilia\n(novel-object contacts)", exp_unique_corners = "Exploration\n(unique corners visited)"))

# =============================================================================
# Manuscript figure settings
# =============================================================================

# Set geographic population order across figures
figure_population_order <- c("gor", "elb", "cap", "cal", "liv")

# Translate population codes to figure labels
population_labels <- c(gor = "Gorgona", elb = "Elba", cap = "Capraia", cal = "Calci", liv = "Livorno")

# Apply established population palette
population_colors <- c(gor = "#7A263A", elb = "#D97A32", cap = "#7A7D3B", cal = "#3C7D8C", liv = "#4C67A1")

# Reorder population factors after model fitting
individual_data <- individual_data %>% dplyr::mutate(pop_id = factor(pop_id, levels = figure_population_order))
repeated_data <- repeated_data %>% dplyr::mutate(pop_id = factor(pop_id, levels = figure_population_order), behavioral_assay = factor(behavioral_assay, 
  levels = assay_order))
population_behavior_summary <- population_behavior_summary %>% dplyr::mutate(pop_id = factor(pop_id, levels = figure_population_order))
repeated_behavior_summary <- repeated_behavior_summary %>% dplyr::mutate(pop_id = factor(pop_id, levels = figure_population_order), 
  behavioral_assay = factor(behavioral_assay, levels = assay_order))
behavior_sex_plot_data <- behavior_sex_plot_data %>% dplyr::mutate(pop_id = factor(pop_id, levels = figure_population_order), 
  sex_binary = factor(sex_binary, levels = c("Female", "Male")), behavior = factor(behavior, levels = c("Aggression\n(probability of biting)", 
    "Boldness\n(emergence latency, s)", "Neophilia\n(novel-object contacts)", "Exploration\n(unique corners visited)")))
repeatability_results <- repeatability_results %>% dplyr::mutate(response = factor(response, levels = c("Climbing", 
  "Sheltering")))

# Apply consistent typography, panel tags, and legend placement
manuscript_theme <- function(show_legend = FALSE) {
  ggplot2::theme_classic(base_size = 12) + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 16, 
    hjust = 0.5), axis.title = ggplot2::element_text(face = "bold", size = 13), axis.text = ggplot2::element_text(size = 11, 
    color = "black"), plot.tag = ggplot2::element_text(face = "bold", size = 18), plot.tag.position = "topleft", 
    legend.position = if (show_legend) 
      "right"
    else "none")
}

# Save each figure as high-resolution PNG and editable PDF
save_figure <- function(plot, filename, width, height, dpi = 600) {
  ggplot2::ggsave(file.path("figures", paste0(filename, ".png")), plot, width = width, height = height, dpi = dpi)
  ggplot2::ggsave(file.path("figures", paste0(filename, ".pdf")), plot, width = width, height = height)
}

# Draw focal population mean and SE panel
plot_population_bar <- function(behavior, y_label, title, tag) {
  plot_data <- population_behavior_summary %>% dplyr::filter(.data$behavior == .env$behavior)
  ggplot2::ggplot(plot_data, ggplot2::aes(x = pop_id, y = mean, fill = pop_id)) + ggplot2::geom_col(color = "black", 
    alpha = 0.85, na.rm = TRUE) + ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - se, ymax = mean + se), 
    width = 0.2, linewidth = 0.7, na.rm = TRUE) + ggplot2::scale_fill_manual(values = population_colors, labels = population_labels) + 
    ggplot2::scale_x_discrete(labels = population_labels) + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 
    0.1))) + ggplot2::labs(title = title, x = NULL, y = y_label, tag = tag) + manuscript_theme()
}

# Draw observed distribution with individual values and mean marker
plot_violin_distribution <- function(raw_data, summary_data, response, y_label, title, tag, y_limits = NULL) {
  plot <- ggplot2::ggplot(raw_data, ggplot2::aes(x = pop_id, y = .data[[response]], fill = pop_id, color = pop_id)) + 
    ggplot2::geom_violin(alpha = 0.25, trim = TRUE, color = "black", na.rm = TRUE) + ggplot2::geom_jitter(width = 0.16, 
    alpha = 0.55, size = 1.7, na.rm = TRUE) + ggplot2::geom_point(data = summary_data, ggplot2::aes(x = pop_id, 
    y = mean_value), inherit.aes = FALSE, shape = 23, size = 3.5, fill = "white", color = "black", na.rm = TRUE) + 
    ggplot2::scale_fill_manual(values = population_colors, labels = population_labels) + ggplot2::scale_color_manual(values = population_colors, 
    labels = population_labels) + ggplot2::scale_x_discrete(labels = population_labels) + ggplot2::labs(title = title, 
    x = NULL, y = y_label, tag = tag) + manuscript_theme()
  if (!is.null(y_limits)) {
    plot <- plot + ggplot2::coord_cartesian(ylim = y_limits)
  }
  plot
}

# Select focal means for white diamonds in violin panels
focal_point_summaries <- population_behavior_summary %>% dplyr::transmute(behavior, response, pop_id, mean_value = mean)

# Calculate population means across assays for repeated-behavior violins
repeated_point_summaries <- dplyr::bind_rows(repeated_data %>% dplyr::group_by(pop_id) %>% dplyr::summarise(mean_value = mean_or_na(climbing_prop), 
  .groups = "drop") %>% dplyr::mutate(behavior = "Climbing", response = "climbing_prop", .before = 1), repeated_data %>% 
  dplyr::group_by(pop_id) %>% dplyr::summarise(mean_value = mean_or_na(shelter_prop), .groups = "drop") %>% dplyr::mutate(behavior = "Sheltering", 
  response = "shelter_prop", .before = 1))

# =============================================================================
# Figure Focal population means
# =============================================================================

# Build four focal-behavior population panels
population_bar_panels <- list(plot_population_bar("Neophilia", "Mean touches", "Neophilia", "A"), plot_population_bar("Boldness", 
  "Mean latency (s)", "Boldness", "B"), plot_population_bar("Aggression", "Mean bites", "Aggression", "C"), plot_population_bar("Exploration", 
  "Mean unique corners", "Exploration", "D"))

# Arrange and save four focal population panels
figure_population_bars <- gridExtra::arrangeGrob(grobs = population_bar_panels, ncol = 2)
# save_figure(figure_population_bars, "Fig1_BarGraphs", 10, 9)

# =============================================================================
# Figure Full observed distributions
# =============================================================================

# Build six-panel violin figure
violin_panels <- list(plot_violin_distribution(individual_data, dplyr::filter(focal_point_summaries, behavior == 
  "Neophilia"), "neo_freq", "Touches", "Neophilia", "A"), plot_violin_distribution(individual_data, dplyr::filter(focal_point_summaries, 
  behavior == "Boldness"), "bold_latency", "Latency (s)", "Boldness", "B"), plot_violin_distribution(individual_data, 
  dplyr::filter(focal_point_summaries, behavior == "Aggression"), "agg_bite_freq", "Bites", "Aggression", "C"), 
  plot_violin_distribution(individual_data, dplyr::filter(focal_point_summaries, behavior == "Exploration"), 
    "exp_unique_corners", "Unique corners", "Exploration", "D", c(0, 3)), plot_violin_distribution(repeated_data, 
    dplyr::filter(repeated_point_summaries, behavior == "Sheltering"), "shelter_prop", "Proportion of time", 
    "Sheltering", "E", c(0, 1)), plot_violin_distribution(repeated_data, dplyr::filter(repeated_point_summaries, 
    behavior == "Climbing"), "climbing_prop", "Proportion of time", "Climbing", "F", c(0, 1)))

# Arrange and save six observed distributions
figure_violin <- gridExtra::arrangeGrob(grobs = violin_panels, ncol = 2)
# save_figure(figure_violin, "Fig2_ViolinPlots_Full", 11, 16)

# =============================================================================
# Figure Repeated behaviors across assays
# =============================================================================

# Draw assay means and SEs for one repeated behavior
plot_repeated_behavior <- function(behavior, title, tag, show_legend = FALSE) {
  plot_data <- repeated_behavior_summary %>% dplyr::filter(summary_scope == "population_by_assay", .data$behavior == 
    .env$behavior)
  ggplot2::ggplot(plot_data, ggplot2::aes(x = behavioral_assay, y = mean_prop, color = pop_id)) + ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.65), 
    size = 3) + ggplot2::geom_errorbar(ggplot2::aes(ymin = mean_prop - se_prop, ymax = mean_prop + se_prop), 
    position = ggplot2::position_dodge(width = 0.65), width = 0.15, linewidth = 0.6) + ggplot2::coord_cartesian(ylim = c(0, 
    1)) + ggplot2::scale_color_manual(values = population_colors, labels = population_labels, name = "Population") + 
    ggplot2::labs(title = title, x = "Behavioral assay", y = "Mean proportion of time", tag = tag) + manuscript_theme(show_legend) + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
}

# Arrange repeated-behavior panels without connecting assay points
figure_flexibility <- gridExtra::arrangeGrob(plot_repeated_behavior("Climbing", "Climbing", "A"), plot_repeated_behavior("Sheltering", 
  "Sheltering", "B", show_legend = TRUE), ncol = 2)

# Save assay-specific climbing and sheltering panels
# save_figure(figure_flexibility, "Fig3_Flexibility_Proportions", 14, 6)

# =============================================================================
# Manuscript Figure 3 
# =============================================================================

# Connect this script's existing summaries and palette to the approved code.
repeated_summary <- repeated_behavior_summary
pop_colors <- population_colors
pop_labels <- population_labels

# Draw assay-specific means and standard errors for one repeated behavior.
plot_flex <- function(behavior, title, tag, show_legend = FALSE, point_size = 3, errorbar_width = 0.15, errorbar_linewidth = 0.6, 
  x_label = "Behavioral assay") {
  plot_data <- repeated_summary %>% dplyr::filter(summary_scope == "population_by_assay", .data$behavior == .env$behavior)
  ggplot2::ggplot(plot_data, ggplot2::aes(x = behavioral_assay, y = mean_prop, color = pop_id)) + ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.65), 
    size = point_size) + ggplot2::geom_errorbar(ggplot2::aes(ymin = mean_prop - se_prop, ymax = mean_prop + 
    se_prop), position = ggplot2::position_dodge(width = 0.65), width = errorbar_width, linewidth = errorbar_linewidth) + 
    ggplot2::coord_cartesian(ylim = c(0, 1)) + ggplot2::scale_color_manual(values = pop_colors, labels = pop_labels, 
    name = "Population") + ggplot2::labs(title = title, x = x_label, y = "Mean proportion\nof time", tag = tag) + 
    manuscript_theme(show_legend) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
}

# ---- Combined Panels --------------------------------------

large_combined_theme <- function(show_legend = FALSE) {
  ggplot2::theme_classic(base_size = 32) + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 40, 
    hjust = 0.5), axis.title.x = ggplot2::element_text(face = "bold", size = 34), axis.title.y = ggplot2::element_text(face = "bold", 
    size = 32), axis.text = ggplot2::element_text(size = 30, color = "black"), axis.text.x = ggplot2::element_text(angle = 35, 
    hjust = 1, vjust = 1, size = 30), plot.tag = ggplot2::element_text(face = "bold", size = 42), legend.title = ggplot2::element_text(face = "bold", 
    size = 28), legend.text = ggplot2::element_text(size = 26), plot.margin = ggplot2::margin(18, 24, 18, 24), 
    legend.position = if (show_legend) 
      "right"
    else "none")
}

# Rebuild Figure panels with enlarged text.
bar_panels_combined <- list(plot_population_bar("Neophilia", "Mean touches", "Neophilia", "A") + large_combined_theme(), 
  plot_population_bar("Boldness", "Mean latency (s)", "Boldness", "B") + large_combined_theme(), plot_population_bar("Aggression", 
    "Mean bites", "Aggression", "C") + large_combined_theme(), plot_population_bar("Exploration", "Mean unique corners", 
    "Exploration", "D") + large_combined_theme())

fig_population_bars_combined <- gridExtra::arrangeGrob(grobs = bar_panels_combined, ncol = 2)

# Rebuild Figure panels with the combined-figure labels E and F.
# Enlarge points, error-bar caps, and error-bar strokes threefold only here.
# Remove the x-axis title from both combined repeated-behavior panels.
fig_flexibility_combined <- gridExtra::arrangeGrob(plot_flex("Climbing", "Climbing", "E", point_size = 9, errorbar_width = 0.45, 
  errorbar_linewidth = 1.8, x_label = NULL) + large_combined_theme(), plot_flex("Sheltering", "Sheltering", "F", 
  point_size = 9, errorbar_width = 0.45, errorbar_linewidth = 1.8, x_label = NULL) + large_combined_theme(), 
  ncol = 2)

# Combine enlarged panels
fig_population_and_flexibility <- gridExtra::arrangeGrob(fig_population_bars_combined, fig_flexibility_combined, 
  ncol = 1, heights = c(1.5, 1))

# Save the six-panel statistical figure embedded as Figure 3 in the manuscript.
save_figure(fig_population_and_flexibility, "Figure3_behavioral_responses", 14, 15)

# Display the combined figure.
# print(fig_population_and_flexibility)

# =============================================================================
# Figure Descriptive focal responses by sex
# =============================================================================

# Build descriptive BR6 population-by-sex layout
figure_sex_population <- ggplot2::ggplot(behavior_sex_plot_data, ggplot2::aes(x = pop_id, y = estimate, fill = pop_id, 
  alpha = sex_binary, group = sex_binary)) + ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), 
  width = 0.7, color = "black", linewidth = 0.25) + ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), 
  position = ggplot2::position_dodge(width = 0.8), width = 0.18, linewidth = 0.45) + ggplot2::facet_wrap(~behavior, 
  scales = "free_y", ncol = 2) + ggplot2::scale_fill_manual(values = population_colors, breaks = figure_population_order, 
  labels = population_labels, name = "Population") + ggplot2::scale_alpha_manual(values = c(Female = 0.55, Male = 1), 
  name = "Sex") + ggplot2::scale_x_discrete(labels = population_labels) + ggplot2::labs(x = NULL, y = "Observed mean ± SE") + 
  manuscript_theme(show_legend = TRUE) + ggplot2::theme(strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_text(face = "bold"), 
  axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))

# Save descriptive population-by-sex figure
# save_figure(figure_sex_population, "Figure_sex_population_focal_behaviors", 9, 7)

# =============================================================================
# Figure Global adjusted repeatability
# =============================================================================

# Plot global repeatability estimates and confidence intervals
figure_repeatability <- ggplot2::ggplot(repeatability_results, ggplot2::aes(x = response, y = R, fill = response)) + 
  ggplot2::geom_col(width = 0.65, color = "black", linewidth = 0.25, na.rm = TRUE) + ggplot2::geom_errorbar(ggplot2::aes(ymin = CI_lower, 
  ymax = CI_upper), width = 0.14, linewidth = 0.45, na.rm = TRUE) + ggplot2::scale_fill_manual(values = c(Climbing = "#4C67A1", 
  Sheltering = "#7A263A")) + ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) + ggplot2::labs(x = NULL, 
  y = "Adjusted repeatability (R)") + manuscript_theme(show_legend = FALSE) + ggplot2::theme(axis.text.x = ggplot2::element_text(face = "bold"))

# Save global repeatability figure
# save_figure(figure_repeatability, "Figure_repeatability_global", 6, 5)

# =============================================================================
# Figure Behavioral-correlation heatmap
# =============================================================================

# Translate analysis variables to concise behavior labels
behavior_labels <- c(agg_bite_freq = "Aggression", boldness_score = "Boldness", neo_freq = "Neophilia", exp_unique_corners = "Exploration")

# Prepare correlation labels and significance symbols from existing tests
correlation_plot_data <- correlation_results %>% dplyr::mutate(Group = dplyr::recode(Group, !!!population_labels), 
  Group = factor(Group, levels = c("All individuals", "Calci", "Capraia", "Elba", "Gorgona", "Livorno")), Behavior_1 = factor(behavior_labels[Var1], 
    levels = c("Aggression", "Boldness", "Neophilia", "Exploration")), Behavior_2 = factor(behavior_labels[Var2], 
    levels = c("Aggression", "Boldness", "Neophilia", "Exploration")), significance = dplyr::case_when(!is.na(P_Holm) & 
    P_Holm < 0.05 ~ "*", !is.na(P_Value) & P_Value < 0.05 ~ "†", TRUE ~ ""), cell_label = ifelse(is.na(Correlation), 
    "NE", paste0(sprintf("%.2f", Correlation), significance)))

# Display correlations as heatmap
figure_correlations <- ggplot2::ggplot(correlation_plot_data, ggplot2::aes(Behavior_2, Behavior_1, fill = Correlation)) + 
  ggplot2::geom_tile(color = "white", linewidth = 1) + ggplot2::geom_text(ggplot2::aes(label = cell_label), size = 4.2) + 
  ggplot2::facet_wrap(~Group, ncol = 2) + ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", 
  midpoint = 0, limits = c(-1, 1), name = expression(Spearman ~ rho)) + ggplot2::coord_fixed() + ggplot2::labs(x = NULL, 
  y = NULL) + manuscript_theme(show_legend = TRUE) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, 
  hjust = 1), strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_text(face = "bold"), 
  panel.spacing = grid::unit(1, "lines"))

# Save supplementary behavioral-correlation heatmap
# save_figure(figure_correlations, "Figure_S_behavioral_syndrome_correlations", 9, 11)

# =============================================================================
# Display final outputs
# =============================================================================

# Display figures
# print(figure_population_bars)
# print(figure_violin)
# print(figure_flexibility)
# print(figure_sex_population)
# print(figure_repeatability)
# print(figure_correlations)

# Close the null device used only in noninteractive runs.
if (!interactive()) grDevices::dev.off()
