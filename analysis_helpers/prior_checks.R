# Author: Daniel Franks

# Priors come from each model specification in Wood_ant_analysis.Rmd.
# Seeds, expected-response curves, main plots and summaries are in the notebook.
# Full simulations here retain shared group effects and observation noise.
# Preserve the draw order so fixed seeds reproduce the displayed prior predictions.

# Model 1: Genetic differentiation simulations ----------------------------------
# Simulate FST with all nest, colony, colony-pair and residual variation.

simulate_genetic_prior <- function(model_spec, draws = 400L) {
  prior_settings <- model_spec$prior_parameters
  design_rows <- unique(round(seq(1, nrow(model_spec$data), length.out = min(800, nrow(model_spec$data)))))
  model_data <- model_spec$data[design_rows, , drop = FALSE]
  replicated_outcomes <- matrix(NA_real_, draws, nrow(model_data))
  for (i in seq_len(draws)) {
    # Population relationship, shared endpoints, colony combinations, then residuals.
    logit_fst <- rnorm(1, prior_settings$intercept[1], prior_settings$intercept[2]) +
      rnorm(1, 0, prior_settings$same_colony) * model_data$same_colony +
      rnorm(1, 0, prior_settings$between_slope) * model_data$distance_between +
      rnorm(1, 0, prior_settings$within_slope) * model_data$distance_within
    nest_effects <- rnorm(nlevels(model_data$from_label), 0, abs(rnorm(1, 0, prior_settings$endpoint_sd)))
    logit_fst <- logit_fst + nest_effects[model_data$from_label] + nest_effects[model_data$to_label]
    colony_effects <- rnorm(nlevels(model_data$from_colony), 0, abs(rnorm(1, 0, prior_settings$endpoint_sd)))
    logit_fst <- logit_fst + colony_effects[model_data$from_colony] + colony_effects[model_data$to_colony]
    colony_pair_intercepts <- rnorm(nlevels(model_data$colony_pair), 0, abs(rnorm(1, 0, prior_settings$pair_intercept_sd)))
    logit_fst <- logit_fst + colony_pair_intercepts[model_data$colony_pair]
    colony_pair_slopes <- rnorm(nlevels(model_data$colony_pair), 0, abs(rnorm(1, 0, prior_settings$pair_slope_sd)))
    logit_fst <- logit_fst + colony_pair_slopes[model_data$colony_pair] * model_data$pair_distance_between
    residual_sd <- exp(rnorm(1, prior_settings$log_sigma[1], prior_settings$log_sigma[2]) +
      rnorm(1, 0, prior_settings$sigma_contrast) * model_data$same_colony)
    replicated_outcomes[i, ] <- plogis(rnorm(nrow(model_data), logit_fst, residual_sd))
  }
  list(data = model_data, yrep = replicated_outcomes)
}

# Models 2 and 3: Trail simulations ---------------------------------------------
# Shared predictor design and correlated colony effects; model-specific steps are below.

simulate_prior <- function(model_spec, draws = 400) {
  model_data <- model_spec$data
  # A fixed systematic subset keeps the prior check inexpensive.
  design_rows <- unique(round(seq(1, nrow(model_data), length.out = min(nrow(model_data), 800))))
  model_data <- model_data[design_rows, , drop = FALSE]
  model_data$Intercept <- 1
  design_matrix <- model.matrix(~ 0 + Intercept + distance_z + fst_z + food_z, model_data)
  replicated_outcomes <- matrix(NA_real_, draws, nrow(model_data))
  investment <- model_spec$response_scale == "lognormal"
  prior_settings <- model_spec$prior_parameters
  slope_sd <- prior_settings$slope_sd
  for (r in seq_len(draws)) {
    # Draw population coefficients, then shared nest and colony effects.
    coefficient_draw <- rnorm(ncol(design_matrix), 0, slope_sd)
    if (investment) coefficient_draw[colnames(design_matrix) == "food_z"] <- rnorm(1, 0, prior_settings$food_sd)
    coefficient_draw[colnames(design_matrix) == "fst_z"] <-
      coefficient_draw[colnames(design_matrix) == "fst_z"] * prior_settings$genetic_sd / slope_sd
    if (!investment) coefficient_draw[colnames(design_matrix) == "distance_z"] <-
      -abs(coefficient_draw[colnames(design_matrix) == "distance_z"]) * prior_settings$distance_sd / slope_sd
    coefficient_draw[colnames(design_matrix) == "Intercept"] <- rnorm(1, prior_settings$intercept[1], prior_settings$intercept[2])
    linear_predictor <- drop(design_matrix %*% coefficient_draw)
    nest_sd <- if (investment) abs(rnorm(1, 0, prior_settings$group_sd)) else rexp(1, prior_settings$group_rate)
    nest_effects <- rnorm(nlevels(model_data$from_label), 0, nest_sd)
    linear_predictor <- linear_predictor + nest_effects[model_data$from_label] + nest_effects[model_data$to_label]
    colony_sds <- if (investment)
      abs(rnorm(2, 0, c(prior_settings$group_sd, prior_settings$colony_slope_sd))) else
      rexp(2, prior_settings$group_rate)
    colony_correlation <- 2 * rbeta(1, prior_settings$correlation, prior_settings$correlation) - 1 # two-dimensional LKJ
    standard_normal_draws <- matrix(rnorm(2 * nlevels(model_data$colony)), ncol = 2)
    colony_intercept <- colony_sds[1] * standard_normal_draws[, 1]
    colony_distance <- colony_sds[2] * (colony_correlation * standard_normal_draws[, 1] +
      sqrt(1 - colony_correlation^2) * standard_normal_draws[, 2])
    linear_predictor <- linear_predictor + colony_intercept[model_data$colony]
    if (investment) {
      # Model 3: additive colony distance slopes on log investment.
      linear_predictor <- linear_predictor + colony_distance[model_data$colony] * model_data$distance_z
    } else {
      # Model 2: multiply the shared decline by each colony's positive multiplier.
      linear_predictor <- linear_predictor + coefficient_draw[colnames(design_matrix) == "distance_z"] *
        (exp(colony_distance[model_data$colony]) - 1) * model_data$distance_z
    }
    # Add response variation: binary trails (Model 2) or positive investment (Model 3).
    replicated_outcomes[r, ] <- switch(model_spec$response_scale,
      bernoulli = rbinom(nrow(model_data), 1, plogis(linear_predictor)),
      lognormal = rlnorm(nrow(model_data), linear_predictor, abs(rnorm(1, 0, prior_settings$residual_sd)))
    )
  }
  stopifnot(all(is.finite(replicated_outcomes)))

  list(data = model_data, yrep = replicated_outcomes)
}

# Reusable expected-response calculations --------------------------------------

long_prior_draws <- function(x, predicted_response) {
  stopifnot(all(is.finite(predicted_response)))
  tibble(draw = rep(seq_len(nrow(predicted_response)), times = length(x)),
    x = rep(x, each = nrow(predicted_response)), y = as.vector(predicted_response))
}

# Integrate FST over residual variation, including unusually broad prior draws.
prior_mean_fst <- function(linear_predictor, residual_sd) {
  quadrature <- normal_quadrature(100L)
  predicted_response <- matrix(0, nrow(linear_predictor), ncol(linear_predictor))
  for (k in seq_along(quadrature$nodes)) {
    predicted_response <- predicted_response + quadrature$weights[k] * plogis(linear_predictor + residual_sd * quadrature$nodes[k])
  }
  for (i in which(residual_sd > 2)) for (j in seq_len(ncol(linear_predictor))) {
    predicted_response[i, j] <- integrate(function(standard_normal_draws)
      plogis(linear_predictor[i, j] + residual_sd[i] * standard_normal_draws) * dnorm(standard_normal_draws),
      -Inf, Inf, rel.tol = 1e-8)$value
  }
  predicted_response
}

# Model 3: tail summaries for the genetic-prior robustness check -----------------

genetic_prior_tail_summary <- function(simulations) {
  imap_dfr(simulations, function(prior_simulation, prior) {
    designs <- list(`All trail designs` = seq_len(ncol(prior_simulation$yrep)),
      `Smallest FST` = which.min(prior_simulation$data$fst_z), `Largest FST` = which.max(prior_simulation$data$fst_z),
      `Largest food predictor` = which.max(prior_simulation$data$food_z))
    imap_dfr(designs, function(design_rows, design) tibble(
      prior = prior, design = design,
      `Investment > 1 (%)` = 100 * mean(prior_simulation$yrep[, design_rows] > 1),
      `Investment > 10 (%)` = 100 * mean(prior_simulation$yrep[, design_rows] > 10)
    ))
  })
}
