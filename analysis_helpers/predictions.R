# Author: Daniel Franks

# Reusable expected-response calculations and marginal effects.
# Prediction grids, main numerical results and plotting calls are in Wood_ant_analysis.Rmd.

# Shared numerical helpers ----------------------------------------------------

# normal_quadrature() is also used by prior_checks.R.
# Tables and curves use posterior means and 95% equal-tailed credible intervals.

summarize_columns <- function(draw_matrix) {
  tibble(
    estimate = apply(draw_matrix, 2, mean),
    lower = apply(draw_matrix, 2, quantile, .025), upper = apply(draw_matrix, 2, quantile, .975)
  )
}

# Nodes and weights integrate over a standard normal variable.
normal_quadrature <- function(n = 30) {
  jacobi_matrix <- matrix(0, n, n)
  jacobi_matrix[cbind(1:(n - 1), 2:n)] <- sqrt(1:(n - 1))
  jacobi_matrix <- jacobi_matrix + t(jacobi_matrix)
  eigendecomposition <- eigen(jacobi_matrix, symmetric = TRUE)
  list(nodes = eigendecomposition$values, weights = eigendecomposition$vectors[1, ]^2)
}

# Shared expected-response prediction -----------------------------------------

# Model 1 integrates over Gaussian residual variation before averaging FST.
# Models 2 and 3 use brms expected responses, including the lognormal mean correction.

expected_response <- function(fit, model_spec, newdata, draw_ids, re_formula = NULL) {
  if (model_spec$response_scale != "logit") {
    return(posterior_epred(fit, newdata = newdata, re_formula = re_formula, draw_ids = draw_ids))
  }
  logit_mean_draws <- posterior_linpred(fit, newdata = newdata, re_formula = re_formula, draw_ids = draw_ids)
  logit_sd_draws <- posterior_linpred(fit,
    newdata = newdata, dpar = "sigma", transform = TRUE,
    re_formula = re_formula, draw_ids = draw_ids
  )
  mean_fst(logit_mean_draws, logit_sd_draws)
}

mean_fst <- function(logit_mean_draws, logit_sd_draws) {
  quadrature <- normal_quadrature()
  mean_response <- matrix(0, nrow(logit_mean_draws), ncol(logit_mean_draws))
  for (k in seq_along(quadrature$nodes)) {
    mean_response <- mean_response + quadrature$weights[k] * plogis(logit_mean_draws + logit_sd_draws * quadrature$nodes[k])
  }
  mean_response
}

# Model 1: genetic differentiation calculations -------------------------------

check_prediction_integration <- function() {
  quadrature <- normal_quadrature()
  # Check the expected FST against an independent adaptive integral.
  grid <- expand.grid(mu = c(-4, -2, -.5), sigma = c(.1, .3, .6, 1))
  error <- apply(grid, 1, function(g) {
    numerical <- sum(quadrature$weights * plogis(g[1] + g[2] * quadrature$nodes))
    reference <- integrate(function(z) plogis(g[1] + g[2] * z) * dnorm(z),
      -Inf, Inf, rel.tol = 1e-10)$value
    abs(numerical - reference)
  })
  stopifnot(max(error) < 1e-7)
  invisible(max(error))
}

# Change distance without re-centring the fitted colony-pair slopes or changing endpoints.
shift_genetic_distance <- function(model_data, increase) {
  model_data$distance_m <- model_data$distance_m + increase
  model_data$distance_100m <- (model_data$distance_m - 50) / 100
  model_data$distance_between <- model_data$distance_100m * (1 - model_data$same_colony)
  model_data$distance_within <- model_data$distance_100m * model_data$same_colony
  model_data$pair_distance_between <- (model_data$distance_m - model_data$pair_mean_distance_m) / 100 * (1 - model_data$same_colony)
  model_data
}

# Sum expected FST and local distance slopes in batches over all observed pairs.
# The notebook divides these draw-wise sums by category counts and forms contrasts.
genetic_row_sums <- function(fit, model_spec, derivative_step = .01, draw_ids = seq_len(ndraws(fit)), re_formula = NULL) {
  model_data <- model_spec$data
  stopifnot(length(derivative_step) == 1, is.finite(derivative_step), derivative_step > 0, all(model_data$distance_m > derivative_step))
  group <- ifelse(model_data$same_colony == 1, "Within colony", "Between colonies")
  groups <- c("Within colony", "Between colonies")
  slope_sums <- matrix(0, length(draw_ids), length(groups), dimnames = list(NULL, groups))
  mean_response_sums <- slope_sums
  # Every fitted row is used. The small step approximates a derivative,
  # not a finite distance contrast; halving it was checked for numerical stability.
  batches <- split(seq_len(nrow(model_data)), ceiling(seq_len(nrow(model_data)) / 1000))
  for (ix in batches) {
    prediction_data <- model_data[ix, , drop = FALSE]
    expected_draws <- expected_response(fit, model_spec, prediction_data, re_formula = re_formula, draw_ids = draw_ids)
    lower_predictions <- expected_response(fit, model_spec,
      shift_genetic_distance(prediction_data, -derivative_step), re_formula = re_formula, draw_ids = draw_ids)
    upper_predictions <- expected_response(fit, model_spec,
      shift_genetic_distance(prediction_data, derivative_step), re_formula = re_formula, draw_ids = draw_ids)
    local_slopes <- (upper_predictions - lower_predictions) / (2 * derivative_step)
    for (g in groups) {
      selected <- group[ix] == g
      if (any(selected)) {
        slope_sums[, g] <- slope_sums[, g] + rowSums(local_slopes[, selected, drop = FALSE])
        mean_response_sums[, g] <- mean_response_sums[, g] + rowSums(expected_draws[, selected, drop = FALSE])
      }
    }
  }
  list(level_sums = mean_response_sums, slope_sums = slope_sums)
}

# Models 2 and 3: shared trail calculations -----------------------------------

# Repeat all pairs at each grid value, keeping their other predictors unchanged.
trail_prediction_data <- function(model_data, trail_transformations, variable, values) {
  predictor_scaling <- trail_transformations[match(variable, trail_transformations$variable), ]
  standardized_column <- c(distance_m = "distance_z", FST_Latter = "fst_z", food_difference = "food_z")[[variable]]
  prediction_data <- model_data[rep(seq_len(nrow(model_data)), times = length(values)), , drop = FALSE]
  prediction_data[[variable]] <- rep(values, each = nrow(model_data))
  prediction_data[[standardized_column]] <- (prediction_data[[variable]] - predictor_scaling$center) / predictor_scaling$scale
  prediction_data
}

# Local slopes at every observed row, retaining the other predictors and fitted groups.
trail_local_slopes <- function(fit, model_spec, trail_transformations, variable, derivative_step, draw_ids = seq_len(ndraws(fit))) {
  stopifnot(length(derivative_step) == 1L, is.finite(derivative_step), derivative_step > 0)
  model_data <- model_spec$data
  predictor_scaling <- trail_transformations[match(variable, trail_transformations$variable), ]
  standardized_column <- c(distance_m = "distance_z", FST_Latter = "fst_z", food_difference = "food_z")[[variable]]
  raw_column <- if (variable == "FST_Latter" && endsWith(model_spec$name, "_WC84")) "FST_WC84" else variable
  stopifnot(all(is.finite(model_data[[raw_column]])), is.finite(predictor_scaling$scale), predictor_scaling$scale > 0)
  predict_at <- function(rows, shift) {
    prediction_data <- model_data[rows, , drop = FALSE]
    prediction_data[[raw_column]] <- prediction_data[[raw_column]] + shift
    prediction_data[[standardized_column]] <- (prediction_data[[raw_column]] - predictor_scaling$center) / predictor_scaling$scale
    expected_response(fit, model_spec, prediction_data, re_formula = NULL, draw_ids = draw_ids)
  }
  boundary <- variable %in% c("distance_m", "food_difference") & model_data[[raw_column]] < derivative_step
  local_slopes <- matrix(NA_real_, length(draw_ids), nrow(model_data))
  rows <- which(!boundary)
  if (length(rows)) local_slopes[, rows] <- (predict_at(rows, derivative_step) - predict_at(rows, -derivative_step)) / (2 * derivative_step)
  rows <- which(boundary)
  # A second-order forward derivative keeps food and distance nonnegative at zero.
  if (length(rows)) local_slopes[, rows] <-
    (-3 * predict_at(rows, 0) + 4 * predict_at(rows, derivative_step) - predict_at(rows, 2 * derivative_step)) / (2 * derivative_step)
  local_slopes
}

# Genetic-effect summaries reused by the trail robustness checks.
genetic_robustness_effects <- function(model_fits, model_specs, trail_transformations) {
  imap_dfr(model_fits, function(fit, model_name) {
    local_slopes <- trail_local_slopes(fit, model_specs[[model_name]], trail_transformations, "FST_Latter", derivative_step = 1e-5)
    multiplier <- if (model_specs[[model_name]]$response_scale == "bernoulli") 100 else 1
    draws <- multiplier * .01 * rowMeans(local_slopes)
    label <- if (endsWith(model_name, "_WC84")) "WC84" else if (endsWith(model_name, "_tighter_genetic_prior"))
      "Tighter genetic prior" else "Main (Latter)"
    bind_cols(tibble(Model = label), summarize_columns(matrix(draws, ncol = 1)))
  })
}
