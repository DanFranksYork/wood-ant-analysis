# Author: Daniel Franks

# Reusable expected-response calculations and marginal effects.
# Prediction grids, main numerical results and plotting calls are in Wood_ant_analysis.Rmd.

# Shared numerical helpers ----------------------------------------------------

# normal_quadrature() is also used by prior_checks.R.
# Tables and curves use posterior means and 95% equal-tailed credible intervals.

summarize_columns <- function(x) {
  tibble(
    estimate = apply(x, 2, mean),
    lower = apply(x, 2, quantile, .025), upper = apply(x, 2, quantile, .975)
  )
}

# Nodes and weights integrate over a standard normal variable.
normal_quadrature <- function(n = 30) {
  j <- matrix(0, n, n)
  j[cbind(1:(n - 1), 2:n)] <- sqrt(1:(n - 1))
  j <- j + t(j)
  e <- eigen(j, symmetric = TRUE)
  list(nodes = e$values, weights = e$vectors[1, ]^2)
}

# Shared expected-response prediction -----------------------------------------

# Model 1 integrates over Gaussian residual variation before averaging FST.
# Models 2 and 3 use brms expected responses, including the lognormal mean correction.

expected_response <- function(fit, s, newdata, ids, re_formula = NULL) {
  if (s$response_scale != "logit") {
    return(posterior_epred(fit, newdata = newdata, re_formula = re_formula, draw_ids = ids))
  }
  mu <- posterior_linpred(fit, newdata = newdata, re_formula = re_formula, draw_ids = ids)
  sigma <- posterior_linpred(fit,
    newdata = newdata, dpar = "sigma", transform = TRUE,
    re_formula = re_formula, draw_ids = ids
  )
  mean_fst(mu, sigma)
}

mean_fst <- function(mu, sigma) {
  quad <- normal_quadrature()
  mean_response <- matrix(0, nrow(mu), ncol(mu))
  for (k in seq_along(quad$nodes)) {
    mean_response <- mean_response + quad$weights[k] * plogis(mu + sigma * quad$nodes[k])
  }
  mean_response
}

# Model 1: genetic differentiation calculations -------------------------------

check_prediction_integration <- function() {
  quad <- normal_quadrature()
  # Check the expected FST against an independent adaptive integral.
  grid <- expand.grid(mu = c(-4, -2, -.5), sigma = c(.1, .3, .6, 1))
  error <- apply(grid, 1, function(g) {
    numerical <- sum(quad$weights * plogis(g[1] + g[2] * quad$nodes))
    reference <- integrate(function(z) plogis(g[1] + g[2] * z) * dnorm(z),
      -Inf, Inf, rel.tol = 1e-10)$value
    abs(numerical - reference)
  })
  stopifnot(max(error) < 1e-7)
  invisible(max(error))
}

# Change distance without re-centring the fitted colony-pair slopes or changing endpoints.
shift_genetic_distance <- function(d, increase) {
  d$distGeo_m <- d$distGeo_m + increase
  d$distance_100m <- (d$distGeo_m - 50) / 100
  d$distance_between <- d$distance_100m * (1 - d$same_colony)
  d$distance_within <- d$distance_100m * d$same_colony
  d$pair_distance_between <- (d$distGeo_m - d$pair_mean_distance_m) / 100 * (1 - d$same_colony)
  d
}

# Sum expected FST and local distance slopes in batches over all observed pairs.
# The notebook divides these draw-wise sums by category counts and forms contrasts.
genetic_row_sums <- function(fit, spec, step = .01, ids = seq_len(ndraws(fit)), re_formula = NULL) {
  d <- spec$data
  stopifnot(length(step) == 1, is.finite(step), step > 0, all(d$distGeo_m > step))
  group <- ifelse(d$same_colony == 1, "Within colony", "Between colonies")
  groups <- c("Within colony", "Between colonies")
  slopes <- matrix(0, length(ids), length(groups), dimnames = list(NULL, groups))
  levels <- slopes
  # Every fitted row is used. The small step approximates a derivative,
  # not a finite distance contrast; halving it was checked for numerical stability.
  batches <- split(seq_len(nrow(d)), ceiling(seq_len(nrow(d)) / 1000))
  for (ix in batches) {
    new <- d[ix, , drop = FALSE]
    ep <- expected_response(fit, spec, new, re_formula = re_formula, ids = ids)
    low <- expected_response(fit, spec, shift_genetic_distance(new, -step), re_formula = re_formula, ids = ids)
    high <- expected_response(fit, spec, shift_genetic_distance(new, step), re_formula = re_formula, ids = ids)
    local <- (high - low) / (2 * step)
    for (g in groups) {
      selected <- group[ix] == g
      if (any(selected)) {
        slopes[, g] <- slopes[, g] + rowSums(local[, selected, drop = FALSE])
        levels[, g] <- levels[, g] + rowSums(ep[, selected, drop = FALSE])
      }
    }
  }
  list(level_sums = levels, slope_sums = slopes)
}

# Models 2 and 3: shared trail calculations -----------------------------------

# Repeat all pairs at each grid value, keeping their other predictors unchanged.
trail_prediction_data <- function(d, transformations, variable, values) {
  ss <- transformations[match(variable, transformations$variable), ]
  target <- c(distGeo_m = "distance_z", FST_Latter = "fst_z", food_difference = "food_z")[[variable]]
  new <- d[rep(seq_len(nrow(d)), times = length(values)), , drop = FALSE]
  new[[variable]] <- rep(values, each = nrow(d))
  new[[target]] <- (new[[variable]] - ss$center) / ss$scale
  new
}

# Local slopes at every observed row, retaining the other predictors and fitted groups.
trail_local_slopes <- function(fit, spec, transformations, variable, step, ids = seq_len(ndraws(fit))) {
  stopifnot(length(step) == 1L, is.finite(step), step > 0)
  d <- spec$data
  ss <- transformations[match(variable, transformations$variable), ]
  target <- c(distGeo_m = "distance_z", FST_Latter = "fst_z", food_difference = "food_z")[[variable]]
  column <- if (variable == "FST_Latter" && endsWith(spec$name, "_WC84")) "FST_WC84" else variable
  stopifnot(all(is.finite(d[[column]])), is.finite(ss$scale), ss$scale > 0)
  predict_at <- function(rows, shift) {
    new <- d[rows, , drop = FALSE]
    new[[column]] <- new[[column]] + shift
    new[[target]] <- (new[[column]] - ss$center) / ss$scale
    expected_response(fit, spec, new, re_formula = NULL, ids = ids)
  }
  boundary <- variable %in% c("distGeo_m", "food_difference") & d[[column]] < step
  local <- matrix(NA_real_, length(ids), nrow(d))
  rows <- which(!boundary)
  if (length(rows)) local[, rows] <- (predict_at(rows, step) - predict_at(rows, -step)) / (2 * step)
  rows <- which(boundary)
  # A second-order forward derivative keeps food and distance nonnegative at zero.
  if (length(rows)) local[, rows] <-
    (-3 * predict_at(rows, 0) + 4 * predict_at(rows, step) - predict_at(rows, 2 * step)) / (2 * step)
  local
}

# Genetic-effect summaries reused by the trail robustness checks.
genetic_robustness_effects <- function(fits, specs, transformations) {
  imap_dfr(fits, function(fit, nm) {
    local <- trail_local_slopes(fit, specs[[nm]], transformations, "FST_Latter", step = 1e-5)
    multiplier <- if (specs[[nm]]$response_scale == "bernoulli") 100 else 1
    draws <- multiplier * .01 * rowMeans(local)
    label <- if (endsWith(nm, "_WC84")) "WC84" else if (endsWith(nm, "_tighter_genetic_prior"))
      "Tighter genetic prior" else "Main (Latter)"
    bind_cols(tibble(Model = label), summarize_columns(matrix(draws, ncol = 1)))
  })
}
