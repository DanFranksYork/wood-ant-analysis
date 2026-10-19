# Author: Daniel Franks

# Priors come from each model specification in Wood_ant_analysis.Rmd.
# Seeds, expected-response curves, main plots and summaries are in the notebook.
# Full simulations here retain shared group effects and observation noise.
# Preserve the draw order so fixed seeds reproduce the displayed prior predictions.

# Model 1: Genetic differentiation simulations ----------------------------------
# Simulate FST with all nest, colony, colony-pair and residual variation.

simulate_genetic_prior <- function(s, draws = 400L) {
  p <- s$prior_parameters
  ii <- unique(round(seq(1, nrow(s$data), length.out = min(800, nrow(s$data)))))
  d <- s$data[ii, , drop = FALSE]
  yrep <- matrix(NA_real_, draws, nrow(d))
  for (i in seq_len(draws)) {
    # Population relationship, shared endpoints, colony combinations, then residuals.
    eta <- rnorm(1, p$intercept[1], p$intercept[2]) + rnorm(1, 0, p$same_colony) * d$same_colony +
      rnorm(1, 0, p$between_slope) * d$distance_between + rnorm(1, 0, p$within_slope) * d$distance_within
    u <- rnorm(nlevels(d$from_label), 0, abs(rnorm(1, 0, p$endpoint_sd)))
    eta <- eta + u[d$from_label] + u[d$to_label]
    v <- rnorm(nlevels(d$from_colony), 0, abs(rnorm(1, 0, p$endpoint_sd)))
    eta <- eta + v[d$from_colony] + v[d$to_colony]
    q <- rnorm(nlevels(d$colony_pair), 0, abs(rnorm(1, 0, p$pair_intercept_sd)))
    eta <- eta + q[d$colony_pair]
    slope <- rnorm(nlevels(d$colony_pair), 0, abs(rnorm(1, 0, p$pair_slope_sd)))
    eta <- eta + slope[d$colony_pair] * d$pair_distance_between
    sigma <- exp(rnorm(1, p$log_sigma[1], p$log_sigma[2]) + rnorm(1, 0, p$sigma_contrast) * d$same_colony)
    yrep[i, ] <- plogis(rnorm(nrow(d), eta, sigma))
  }
  list(data = d, yrep = yrep)
}

# Models 2 and 3: Trail simulations ---------------------------------------------
# Shared predictor design and correlated colony effects; model-specific steps are below.

simulate_prior <- function(s, draws = 400) {
  d <- s$data
  # A fixed systematic subset keeps the prior check inexpensive.
  ii <- unique(round(seq(1, nrow(d), length.out = min(nrow(d), 800))))
  d <- d[ii, , drop = FALSE]
  d$Intercept <- 1
  X <- model.matrix(~ 0 + Intercept + distance_z + fst_z + food_z, d)
  yrep <- matrix(NA_real_, draws, nrow(d))
  investment <- s$response_scale == "lognormal"
  p <- s$prior_parameters
  slope_sd <- p$slope_sd
  for (r in seq_len(draws)) {
    # Draw population coefficients, then shared nest and colony effects.
    b <- rnorm(ncol(X), 0, slope_sd)
    if (investment) b[colnames(X) == "food_z"] <- rnorm(1, 0, p$food_sd)
    b[colnames(X) == "fst_z"] <- b[colnames(X) == "fst_z"] * p$genetic_sd / slope_sd
    if (!investment) b[colnames(X) == "distance_z"] <-
      -abs(b[colnames(X) == "distance_z"]) * p$distance_sd / slope_sd
    b[colnames(X) == "Intercept"] <- rnorm(1, p$intercept[1], p$intercept[2])
    mu <- drop(X %*% b)
    nest_sd <- if (investment) abs(rnorm(1, 0, p$group_sd)) else rexp(1, p$group_rate)
    nest <- rnorm(nlevels(d$from_label), 0, nest_sd)
    mu <- mu + nest[d$from_label] + nest[d$to_label]
    tau <- if (investment) abs(rnorm(2, 0, c(p$group_sd, p$colony_slope_sd))) else rexp(2, p$group_rate)
    rho <- 2 * rbeta(1, p$correlation, p$correlation) - 1 # two-dimensional LKJ
    z <- matrix(rnorm(2 * nlevels(d$colony)), ncol = 2)
    colony_intercept <- tau[1] * z[, 1]
    colony_distance <- tau[2] * (rho * z[, 1] + sqrt(1 - rho^2) * z[, 2])
    mu <- mu + colony_intercept[d$colony]
    if (investment) {
      # Model 3: additive colony distance slopes on log investment.
      mu <- mu + colony_distance[d$colony] * d$distance_z
    } else {
      # Model 2: replace the common slope by -magnitude * exp(colony variation).
      mu <- mu + b[colnames(X) == "distance_z"] *
        (exp(colony_distance[d$colony]) - 1) * d$distance_z
    }
    # Add response variation: binary trails (Model 2) or positive investment (Model 3).
    yrep[r, ] <- switch(s$response_scale,
      bernoulli = rbinom(nrow(d), 1, plogis(mu)),
      lognormal = rlnorm(nrow(d), mu, abs(rnorm(1, 0, p$residual_sd)))
    )
  }
  stopifnot(all(is.finite(yrep)))

  list(data = d, yrep = yrep)
}

# Reusable expected-response calculations --------------------------------------

long_prior_draws <- function(x, y) {
  stopifnot(all(is.finite(y)))
  tibble(draw = rep(seq_len(nrow(y)), times = length(x)),
    x = rep(x, each = nrow(y)), y = as.vector(y))
}

# Integrate FST over residual variation, including unusually broad prior draws.
prior_mean_fst <- function(mu, sigma) {
  quad <- normal_quadrature(100L)
  y <- matrix(0, nrow(mu), ncol(mu))
  for (k in seq_along(quad$nodes)) {
    y <- y + quad$weights[k] * plogis(mu + sigma * quad$nodes[k])
  }
  for (i in which(sigma > 2)) for (j in seq_len(ncol(mu))) {
    y[i, j] <- integrate(function(z) plogis(mu[i, j] + sigma[i] * z) * dnorm(z),
      -Inf, Inf, rel.tol = 1e-8)$value
  }
  y
}

# Model 3: tail summaries for the genetic-prior robustness check -----------------

genetic_prior_tail_summary <- function(simulations) {
  imap_dfr(simulations, function(pr, prior) {
    designs <- list(`All trail designs` = seq_len(ncol(pr$yrep)),
      `Smallest FST` = which.min(pr$data$fst_z), `Largest FST` = which.max(pr$data$fst_z),
      `Largest food predictor` = which.max(pr$data$food_z))
    imap_dfr(designs, function(ii, design) tibble(
      prior = prior, design = design,
      `Investment > 1 (%)` = 100 * mean(pr$yrep[, ii] > 1),
      `Investment > 10 (%)` = 100 * mean(pr$yrep[, ii] > 10)
    ))
  })
}
