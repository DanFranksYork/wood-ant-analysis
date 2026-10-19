# Author: Daniel Franks

# Sampling diagnostics and reusable posterior predictive summaries.
# Model identifiers: m1 = Model 1; m2 = Model 2; m2b = Model 3.

# Shared sampling diagnostics ---------------------------------------------------

diagnose_fit <- function(fit, name) {
  dr <- as_draws_array(fit)
  sm <- summarise_draws(dr)
  # Constant parameters have undefined R-hat and are excluded explicitly.
  sm <- sm |> filter(variable != "lp__", is.finite(sd), sd > 0)

  np <- nuts_params(fit)
  energy <- np |>
    filter(Parameter == "energy__") |>
    group_by(Chain) |>
    summarise(ebfmi = mean(diff(Value)^2) / var(Value), .groups = "drop")

  depth_limit <- fit$fit@stan_args[[1]]$control$max_treedepth
  if (is.null(depth_limit)) depth_limit <- 12
  tibble(
    model = name, max_rhat = max(sm$rhat, na.rm = TRUE),
    min_ess_bulk = min(sm$ess_bulk, na.rm = TRUE),
    min_ess_tail = min(sm$ess_tail, na.rm = TRUE),
    divergences = sum(np$Parameter == "divergent__" & np$Value == 1),
    initial_divergences = if (is.null(attr(fit, "initial_divergences"))) 0L else attr(fit, "initial_divergences"),
    initial_depth_hits = if (is.null(attr(fit, "initial_depth_hits"))) 0L else attr(fit, "initial_depth_hits"),
    treedepth_hits = sum(np$Parameter == "treedepth__" & np$Value >= depth_limit),
    max_treedepth = depth_limit,
    min_ebfmi = min(energy$ebfmi),
    pass = all(is.finite(sm$rhat), is.finite(sm$ess_bulk), is.finite(sm$ess_tail)) &
      max_rhat < 1.01 & min_ess_bulk > 400 & min_ess_tail > 400 &
      divergences == 0 & treedepth_hits == 0 & min_ebfmi > .3
  )
}

# Shared posterior predictive summaries ----------------------------------------
# Predictions, response transformations and grouping choices are in the notebook.

posterior_predictive_summary <- function(y, yr, name, statistics) {
  stopifnot(ncol(yr) == length(y), all(is.finite(yr)))
  ppc_summary <- tibble()
  for (stat in names(statistics)) {
    f <- statistics[[stat]]
    observed <- f(y)
    simulated <- apply(yr, 1, f)
    ppc_summary <- bind_rows(ppc_summary, tibble(
      model = name, statistic = stat,
      observed = observed, lower_rep = quantile(simulated, .025),
      median_rep = median(simulated), upper_rep = quantile(simulated, .975),
      probability_rep_above_observed = mean(simulated >= observed)
    ))
  }
  ppc_summary
}

posterior_group_checks <- function(y, yr, strata) {
  imap_dfr(strata, function(stratum, name) {
    groups <- split(seq_along(y), stratum)
    outside <- vapply(groups, function(rows) {
      interval <- quantile(rowMeans(yr[, rows, drop = FALSE]), c(.025, .975))
      mean(y[rows]) < interval[1] || mean(y[rows]) > interval[2]
    }, logical(1))
    tibble(
      Grouping = if (name == "colony") "Colony" else "Distance quartile",
      `Groups checked` = length(groups), `Outside 95% predictive interval` = sum(outside)
    )
  })
}

# Model 1: Genetic distributions and residual skewness --------------------------

genetic_predictive_summary <- function(zr, mu, sigma, d, name) {
  y <- plogis(d$response)
  yr <- plogis(zr)
  category <- ifelse(d$same_colony == 1, "Within", "Between")
  skewness <- function(x) mean((x - mean(x))^3) / stats::sd(x)^3
  groups <- split(seq_len(nrow(d)), category)
  results <- list()
  for (group in names(groups)) {
    ii <- groups[[group]]
    statistics <- list(
      minimum = min, median = median, maximum = max, sd = stats::sd,
      skewness = skewness
    )
    out <- imap_dfr(statistics, function(f, statistic) {
      simulated <- apply(yr[, ii, drop = FALSE], 1, f)
      tibble(
        model = name, group = group, statistic = statistic,
        observed = f(y[ii]), lower = quantile(simulated, .025), estimate = median(simulated),
        upper = quantile(simulated, .975), probability_rep_above = mean(simulated >= f(y[ii]))
      )
    })
    # Compare observed and replicated residuals using the same posterior draw.
    observed_residual <- sweep(-mu[, ii, drop = FALSE], 2, d$response[ii], "+") / sigma[, ii, drop = FALSE]
    replicated_residual <- (zr[, ii, drop = FALSE] - mu[, ii, drop = FALSE]) / sigma[, ii, drop = FALSE]
    observed_skewness <- apply(observed_residual, 1, skewness)
    simulated <- apply(replicated_residual, 1, skewness)
    out <- bind_rows(out, tibble(
      model = name, group = group, statistic = "residual_skewness",
      observed = median(observed_skewness),
      lower = quantile(simulated, .025), estimate = median(simulated), upper = quantile(simulated, .975),
      probability_rep_above = mean(simulated >= observed_skewness)
    ))
    results[[group]] <- out
  }
  bind_rows(results)
}
