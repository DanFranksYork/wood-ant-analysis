# Author: Daniel Franks

# Sampling diagnostics and reusable posterior predictive summaries.
# Model identifiers: m1 = Model 1; m2 = Model 2; m2b = Model 3.

# Shared sampling diagnostics ---------------------------------------------------

diagnose_fit <- function(fit, name) {
  posterior_draws <- as_draws_array(fit)
  draw_summary <- summarise_draws(posterior_draws)
  # Constant parameters have undefined R-hat and are excluded explicitly.
  draw_summary <- draw_summary |> filter(variable != "lp__", is.finite(sd), sd > 0)

  sampler_diagnostics <- nuts_params(fit)
  energy <- sampler_diagnostics |>
    filter(Parameter == "energy__") |>
    group_by(Chain) |>
    summarise(ebfmi = mean(diff(Value)^2) / var(Value), .groups = "drop")

  depth_limit <- fit$fit@stan_args[[1]]$control$max_treedepth
  if (is.null(depth_limit)) depth_limit <- 12
  tibble(
    model = name, max_rhat = max(draw_summary$rhat, na.rm = TRUE),
    min_ess_bulk = min(draw_summary$ess_bulk, na.rm = TRUE),
    min_ess_tail = min(draw_summary$ess_tail, na.rm = TRUE),
    divergences = sum(sampler_diagnostics$Parameter == "divergent__" & sampler_diagnostics$Value == 1),
    initial_divergences = if (is.null(attr(fit, "initial_divergences"))) 0L else attr(fit, "initial_divergences"),
    initial_depth_hits = if (is.null(attr(fit, "initial_depth_hits"))) 0L else attr(fit, "initial_depth_hits"),
    treedepth_hits = sum(sampler_diagnostics$Parameter == "treedepth__" & sampler_diagnostics$Value >= depth_limit),
    max_treedepth = depth_limit,
    min_ebfmi = min(energy$ebfmi),
    pass = all(is.finite(draw_summary$rhat), is.finite(draw_summary$ess_bulk), is.finite(draw_summary$ess_tail)) &
      max_rhat < 1.01 & min_ess_bulk > 400 & min_ess_tail > 400 &
      divergences == 0 & treedepth_hits == 0 & min_ebfmi > .3
  )
}

# Shared posterior predictive summaries ----------------------------------------
# Predictions, response transformations and grouping choices are in the notebook.

posterior_predictive_summary <- function(observed_outcomes, replicated_outcomes, name, statistics) {
  stopifnot(ncol(replicated_outcomes) == length(observed_outcomes), all(is.finite(replicated_outcomes)))
  ppc_summary <- tibble()
  for (statistic_name in names(statistics)) {
    statistic_function <- statistics[[statistic_name]]
    observed <- statistic_function(observed_outcomes)
    simulated <- apply(replicated_outcomes, 1, statistic_function)
    ppc_summary <- bind_rows(ppc_summary, tibble(
      model = name, statistic = statistic_name,
      observed = observed, lower_rep = quantile(simulated, .025),
      median_rep = median(simulated), upper_rep = quantile(simulated, .975),
      probability_rep_above_observed = mean(simulated >= observed)
    ))
  }
  ppc_summary
}

posterior_group_checks <- function(observed_outcomes, replicated_outcomes, strata) {
  imap_dfr(strata, function(stratum, name) {
    groups <- split(seq_along(observed_outcomes), stratum)
    outside <- vapply(groups, function(rows) {
      interval <- quantile(rowMeans(replicated_outcomes[, rows, drop = FALSE]), c(.025, .975))
      mean(observed_outcomes[rows]) < interval[1] || mean(observed_outcomes[rows]) > interval[2]
    }, logical(1))
    tibble(
      Grouping = if (name == "colony") "Colony" else "Distance quartile",
      `Groups checked` = length(groups), `Outside 95% predictive interval` = sum(outside)
    )
  })
}

# Model 1: Genetic distributions and residual skewness --------------------------

genetic_predictive_summary <- function(replicated_logit_fst, logit_mean_draws, logit_sd_draws, model_data, name) {
  observed_outcomes <- plogis(model_data$response)
  replicated_outcomes <- plogis(replicated_logit_fst)
  category <- ifelse(model_data$same_colony == 1, "Within", "Between")
  skewness <- function(x) mean((x - mean(x))^3) / stats::sd(x)^3
  groups <- split(seq_len(nrow(model_data)), category)
  results <- list()
  for (group in names(groups)) {
    group_rows <- groups[[group]]
    statistics <- list(
      minimum = min, median = median, maximum = max, sd = stats::sd,
      skewness = skewness
    )
    group_summary <- imap_dfr(statistics, function(statistic_function, statistic) {
      simulated <- apply(replicated_outcomes[, group_rows, drop = FALSE], 1, statistic_function)
      tibble(
        model = name, group = group, statistic = statistic,
        observed = statistic_function(observed_outcomes[group_rows]), lower = quantile(simulated, .025), estimate = median(simulated),
        upper = quantile(simulated, .975), probability_rep_above = mean(simulated >= statistic_function(observed_outcomes[group_rows]))
      )
    })
    # Compare observed and replicated residuals using the same posterior draw.
    observed_residual <- sweep(-logit_mean_draws[, group_rows, drop = FALSE],
      2, model_data$response[group_rows], "+") / logit_sd_draws[, group_rows, drop = FALSE]
    replicated_residual <- (replicated_logit_fst[, group_rows, drop = FALSE] -
      logit_mean_draws[, group_rows, drop = FALSE]) / logit_sd_draws[, group_rows, drop = FALSE]
    observed_skewness <- apply(observed_residual, 1, skewness)
    simulated <- apply(replicated_residual, 1, skewness)
    group_summary <- bind_rows(group_summary, tibble(
      model = name, group = group, statistic = "residual_skewness",
      observed = median(observed_skewness),
      lower = quantile(simulated, .025), estimate = median(simulated), upper = quantile(simulated, .975),
      probability_rep_above = mean(simulated >= observed_skewness)
    ))
    results[[group]] <- group_summary
  }
  bind_rows(results)
}
