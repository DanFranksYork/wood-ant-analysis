# Author: Daniel Franks

# Shared predictor transformations and conditional sampling refinement.
# Formulas, priors and brm() calls are in Wood_ant_analysis.Rmd.

normal_prior <- function(mean, sd) sprintf("normal(%.16g, %.16g)", mean, sd)

prepare_groups <- function(model_data) {
  nest_levels <- sort(unique(c(as.character(model_data$from_label), as.character(model_data$to_label))))
  colony_levels <- sort(unique(na.omit(c(as.character(model_data$from_colony), as.character(model_data$to_colony)))))
  model_data$from_label <- factor(model_data$from_label, levels = nest_levels)
  model_data$to_label <- factor(model_data$to_label, levels = nest_levels)
  model_data$from_colony <- factor(model_data$from_colony, levels = colony_levels)
  model_data$to_colony <- factor(model_data$to_colony, levels = colony_levels)
  model_data$colony <- model_data$from_colony
  model_data$w1 <- model_data$w2 <- 1
  model_data
}

standardizer <- function(x, variable) {
  predictor_sd <- sd(x)
  if (!is.finite(predictor_sd) || predictor_sd <= 0) stop("No variation in ", variable)
  tibble(variable = variable, center = mean(x), scale = predictor_sd)
}

apply_trail_scales <- function(model_data, scales, fst_column = "FST_Latter") {
  for (j in seq_len(nrow(scales))) {
    source <- if (scales$variable[j] == "FST_Latter") fst_column else scales$variable[j]
    target <- c(
      distance_m = "distance_z", FST_Latter = "fst_z",
      food_difference = "food_z"
    )[[scales$variable[j]]]
    model_data[[target]] <- (model_data[[source]] - scales$center[j]) / scales$scale[j]
  }
  prepare_groups(model_data)
}

# Pair-specific deviations are centred within each unordered colony combination.
centre_colony_pair_distance <- function(model_data, reference = model_data) {
  if (!"colony_pair" %in% names(reference)) {
    reference$colony_pair <- with(reference, factor(paste(
      pmin(as.character(from_colony), as.character(to_colony)),
      pmax(as.character(from_colony), as.character(to_colony)), sep = " / "),
      levels = levels(model_data$colony_pair)))
  }
  centres <- reference |>
    group_by(colony_pair) |>
    summarise(pair_mean_distance_m = mean(distance_m), .groups = "drop")
  model_data <- model_data |>
    select(-any_of(c("pair_mean_distance_m", "pair_distance_between"))) |>
    left_join(centres, by = "colony_pair") |>
    mutate(pair_distance_between = (distance_m - pair_mean_distance_m) / 100 * (1 - same_colony))
  stopifnot(
    all(model_data$pair_distance_between[model_data$same_colony == 1] == 0),
    all(is.finite(model_data$pair_mean_distance_m))
  )
  model_data
}
prepare_genetic <- function(model_data) {
  stopifnot(all(is.finite(model_data$FST_Latter) & model_data$FST_Latter > 0 & model_data$FST_Latter < 1))
  model_data <- prepare_groups(model_data)
  model_data |> mutate(
    distance_100m = (distance_m - 50) / 100,
    same_colony = as.integer(same_colony),
    distance_between = distance_100m * (1 - same_colony),
    distance_within = distance_100m * same_colony, response = qlogis(FST_Latter)
  )
}

make_trail_transformations <- function(analysis_data) {
  bind_rows(
    # Keep the calibrated distance prior on its original metre scale.
    standardizer(analysis_data$m2$gps_distance_m, "distance_m"),
    standardizer(analysis_data$m2$FST_Latter, "FST_Latter"),
    standardizer(analysis_data$m2$food_difference, "food_difference")
  )
}

# Sampling refinement ----------------------------------------------------------

refine_sampling_fit <- function(fit, seed, csv_dir) {
  on.exit(unlink(csv_dir, recursive = TRUE))
  sampler_diagnostics <- nuts_params(fit)
  initial_divergences <- sum(sampler_diagnostics$Parameter == "divergent__" & sampler_diagnostics$Value == 1)
  control <- fit$fit@stan_args[[1]]$control
  initial_depth_hits <- sum(sampler_diagnostics$Parameter == "treedepth__" & sampler_diagnostics$Value >= control$max_treedepth)
  if (initial_divergences == 0L && initial_depth_hits == 0L) return(fit)

  # Retain the geometry learned in this run and readapt the step size.
  # Only the posterior from the final sampling run is reported.
  raw_draws <- cmdstanr::read_cmdstan_csv(
    list.files(csv_dir, pattern = "\\.csv$", full.names = TRUE)
  )$post_warmup_draws
  sampler_metadata <- attr(fit$fit, "metadata")
  stopifnot(!is.null(sampler_metadata$inv_metric))
  draws <- posterior::as_draws_array(raw_draws)
  last_draws <- posterior::subset_draws(draws, iteration = dim(draws)[1])
  set.seed(seed + 1L)
  fit <- update(fit,
    backend = "cmdstanr", cores = 4,
    iter = 500L + posterior::niterations(fit), warmup = 500L, seed = seed + 1L,
    init = last_draws, inv_metric = sampler_metadata$inv_metric,
    step_size = unlist(sampler_metadata$step_size),
    control = list(
      adapt_delta = if (initial_divergences > 0L) .9999 else control$adapt_delta,
      max_treedepth = control$max_treedepth + as.integer(initial_depth_hits > 0L),
      metric = sampler_metadata$metadata$metric,
      init_buffer = 0L, window = 0L, term_buffer = 500L
    )
  )
  attr(fit, "initial_divergences") <- initial_divergences
  attr(fit, "initial_depth_hits") <- initial_depth_hits
  fit
}
