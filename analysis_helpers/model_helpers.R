# Author: Daniel Franks

# Shared predictor transformations and conditional sampling refinement.
# Formulas, priors and brm() calls are in Wood_ant_analysis.Rmd.

normal_prior <- function(mean, sd) sprintf("normal(%.16g, %.16g)", mean, sd)

prepare_groups <- function(d) {
  nl <- sort(unique(c(as.character(d$from_label), as.character(d$to_label))))
  cl <- sort(unique(na.omit(c(as.character(d$from_colony), as.character(d$to_colony)))))
  d$from_label <- factor(d$from_label, levels = nl)
  d$to_label <- factor(d$to_label, levels = nl)
  d$from_colony <- factor(d$from_colony, levels = cl)
  d$to_colony <- factor(d$to_colony, levels = cl)
  d$colony <- d$from_colony
  d$w1 <- d$w2 <- 1
  d
}

standardizer <- function(x, variable) {
  s <- sd(x)
  if (!is.finite(s) || s <= 0) stop("No variation in ", variable)
  tibble(variable = variable, center = mean(x), scale = s)
}

apply_trail_scales <- function(d, scales, fst_column = "FST_Latter") {
  for (j in seq_len(nrow(scales))) {
    source <- if (scales$variable[j] == "FST_Latter") fst_column else scales$variable[j]
    target <- c(
      distGeo_m = "distance_z", FST_Latter = "fst_z",
      food_difference = "food_z"
    )[[scales$variable[j]]]
    d[[target]] <- (d[[source]] - scales$center[j]) / scales$scale[j]
  }
  prepare_groups(d)
}

# Pair-specific deviations are centred within each unordered colony combination.
centre_colony_pair_distance <- function(d) {
  centres <- d |>
    group_by(colony_pair) |>
    summarise(pair_mean_distance_m = mean(distGeo_m), .groups = "drop")
  d <- d |>
    select(-any_of(c("pair_mean_distance_m", "pair_distance_between"))) |>
    left_join(centres, by = "colony_pair") |>
    mutate(pair_distance_between = (distGeo_m - pair_mean_distance_m) / 100 * (1 - same_colony))
  stopifnot(
    all(d$pair_distance_between[d$same_colony == 1] == 0),
    max(abs(tapply(d$pair_distance_between, d$colony_pair, mean))) < 1e-12
  )
  d
}
prepare_genetic <- function(d) {
  stopifnot(all(is.finite(d$FST_Latter) & d$FST_Latter > 0 & d$FST_Latter < 1))
  d <- prepare_groups(d)
  d |> mutate(
    distance_100m = (distGeo_m - 50) / 100,
    same_colony = as.integer(same_colony),
    distance_between = distance_100m * (1 - same_colony),
    distance_within = distance_100m * same_colony, response = qlogis(FST_Latter)
  )
}

make_trail_transformations <- function(analysis_data) {
  bind_rows(
    standardizer(analysis_data$m2$distGeo_m, "distGeo_m"),
    standardizer(analysis_data$m2$FST_Latter, "FST_Latter"),
    standardizer(analysis_data$m2$food_difference, "food_difference")
  )
}

# Sampling refinement ----------------------------------------------------------

refine_sampling_fit <- function(fit, seed, csv_dir) {
  on.exit(unlink(csv_dir, recursive = TRUE))
  np <- nuts_params(fit)
  initial_divergences <- sum(np$Parameter == "divergent__" & np$Value == 1)
  control <- fit$fit@stan_args[[1]]$control
  initial_depth_hits <- sum(np$Parameter == "treedepth__" & np$Value >= control$max_treedepth)
  if (initial_divergences == 0L && initial_depth_hits == 0L) return(fit)

  # Retain the geometry learned in this run and readapt the step size.
  # Only the posterior from the final sampling run is reported.
  raw_draws <- cmdstanr::read_cmdstan_csv(
    list.files(csv_dir, pattern = "\\.csv$", full.names = TRUE)
  )$post_warmup_draws
  sampler <- attr(fit$fit, "metadata")
  stopifnot(!is.null(sampler$inv_metric))
  draws <- posterior::as_draws_array(raw_draws)
  last_draws <- posterior::subset_draws(draws, iteration = dim(draws)[1])
  set.seed(seed + 1L)
  fit <- update(fit,
    backend = "cmdstanr", cores = 4,
    iter = 500L + posterior::niterations(fit), warmup = 500L, seed = seed + 1L,
    init = last_draws, inv_metric = sampler$inv_metric,
    step_size = unlist(sampler$step_size),
    control = list(
      adapt_delta = if (initial_divergences > 0L) .9999 else control$adapt_delta,
      max_treedepth = control$max_treedepth + as.integer(initial_depth_hits > 0L),
      metric = sampler$metadata$metric,
      init_buffer = 0L, window = 0L, term_buffer = 500L
    )
  )
  attr(fit, "initial_divergences") <- initial_divergences
  attr(fit, "initial_depth_hits") <- initial_depth_hits
  fit
}
