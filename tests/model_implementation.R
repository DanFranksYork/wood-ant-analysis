# Author: Daniel Franks

# Run these implementation tests after changing model formulas or endpoint indexing.
# After defining model_specs in the notebook, run:
# source("tests/model_implementation.R")
# check_model_implementation(model_specs)

# Shared model-implementation checks --------------------------------------------
# Models 1 and 3 have linear predictors; Model 2 uses the nonlinear check below.

check_model_implementation <- function(model_specs) {
  endpoint_checks <- imap_dfr(model_specs, function(model_spec, model_name) mutate(verify_mm(model_spec), model = model_name, .before = 1))

  if ("m1" %in% names(model_specs)) {
    stopifnot(!grepl("cholesky_factor_corr|lkj_corr", make_stancode(model_specs$m1$formula,
      data = model_specs$m1$data, family = model_specs$m1$family, prior = model_specs$m1$prior
    )))
  }
  endpoint_checks
}

verify_mm <- function(model_spec, tolerance = 1e-12) {
  if (model_spec$response_scale == "bernoulli") return(verify_presence_mm(model_spec, tolerance))
  model_data <- as.data.frame(model_spec$data)
  stopifnot(nrow(model_data) > 0L, all(model_data$w1 == 1), all(model_data$w2 == 1))
  parsed <- brms::brmsterms(model_spec$formula, family = model_spec$family)$dpars$mu
  groups <- parsed$re
  stopifnot(nrow(groups) > 0L, any(groups$gtype == "mm"))
  standata <- function(data) {
    brms::make_standata(
      model_spec$formula,
      data = data, family = model_spec$family, prior = model_spec$prior
    )
  }
  original_standata <- standata(model_data)
  stopifnot(original_standata$N == nrow(model_data))

  # Give every endpoint a distinct, fixed effect; compare two independent sums.
  effects <- lapply(seq_len(nrow(groups)), function(g) {
    cols <- groups$gcall[[g]]$groups
    labels <- sort(unique(as.character(unlist(model_data[cols], use.names = FALSE))))
    z <- stats::model.matrix(groups$form[[g]], model_data)
    ans <- outer(
      seq_along(labels), seq_len(ncol(z)),
      function(i, j) 0.3 * sin(0.7 * i + j)
    )
    dimnames(ans) <- list(labels, colnames(z))
    ans
  })
  fixed_names <- colnames(original_standata$X)
  fixed_coefficients <- setNames(seq_along(fixed_names) / 7 - 0.3, fixed_names)

  evaluate <- function(data, stan_data) {
    fixed_data <- data
    fixed_data$Intercept <- 1
    fixed_x <- stats::model.matrix(parsed$fe, fixed_data)
    stopifnot(setequal(colnames(fixed_x), fixed_names))
    fixed_x <- fixed_x[, fixed_names, drop = FALSE]
    stopifnot(max(abs(fixed_x - stan_data$X)) < tolerance)
    manual <- drop(fixed_x %*% fixed_coefficients)
    generated <- drop(stan_data$X %*% fixed_coefficients)
    doubled_colony_rows <- 0L
    for (g in seq_len(nrow(groups))) {
      cols <- groups$gcall[[g]]$groups
      z <- stats::model.matrix(groups$form[[g]], data)
      effect <- effects[[g]]
      is_mm <- groups$gtype[g] == "mm"
      if (is_mm) {
        stopifnot(
          length(cols) == 2L, ncol(z) == 1L,
          identical(all.vars(groups$gcall[[g]]$weights), c("w1", "w2")),
          identical(attr(groups$gcall[[g]]$weights, "scale"), FALSE)
        )
      }
      index_names <- if (is_mm) paste0("J_", g, "_", seq_along(cols)) else paste0("J_", g)
      stopifnot(all(index_names %in% names(stan_data)))
      indices <- unlist(stan_data[index_names], use.names = FALSE)
      labels <- as.character(unlist(data[cols], use.names = FALSE))
      mapping <- unique(data.frame(index = indices, label = labels))
      stopifnot(
        !anyDuplicated(mapping$index), !anyDuplicated(mapping$label),
        nrow(mapping) == stan_data[[paste0("N_", g)]]
      )
      mapping <- mapping[order(mapping$index), ]
      stopifnot(identical(mapping$index, seq_len(nrow(mapping))))
      indexed_effect <- effect[mapping$label, , drop = FALSE]
      generated_group <- manual_group <- numeric(nrow(data))
      for (member in seq_along(cols)) {
        j <- stan_data[[index_names[member]]]
        w <- if (is_mm) stan_data[[paste0("W_", g, "_", member)]] else rep(1, nrow(data))
        stopifnot(length(w) == nrow(data), all(w == 1))
        for (k in seq_len(ncol(z))) {
          zname <- if (is_mm) paste0("Z_", g, "_", k, "_", member) else paste0("Z_", g, "_", k)
          stopifnot(max(abs(stan_data[[zname]] - z[, k])) < tolerance)
          generated_group <- generated_group + w * stan_data[[zname]] * indexed_effect[j, k]
          manual_group <- manual_group + z[, k] * effect[as.character(data[[cols[member]]]), k]
        }
      }
      if (is_mm && setequal(cols, c("from_colony", "to_colony"))) {
        same <- as.character(data$from_colony) == as.character(data$to_colony)
        doubled_colony_rows <- sum(same)
        stopifnot(
          doubled_colony_rows > 0L,
          max(abs(generated_group[same] - 2 * effect[as.character(data$from_colony[same]), 1])) < tolerance
        )
      }
      generated <- generated + generated_group
      manual <- manual + manual_group
    }
    stopifnot(max(abs(generated - manual)) < tolerance)
    list(eta = generated, error = max(abs(generated - manual)), doubled = doubled_colony_rows)
  }
  a <- evaluate(model_data, original_standata)
  swapped <- model_data
  for (pair in list(c("from_label", "to_label"), c("from_colony", "to_colony"))) {
    if (all(pair %in% names(model_data))) swapped[pair] <- model_data[rev(pair)]
  }
  swapped[c("w1", "w2")] <- model_data[c("w2", "w1")]
  if ("colony_pair" %in% names(model_data)) {
    first_colony <- as.character(swapped$from_colony)
    second_colony <- as.character(swapped$to_colony)
    swapped$colony_pair <- factor(
      paste(pmin(first_colony, second_colony),
        pmax(first_colony, second_colony),
        sep = " / "
      ),
      levels = levels(model_data$colony_pair)
    )
    stopifnot(identical(swapped$colony_pair, model_data$colony_pair))
  }
  has_colony_term <- any(vapply(groups$gcall, function(x) identical(x$groups, "colony"), logical(1)))
  if (has_colony_term) {
    stopifnot(
      all(as.character(model_data$colony) == as.character(model_data$from_colony)),
      all(as.character(model_data$colony) == as.character(model_data$to_colony))
    )
    swapped$colony <- swapped$from_colony
  }
  swapped_standata <- standata(swapped)
  b <- evaluate(swapped, swapped_standata)
  stopifnot(identical(original_standata$Y, swapped_standata$Y), max(abs(a$eta - b$eta)) < tolerance)
  if (!is.null(original_standata$X_sigma)) {
    stopifnot(identical(original_standata$X_sigma, swapped_standata$X_sigma))
  }
  data.frame(
    manual_predictor_error = a$error,
    exchanged_predictor_error = max(abs(a$eta - b$eta)),
    same_colony_double_contribution_rows = a$doubled,
    unit_endpoint_weights = TRUE
  )
}

# Model 2: Nonlinear presence predictor and endpoint invariance -----------------
# Model 3's linear predictor is covered by verify_mm above.

verify_presence_mm <- function(model_spec, tolerance = 1e-12) {
  model_data <- as.data.frame(model_spec$data)
  model_data$Intercept <- 1
  standata <- function(data) make_standata(model_spec$formula, data = data,
    family = model_spec$family, prior = model_spec$prior)
  original_standata <- standata(model_data)
  fixed <- model.matrix(~ 0 + Intercept + fst_z + food_z, model_data)
  stopifnot(max(abs(original_standata$X_baselineLogOdds - fixed)) < tolerance,
    all(original_standata$X_distanceDecline == 1), all(original_standata$C_1 == model_data$distance_z))
  fixed_coefficients <- c(-2, .4, -.3)
  distanceDecline <- .7
  colonies <- sort(unique(as.character(model_data$colony)))
  nests <- sort(unique(c(as.character(model_data$from_label), as.character(model_data$to_label))))
  colony_intercept <- setNames(.3 * sin(seq_along(colonies)), colonies)
  colony_log_multiplier <- setNames(.2 * cos(seq_along(colonies)), colonies)
  nest <- setNames(.2 * sin(seq_along(nests)), nests)
  # brms may reorder levels; validate its index-to-label mapping explicitly.
  indexed <- function(effect, indices, labels) {
    mapping <- unique(data.frame(index = indices, label = as.character(labels)))
    stopifnot(!anyDuplicated(mapping$index), !anyDuplicated(mapping$label),
      setequal(mapping$label, names(effect)))
    mapping <- mapping[order(mapping$index), ]
    stopifnot(identical(mapping$index, seq_len(nrow(mapping))))
    unname(effect[mapping$label])
  }
  generated <- function(stan_data, data) {
    stopifnot(all(stan_data$W_2_1 == 1), all(stan_data$W_2_2 == 1),
      all(stan_data$Z_1_baselineLogOdds_1 == 1), all(stan_data$Z_1_logColonyMultiplier_2 == 1),
      all(stan_data$Z_2_baselineLogOdds_1_1 == 1), all(stan_data$Z_2_baselineLogOdds_1_2 == 1))
    a <- indexed(colony_intercept, stan_data$J_1, data$colony)
    v <- indexed(colony_log_multiplier, stan_data$J_1, data$colony)
    u <- indexed(nest, c(stan_data$J_2_1, stan_data$J_2_2),
      c(as.character(data$from_label), as.character(data$to_label)))
    drop(stan_data$X_baselineLogOdds %*% fixed_coefficients) + a[stan_data$J_1] +
      stan_data$W_2_1 * u[stan_data$J_2_1] + stan_data$W_2_2 * u[stan_data$J_2_2] -
      drop(stan_data$X_distanceDecline * distanceDecline) * exp(v[stan_data$J_1]) * stan_data$C_1
  }
  manual <- -2 + .4 * model_data$fst_z - .3 * model_data$food_z +
    colony_intercept[as.character(model_data$colony)] +
    nest[as.character(model_data$from_label)] + nest[as.character(model_data$to_label)] -
    distanceDecline * exp(colony_log_multiplier[as.character(model_data$colony)]) * model_data$distance_z
  a <- generated(original_standata, model_data)
  swapped <- model_data
  swapped[c("from_label", "to_label")] <- model_data[c("to_label", "from_label")]
  swapped[c("from_colony", "to_colony")] <- model_data[c("to_colony", "from_colony")]
  swapped[c("w1", "w2")] <- model_data[c("w2", "w1")]
  swapped$colony <- swapped$from_colony
  stopifnot(all(as.character(model_data$from_colony) == as.character(model_data$to_colony)))
  swapped_standata <- standata(swapped)
  b <- generated(swapped_standata, swapped)
  stopifnot(identical(original_standata$Y, swapped_standata$Y),
    max(abs(a - manual)) < tolerance, max(abs(a - b)) < tolerance,
    all(-distanceDecline * exp(colony_log_multiplier) <= 0))
  data.frame(manual_predictor_error = max(abs(a - manual)),
    exchanged_predictor_error = max(abs(a - b)),
    same_colony_double_contribution_rows = 0L, unit_endpoint_weights = TRUE)
}
