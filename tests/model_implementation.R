# Author: Daniel Franks

# Run these implementation tests after changing model formulas or endpoint indexing.
# After defining specs in the notebook, run:
# source("tests/model_implementation.R")
# check_model_implementation(specs)

# Shared model-implementation checks --------------------------------------------
# Models 1 and 3 have linear predictors; Model 2 uses the nonlinear check below.

check_model_implementation <- function(specs) {
  endpoint_checks <- imap_dfr(specs, function(s, nm) mutate(verify_mm(s), model = nm, .before = 1))

  if ("m1" %in% names(specs)) {
    stopifnot(!grepl("cholesky_factor_corr|lkj_corr", make_stancode(specs$m1$formula,
      data = specs$m1$data, family = specs$m1$family, prior = specs$m1$prior
    )))
  }
  endpoint_checks
}

verify_mm <- function(spec, tolerance = 1e-12) {
  if (spec$response_scale == "bernoulli") return(verify_presence_mm(spec, tolerance))
  d <- as.data.frame(spec$data)
  stopifnot(nrow(d) > 0L, all(d$w1 == 1), all(d$w2 == 1))
  parsed <- brms::brmsterms(spec$formula, family = spec$family)$dpars$mu
  groups <- parsed$re
  stopifnot(nrow(groups) > 0L, any(groups$gtype == "mm"))
  standata <- function(data) {
    brms::make_standata(
      spec$formula,
      data = data, family = spec$family, prior = spec$prior
    )
  }
  original <- standata(d)
  stopifnot(original$N == nrow(d))

  # Give every endpoint a distinct, fixed effect; compare two independent sums.
  effects <- lapply(seq_len(nrow(groups)), function(g) {
    cols <- groups$gcall[[g]]$groups
    labels <- sort(unique(as.character(unlist(d[cols], use.names = FALSE))))
    z <- stats::model.matrix(groups$form[[g]], d)
    ans <- outer(
      seq_along(labels), seq_len(ncol(z)),
      function(i, j) 0.3 * sin(0.7 * i + j)
    )
    dimnames(ans) <- list(labels, colnames(z))
    ans
  })
  fixed_names <- colnames(original$X)
  beta <- setNames(seq_along(fixed_names) / 7 - 0.3, fixed_names)

  evaluate <- function(data, sd) {
    fixed_data <- data
    fixed_data$Intercept <- 1
    fixed_x <- stats::model.matrix(parsed$fe, fixed_data)
    stopifnot(setequal(colnames(fixed_x), fixed_names))
    fixed_x <- fixed_x[, fixed_names, drop = FALSE]
    stopifnot(max(abs(fixed_x - sd$X)) < tolerance)
    manual <- drop(fixed_x %*% beta)
    generated <- drop(sd$X %*% beta)
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
      stopifnot(all(index_names %in% names(sd)))
      indices <- unlist(sd[index_names], use.names = FALSE)
      labels <- as.character(unlist(data[cols], use.names = FALSE))
      mapping <- unique(data.frame(index = indices, label = labels))
      stopifnot(
        !anyDuplicated(mapping$index), !anyDuplicated(mapping$label),
        nrow(mapping) == sd[[paste0("N_", g)]]
      )
      mapping <- mapping[order(mapping$index), ]
      stopifnot(identical(mapping$index, seq_len(nrow(mapping))))
      indexed_effect <- effect[mapping$label, , drop = FALSE]
      generated_group <- manual_group <- numeric(nrow(data))
      for (member in seq_along(cols)) {
        j <- sd[[index_names[member]]]
        w <- if (is_mm) sd[[paste0("W_", g, "_", member)]] else rep(1, nrow(data))
        stopifnot(length(w) == nrow(data), all(w == 1))
        for (k in seq_len(ncol(z))) {
          zname <- if (is_mm) paste0("Z_", g, "_", k, "_", member) else paste0("Z_", g, "_", k)
          stopifnot(max(abs(sd[[zname]] - z[, k])) < tolerance)
          generated_group <- generated_group + w * sd[[zname]] * indexed_effect[j, k]
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
  a <- evaluate(d, original)
  swapped <- d
  for (pair in list(c("from_label", "to_label"), c("from_colony", "to_colony"))) {
    if (all(pair %in% names(d))) swapped[pair] <- d[rev(pair)]
  }
  swapped[c("w1", "w2")] <- d[c("w2", "w1")]
  if ("colony_pair" %in% names(d)) {
    first_colony <- as.character(swapped$from_colony)
    second_colony <- as.character(swapped$to_colony)
    swapped$colony_pair <- factor(
      paste(pmin(first_colony, second_colony),
        pmax(first_colony, second_colony),
        sep = " / "
      ),
      levels = levels(d$colony_pair)
    )
    stopifnot(identical(swapped$colony_pair, d$colony_pair))
  }
  has_colony_term <- any(vapply(groups$gcall, function(x) identical(x$groups, "colony"), logical(1)))
  if (has_colony_term) {
    stopifnot(
      all(as.character(d$colony) == as.character(d$from_colony)),
      all(as.character(d$colony) == as.character(d$to_colony))
    )
    swapped$colony <- swapped$from_colony
  }
  exchanged <- standata(swapped)
  b <- evaluate(swapped, exchanged)
  stopifnot(identical(original$Y, exchanged$Y), max(abs(a$eta - b$eta)) < tolerance)
  if (!is.null(original$X_sigma)) {
    stopifnot(identical(original$X_sigma, exchanged$X_sigma))
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

verify_presence_mm <- function(spec, tolerance = 1e-12) {
  d <- as.data.frame(spec$data)
  d$Intercept <- 1
  standata <- function(data) make_standata(spec$formula, data = data,
    family = spec$family, prior = spec$prior)
  original <- standata(d)
  fixed <- model.matrix(~ 0 + Intercept + fst_z + food_z, d)
  stopifnot(max(abs(original$X_eta - fixed)) < tolerance,
    all(original$X_magnitude == 1), all(original$C_1 == d$distance_z))
  beta <- c(-2, .4, -.3)
  magnitude <- .7
  colonies <- sort(unique(as.character(d$colony)))
  nests <- sort(unique(c(as.character(d$from_label), as.character(d$to_label))))
  colony_intercept <- setNames(.3 * sin(seq_along(colonies)), colonies)
  colony_variation <- setNames(.2 * cos(seq_along(colonies)), colonies)
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
  generated <- function(sd, data) {
    stopifnot(all(sd$W_2_1 == 1), all(sd$W_2_2 == 1),
      all(sd$Z_1_eta_1 == 1), all(sd$Z_1_variation_2 == 1),
      all(sd$Z_2_eta_1_1 == 1), all(sd$Z_2_eta_1_2 == 1))
    a <- indexed(colony_intercept, sd$J_1, data$colony)
    v <- indexed(colony_variation, sd$J_1, data$colony)
    u <- indexed(nest, c(sd$J_2_1, sd$J_2_2),
      c(as.character(data$from_label), as.character(data$to_label)))
    drop(sd$X_eta %*% beta) + a[sd$J_1] +
      sd$W_2_1 * u[sd$J_2_1] + sd$W_2_2 * u[sd$J_2_2] -
      drop(sd$X_magnitude * magnitude) * exp(v[sd$J_1]) * sd$C_1
  }
  manual <- -2 + .4 * d$fst_z - .3 * d$food_z +
    colony_intercept[as.character(d$colony)] +
    nest[as.character(d$from_label)] + nest[as.character(d$to_label)] -
    magnitude * exp(colony_variation[as.character(d$colony)]) * d$distance_z
  a <- generated(original, d)
  swapped <- d
  swapped[c("from_label", "to_label")] <- d[c("to_label", "from_label")]
  swapped[c("from_colony", "to_colony")] <- d[c("to_colony", "from_colony")]
  swapped[c("w1", "w2")] <- d[c("w2", "w1")]
  swapped$colony <- swapped$from_colony
  stopifnot(all(as.character(d$from_colony) == as.character(d$to_colony)))
  exchanged <- standata(swapped)
  b <- generated(exchanged, swapped)
  stopifnot(identical(original$Y, exchanged$Y),
    max(abs(a - manual)) < tolerance, max(abs(a - b)) < tolerance,
    all(-magnitude * exp(colony_variation) <= 0))
  data.frame(manual_predictor_error = max(abs(a - manual)),
    exchanged_predictor_error = max(abs(a - b)),
    same_colony_double_contribution_rows = 0L, unit_endpoint_weights = TRUE)
}
