# Author: Daniel Franks

# Shared notebook display and figure export for all three models.
# Prediction and main plotting code is in Wood_ant_analysis.Rmd.
# Reusable calculations are in predictions.R, prior_checks.R and diagnostics.R.

# Model labels and headline filenames -------------------------------------------

analysis_names <- c(
  m1 = "Genetic differentiation",
  m2 = "Trail presence", m2b = "Trail investment"
)
model_label <- function(name) if (name %in% names(analysis_names)) unname(analysis_names[name]) else name
headline_files <- c(
  m1_distance_comparison = "genetic_differentiation_by_distance",
  m2_effects = "trail_presence_effects", m2b_effects = "trail_investment_effects"
)

# Shared figure style and export ------------------------------------------------
# Within/between colours are used for Model 1.
figure_colours <- c("Within colony" = "#D55E00", "Between colonies" = "#0072B2")
publication_theme <- function() {
  theme_classic(base_size = 8.5, base_family = "Arial") +
    theme(
      text = element_text(colour = "#222222"),
      axis.text = element_text(size = 8, colour = "#333333"),
      axis.title = element_text(size = 8.5),
      axis.title.x = element_text(margin = margin(t = 5)),
      axis.title.y = element_text(margin = margin(r = 5)),
      axis.line = element_line(linewidth = .3, colour = "#555555"),
      axis.ticks = element_line(linewidth = .3, colour = "#555555"),
      axis.ticks.length = grid::unit(1.1, "mm"),
      panel.grid.major.y = element_line(linewidth = .2, colour = "#EAEAEA"),
      plot.title = element_text(size = 9, face = "bold", margin = margin(b = 5)),
      plot.subtitle = element_text(size = 8, margin = margin(b = 6)),
      strip.background = element_blank(),
      strip.text = element_text(size = 8.5, face = "bold", margin = margin(b = 5, t = 3)),
      legend.position = "bottom", legend.title = element_blank(),
      legend.text = element_text(size = 8),
      legend.key.width = grid::unit(7, "mm"),
      legend.key.height = grid::unit(3, "mm"),
      legend.margin = margin(0, 0, 3, 0),
      panel.spacing = grid::unit(5, "mm"),
      plot.margin = margin(7, 8, 6, 6)
    )
}

# Posterior figure formatting --------------------------------------------------

style_genetic_plot <- function(p, groups) {
  p + scale_colour_manual(values = figure_colours, breaks = groups) +
    scale_fill_manual(values = figure_colours, breaks = groups, guide = "none") +
    scale_linetype_manual(values = setNames(c("solid", "22")[seq_along(groups)], groups), breaks = groups) +
    scale_x_continuous(breaks = scales::breaks_pretty(4)) +
    scale_y_continuous(breaks = scales::breaks_pretty(5)) +
    labs(x = "Geographic distance (m)", y = expression(F[ST]), colour = NULL, linetype = NULL)
}

style_trail_panel <- function(p, variable, j, presence, y_range) {
  titles <- c("A   Distance", "B   Genetic differentiation", "C   Food difference")
  xlabels <- list("Geographic distance (m)", expression(F[ST]),
    "Food-income difference\n(10⁻⁴ cm⁻¹)")
  p <- p + scale_x_continuous(breaks = scales::breaks_pretty(4),
    labels = if (variable == "food_difference")
      function(x) scales::number(x * 1e4, accuracy = .1) else scales::label_number(),
    expand = expansion(mult = c(.025, .04))) +
    labs(x = xlabels[[j]], y = if (j == 1L) {
      if (presence) "Probability of a trail" else "Trail investment (log scale)"
    } else NULL, title = titles[j]) +
    theme(plot.margin = margin(6, 5, 6, 4))
  if (presence) {
    p + scale_y_continuous(breaks = seq(0, 1, .25), limits = c(0, 1),
      expand = expansion(mult = c(.025, .04)))
  } else {
    p + scale_y_log10(limits = y_range, labels = function(x)
      vapply(x, function(y) format(y, scientific = FALSE, trim = TRUE), character(1)),
      expand = expansion(mult = c(.035, .05)))
  }
}

arrange_panels <- function(panels, width, height, legend_plot = NULL, align_widths = FALSE) {
  ragg::agg_capture(width = width, height = height, units = "mm", res = 250)
  on.exit(grDevices::dev.off(), add = TRUE)
  if (!is.null(legend_plot)) {
    legend_plot <- ggplotGrob(legend_plot)
    legend <- legend_plot$grobs[[which(legend_plot$layout$name == "guide-box-bottom")]]
    panels <- lapply(panels, function(p) p + theme(legend.position = "none"))
  }
  grobs <- lapply(panels, ggplotGrob)
  if (align_widths) {
    widths <- do.call(grid::unit.pmax, lapply(grobs, function(g) g$widths))
    grobs <- lapply(grobs, function(g) { g$widths <- widths; g })
  }
  heights <- do.call(grid::unit.pmax, lapply(grobs, function(g) g$heights))
  grobs <- lapply(grobs, function(g) { g$heights <- heights; g })
  p <- gridExtra::arrangeGrob(grobs = grobs, nrow = 1)
  if (!is.null(legend_plot)) {
    p <- gridExtra::arrangeGrob(p, legend, ncol = 1,
      heights = grid::unit.c(grid::unit(1, "null"), grid::grobHeight(legend)))
  }
  p
}

# Prior figure formatting ------------------------------------------------------

style_prior_panel <- function(p, variable, j, genetic, investment, y_range) {
  titles <- if (genetic) c("A   Within colony", "B   Between colonies") else
    c("A   Distance", "B   Genetic differentiation", "C   Food difference")
  xlabels <- if (genetic) rep(list("Geographic distance (m)"), 2) else
    list("Geographic distance (m)", expression(F[ST]), "Food-income difference\n(10⁻⁴ cm⁻¹)")
  p <- p + scale_x_continuous(breaks = scales::breaks_pretty(4),
    labels = if (variable == "food_difference") function(x) scales::number(x * 1e4, accuracy = .1)
      else scales::label_number(), expand = expansion(mult = c(.02, .03))) +
    labs(title = titles[j], x = xlabels[[j]], y = if (j == 1L) {
      if (genetic) expression(F[ST]) else if (investment) "Trail investment (log scale)"
        else "Probability of a trail"
    } else NULL)
  if (investment) p + scale_y_log10(limits = y_range, labels = scales::label_number()) else
    p + scale_y_continuous(limits = y_range, breaks = seq(0, 1, .25))
}

initialize_output <- function() {
  dir.create("figures", showWarnings = FALSE)
  theme_set(publication_theme())
  bayesplot::bayesplot_theme_set(publication_theme())
  bayesplot::color_scheme_set("blue")
  figure_files <- new.env(parent = emptyenv())
  figure_temp <- tempfile("wood-ant-figures-")
  dir.create(figure_temp)
  # Embed screen-sized previews; publication files stay in the figures folder.
  save_figure <- function(key, plot, width = 90, height = 78) {
    headline <- key %in% names(headline_files)
    preview <- file.path(figure_temp, paste0(key, ".png"))
    ggsave(preview, plot,
      device = ragg::agg_png, width = width, height = height,
      units = "mm", dpi = 250, bg = "white"
    )
    paths <- NULL
    if (headline) {
      paths <- paste0(file.path("figures", headline_files[[key]]), c(".svg", ".png"))
      ggsave(paths[1], plot, device = svglite::svglite,
        width = width, height = height, units = "mm", bg = "white"
      )
      ggsave(paths[2], plot, device = ragg::agg_png,
        width = width, height = height, units = "mm", dpi = 600, bg = "white"
      )
    }
    figure_files[[key]] <- list(preview = preview, paths = paths, width = width)
  }
  show_figure <- function(key, caption = "") {
    f <- figure_files[[key]]
    files <- if (length(f$paths)) htmltools::tagList(
      htmltools::tags$br(), "Files: ",
      htmltools::tags$a(href = f$paths[1], download = basename(f$paths[1]), basename(f$paths[1])),
      " · ", htmltools::tags$a(href = f$paths[2], download = basename(f$paths[2]), basename(f$paths[2]))
    )
    # RStudio captures printed HTML objects as inline notebook output.
    print(htmltools::tags$figure(
      style = paste0("max-width:", round(f$width / 180 * 880), "px;margin:24px auto"),
      htmltools::tags$img(
        style = "width:100%", alt = caption,
        src = base64enc::dataURI(file = f$preview, mime = "image/png")
      ),
      htmltools::tags$figcaption(
        caption, files
      )
    ))
  }
  list(save = save_figure, show = show_figure)
}

# Shared table display ---------------------------------------------------------

show_table <- function(x, digits = 3, caption = NULL) {
  print(knitr::kable(x, digits = digits, caption = caption, format = "html", row.names = FALSE))
}

# Model 1: genetic result tables ----------------------------------------------

show_genetic_results <- function(predictions) {
  shared <- sprintf("%.1f–%.1f m", predictions$shared_range[1], predictions$shared_range[2])
  show_table(predictions$level_estimates |> transmute(
    Quantity = recode(Quantity, `Within colony` = "Mean FST: within colony",
      `Between colonies` = "Mean FST: between colonies",
      `Within minus between` = "FST contrast: within minus between"),
    `FST (95% credible interval)` = sprintf("%.4f (%.4f to %.4f)", estimate, lower, upper)
  ), caption = paste0("Distance-standardised mean FST and its contrast (shared range: ", shared, ")"))
  show_table(predictions$estimates |> transmute(
    Quantity = recode(Quantity, `Within colony` = "Mean slope: within colony",
      `Between colonies` = "Mean slope: between colonies",
      `Within minus between` = "Slope contrast: within minus between"),
    `FST per m (95% credible interval)` = sprintf("%.5f (%.5f to %.5f)", estimate, lower, upper)
  ), caption = paste0("Average local distance slopes and their contrast (shared range: ", shared, ")"))
  invisible(predictions$estimates)
}

# Models 2 and 3: trail result and robustness tables ---------------------------

show_trail_results <- function(predictions, presence) {
  column <- paste0(if (presence) "AME in trail probability (percentage points)" else
    "AME in mean investment", " (95% credible interval)")
  show_table(predictions$estimates |> transmute(Predictor, Pairs,
    !!column := sprintf(if (presence) "%.2f (%.2f to %.2f)" else "%.5f (%.5f to %.5f)",
      estimate, lower, upper)), caption = "Average marginal effects at the observed predictor values")
  invisible(NULL)
}

show_presence_distance_slopes <- function(predictions) {
  show_table(predictions$estimates |> transmute(
    `Observed distance`, Pairs,
    `Mean local rate, percentage points per m (95% credible interval)` =
      sprintf("%.2f (%.2f to %.2f)", estimate, lower, upper)
  ), caption = "Local distance slopes averaged over observed pairs")
  invisible(predictions$estimates)
}

show_genetic_robustness <- function(fits, specs, transformations) {
  presence <- specs[[1]]$response_scale == "bernoulli"
  effects <- genetic_robustness_effects(fits, specs, transformations)
  column <- paste0(if (presence) "AME in trail probability (percentage points per 0.01 FST)" else
    "AME in mean investment per 0.01 FST", " (95% credible interval)")
  show_table(effects |> transmute(Model,
    !!column := sprintf(if (presence) "%.2f (%.2f to %.2f)" else "%.5f (%.5f to %.5f)",
      estimate, lower, upper)))
  cat("\n\nLocal FST slopes averaged over every analysed pair at its observed values, retaining other predictors and fitted group effects. Zero denotes no local association.\n\n")
  invisible(effects)
}

# Data exclusions and FST checks (all three models) -------------------------------

exclusion_tables <- function(processed) {
  get_count <- function(label) processed$sample_counts$n[match(label, processed$sample_counts$quantity)]
  genotypes <- processed$genotype_exclusions
  usable_nests <- n_distinct(genotypes$label[genotypes$included])
  input <- tibble(
    Records = c("Diploid ant–locus calls", "Ants", "Nests", "Internest trails", "Foraging trails"),
    Available = c(sum(processed$locus_audit$diploid_samples), nrow(genotypes),
      nrow(processed$nests), get_count("Internest trail records"), get_count("Foraging trail records")),
    Retained = c(sum(processed$locus_audit$typed_diploid_calls), sum(genotypes$included),
      usable_nests, get_count("Mapped internest trails"),
      get_count("Foraging trail records") - get_count("Unmapped foraging source nests")),
    Reason = c("At least one failed or missing allele", "No typed locus",
      "No usable genotype sample", "Nest identifier could not be matched",
      "Source nest identifier could not be matched")
  ) |> mutate(Excluded = Available - Retained) |>
    select(Records, Available, Excluded, Retained, Reason)

  reason_labels <- c(
    missing_genetic_estimate = "Unavailable FST", invalid_geography = "Unavailable geographic distance",
    unknown_trail_status = "Unknown trail status", unknown_food_difference = "Unavailable food difference",
    missing_or_nonpositive_endpoint_size = "Missing or nonpositive nest size",
    invalid_positive_investment = "Missing or nonpositive investment"
  )
  size_pairs <- with(processed$dyads, pair_id[!from_valid_size | !to_valid_size])
  missing_size_pairs <- with(processed$dyads,
    pair_id[!is.finite(from_nest_size) | !is.finite(to_nest_size)])
  reasons <- vapply(processed$sample_sizes$analysis, function(name) {
    rows <- filter(processed$exclusion_rows, analysis == name)
    # Separate the analysis's scope from missing measurements; count each pair once.
    outside_scope <- rows$pair_id[rows$reason %in%
      c("unknown_colony", "not_known_within_colony", "no_observed_trail")]
    rows <- filter(rows, !pair_id %in% outside_scope)
    size_related <- rows$pair_id %in% size_pairs & rows$reason %in%
      c("unknown_food_difference", "missing_or_nonpositive_endpoint_size", "invalid_positive_investment")
    rows$reason <- unname(reason_labels[rows$reason])
    rows$reason[size_related] <- "Nonpositive nest size"
    rows$reason[size_related & rows$pair_id %in% missing_size_pairs] <- "Missing nest size"
    stopifnot(!anyNA(rows$reason))
    counts <- rows |>
      group_by(pair_id) |>
      summarise(reason = paste(sort(unique(reason)), collapse = "; "), .groups = "drop") |>
      count(reason, name = "pairs")
    stopifnot(sum(counts$pairs) == processed$sample_sizes$measurement_exclusions[
      processed$sample_sizes$analysis == name])
    if (!nrow(counts)) "None" else if (nrow(counts) == 1L) counts$reason else
      paste(paste0(counts$pairs, ": ", counts$reason), collapse = "; ")
  }, character(1))
  models <- processed$sample_sizes |> transmute(
    Analysis = vapply(analysis, model_label, character(1)),
    `Eligible pairs` = design_eligible_pairs,
    Excluded = measurement_exclusions, `Analysed pairs` = included_pairs,
    Nests = nests, Reason = reasons
  )
  list(input = input, models = models)
}

show_data_audit <- function(processed) {
  exclusions <- exclusion_tables(processed)
  show_table(exclusions$input, caption = "Pre-processing exclusions, by type of record")
  show_table(exclusions$models, caption = "Model samples after applying each analysis's scope")
}

show_fst_audit <- function(processed) {
  show_table(processed$fst_support |> transmute(
    Estimator = estimator, Minimum = minimum, Maximum = maximum,
    Negative = negative, Zero = zero, `At least one` = one_or_above,
    Missing = missing_or_nonfinite
  ), digits = 4, caption = "FST estimates: support and missing values")
  show_table(processed$fst_agreement |> transmute(
    Estimator = recode(estimator, FST_WC84 = "Weir–Cockerham", FST_Nei = "Nei"),
    Pairs = recode(pair_group, all = "All", within = "Within colony", between = "Between colonies"),
    `Correlation with Latter` = pearson,
    `Mean difference from Latter` = mean_difference
  ), digits = 4, caption = "Agreement of alternative FST estimates with Latter")
}
