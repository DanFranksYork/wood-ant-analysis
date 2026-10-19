# Authors: Sandor Piross and Daniel Franks

# Process the workbook data and select the samples for each model.
# Corrections by Daniel Franks: pair both allele columns before df2genind();
# use both nest sizes in the investment denominator (from_nest_size was duplicated).
# Separately, normalise investment by each pair's mean nest size, not a column-wide mean.
# FST is the response in Model 1 and a predictor in Models 2 and 3.

process_wood_ant_data <- function(input_file) {
  # 1. Workbook inputs and identifier checks -------------------------------------
  stopifnot(file.exists(input_file))
  microsat_raw <- readxl::read_excel(input_file, sheet = "microsat_data")
  nest_raw <- readxl::read_excel(input_file, sheet = "nest_data")
  trails_raw <- readxl::read_excel(input_file, sheet = "trails")
  assert_unique <- function(x, columns, description) {
    if (anyDuplicated(x[columns])) stop("Duplicate identifiers: ", description)
  }
  pair_key <- function(a, b) {
    paste(pmin(as.integer(a), as.integer(b)),
      pmax(as.integer(a), as.integer(b)),
      sep = "_"
    )
  }
  complete_sum <- function(x) if (any(!is.finite(x))) NA_real_ else sum(x)
  # Calculate only valid points: distGeo rejects out-of-range coordinates.
  geographic_distance <- function(lon1, lat1, lon2, lat2, valid) {
    distance <- rep(NA_real_, length(lon1))
    if (any(valid)) {
      distance[valid] <- geosphere::distGeo(
        cbind(lon1[valid], lat1[valid]), cbind(lon2[valid], lat2[valid])
      )
    }
    distance
  }
  assert_unique(microsat_raw, "microsat_sample_ID", "ant samples")
  assert_unique(microsat_raw, c("label", "ant_ID"), "nest and ant")
  assert_unique(nest_raw, "label", "nest labels")
  stopifnot(
    !anyNA(microsat_raw[c("microsat_sample_ID", "label", "ant_ID")]),
    !anyNA(nest_raw$label),
    all(microsat_raw$microsat_sample_ID ==
      paste(microsat_raw$label, microsat_raw$ant_ID, sep = "-")),
    all(microsat_raw$label %in% nest_raw$label)
  )

  loci <- sub("_1$", "", grep("_1$", names(microsat_raw), value = TRUE))
  allele_columns <- as.vector(rbind(paste0(loci, "_1"), paste0(loci, "_2")))
  stopifnot(
    length(loci) == 11L, length(allele_columns) == 22L,
    all(allele_columns %in% names(microsat_raw)),
    all(vapply(microsat_raw[allele_columns], function(x) {
      all(is.na(x) | grepl("^[0-9]{3}$", x))
    }, logical(1)))
  )
  locus_audit <- dplyr::bind_rows(lapply(loci, function(locus) {
    a <- microsat_raw[[paste0(locus, "_1")]]
    b <- microsat_raw[[paste0(locus, "_2")]]
    tibble::tibble(
      locus = locus, diploid_samples = length(a),
      zero_pairs = sum(a == "000" & b == "000", na.rm = TRUE),
      half_zero_pairs = sum(xor(a == "000", b == "000"), na.rm = TRUE),
      other_missing_pairs = sum(is.na(a) | is.na(b)),
      nonzero_alleles = length(setdiff(unique(c(a, b)), c("000", NA_character_))),
      allele_values = paste(sort(unique(c(a, b))), collapse = ", ")
    )
  }))

  # 2. Nest metadata, trails and foraging income ---------------------------------
  # Trail status and food income are needed for Models 2 and 3.
  nests <- nest_raw |>
    dplyr::mutate(
      label = as.character(label),
      colony = dplyr::na_if(as.character(colony), ""),
      ID = dplyr::na_if(as.character(ID), ""),
      known_colony = !is.na(colony),
      nest_key = dplyr::if_else(known_colony & !is.na(ID), paste(colony, ID, sep = "_"), NA_character_),
      has_genotypes = label %in% as.character(microsat_raw$label),
      valid_size = is.finite(nest_size) & nest_size > 0,
      valid_gps = is.finite(lat) & is.finite(lon) & abs(lat) <= 90 & abs(lon) <= 180,
      coordinate_proxy = label %in% c("129", "190")
    )
  assert_unique(dplyr::filter(nests, !is.na(nest_key)), "nest_key", "colony and nest ID")
  stopifnot(
    all(nests$known_colony == (nests$nest_type == "Colony nest")),
    all(!is.na(nests$ID[nests$known_colony]))
  )

  trails <- trails_raw |>
    dplyr::mutate(
      trail_row = dplyr::row_number(),
      colony = as.character(colony),
      from_key = paste(colony, from, sep = "_"),
      to_key = paste(colony, to, sep = "_"),
      unordered_trail_key = paste(colony, pmin(from, to), pmax(from, to), sep = "_")
    )
  stopifnot(
    all(trails$trail_type %in% c("internest", "foraging")),
    !anyNA(trails[c("from", "to", "colony")]),
    all(trails$from != trails$to),
    all(is.na(trails$strength) | trails$strength >= 0)
  )
  assert_unique(trails, "unordered_trail_key", "unordered trail endpoints")
  endpoint_lookup <- nests |>
    dplyr::filter(!is.na(nest_key)) |>
    dplyr::select(nest_key, label, nest_x, nest_y, coordinate_proxy)
  trails <- trails |>
    dplyr::left_join(
      dplyr::rename(endpoint_lookup,
        from_label = label, from_x = nest_x, from_y = nest_y,
        from_coordinate_proxy = coordinate_proxy
      ),
      by = c("from_key" = "nest_key"), relationship = "many-to-one"
    ) |>
    dplyr::left_join(
      dplyr::rename(endpoint_lookup,
        to_label = label, to_x = nest_x, to_y = nest_y,
        to_coordinate_proxy = coordinate_proxy
      ),
      by = c("to_key" = "nest_key"), relationship = "many-to-one"
    ) |>
    dplyr::mutate(
      trail_length_cm = dplyr::if_else(trail_type == "internest",
        sqrt((from_x - to_x)^2 + (from_y - to_y)^2), NA_real_
      ),
      mapped_endpoints = !is.na(from_label) & !is.na(to_label)
    )
  stopifnot(nrow(trails) == nrow(trails_raw))
  trail_mapping_issues <- trails |>
    dplyr::filter((trail_type == "internest" & !mapped_endpoints) |
      (trail_type == "foraging" & is.na(from_label)))

  foraging <- trails |>
    dplyr::filter(trail_type == "foraging", !is.na(from_label)) |>
    dplyr::group_by(from_label) |>
    dplyr::summarise(
      recorded_foraging_trails = dplyr::n(),
      recorded_foraging_sum = complete_sum(strength), .groups = "drop"
    )
  nests <- nests |>
    dplyr::left_join(foraging, by = c("label" = "from_label"), relationship = "one-to-one") |>
    dplyr::mutate(
      recorded_foraging_trails = dplyr::coalesce(recorded_foraging_trails, 0L),
      foraging_status = dplyr::case_when(
        !known_colony ~ "unobserved",
        recorded_foraging_trails == 0L ~ "confirmed_absent",
        is.na(recorded_foraging_sum) ~ "recorded_strength_missing",
        TRUE ~ "recorded"
      ),
      foraging_strength = dplyr::case_when(
        foraging_status == "confirmed_absent" ~ 0,
        foraging_status == "recorded" ~ recorded_foraging_sum,
        TRUE ~ NA_real_
      ),
      foraging_income = dplyr::if_else(valid_size, foraging_strength / nest_size, NA_real_)
    )

  sample_counts <- tibble::tibble(
    quantity = c(
      "Ant samples", "Genotyped nests", "Nests in workbook", "Known-colony nests",
      "Known colonies", "Ungenotyped nests", "Internest trail records", "Mapped internest trails",
      "Foraging trail records", "Unmapped foraging source nests", "Known nests without recorded foraging",
      "Known nests with missing or nonpositive size"
    ),
    n = c(
      nrow(microsat_raw), dplyr::n_distinct(microsat_raw$label), nrow(nests),
      sum(nests$known_colony), dplyr::n_distinct(nests$colony[nests$known_colony]),
      sum(!nests$has_genotypes), sum(trails$trail_type == "internest"),
      sum(trails$trail_type == "internest" & trails$mapped_endpoints),
      sum(trails$trail_type == "foraging"),
      sum(trails$trail_type == "foraging" & is.na(trails$from_label)),
      sum(nests$known_colony & nests$recorded_foraging_trails == 0L),
      sum(nests$known_colony & !nests$valid_size)
    )
  )

  # 3. Unordered nest pairs and predictors ---------------------------------------
  genotyped_labels <- as.character(sort(unique(microsat_raw$label)))
  pair_index <- t(utils::combn(genotyped_labels, 2L))
  dyads <- tibble::tibble(from_label = pair_index[, 1L], to_label = pair_index[, 2L]) |>
    dplyr::mutate(pair_id = pair_key(from_label, to_label))
  endpoint_variables <- c(
    "colony", "ID", "nest_size", "valid_size", "valid_gps", "lat", "lon",
    "foraging_income", "foraging_status", "coordinate_proxy"
  )
  for (endpoint in c("from", "to")) {
    lookup <- dplyr::select(nests, label, dplyr::all_of(endpoint_variables))
    names(lookup) <- paste0(endpoint, "_", names(lookup))
    dyads <- dplyr::left_join(dyads, lookup,
      by = paste0(endpoint, "_label"),
      relationship = "many-to-one"
    )
  }
  expected_pairs <- choose(length(genotyped_labels), 2L)
  stopifnot(nrow(dyads) == expected_pairs)
  internest <- trails |>
    dplyr::filter(trail_type == "internest", mapped_endpoints) |>
    dplyr::mutate(pair_id = pair_key(from_label, to_label)) |>
    dplyr::select(pair_id, trail_row, strength, trail_length_cm)
  assert_unique(internest, "pair_id", "mapped internest pairs")
  dyads <- dyads |>
    dplyr::left_join(internest, by = "pair_id", relationship = "one-to-one") |>
    dplyr::mutate(
      known_colony_pair = !is.na(from_colony) & !is.na(to_colony),
      same_colony = dplyr::if_else(known_colony_pair, from_colony == to_colony, NA),
      listed_trail = !is.na(trail_row),
      trail_status = dplyr::case_when(
        listed_trail ~ "observed_present",
        !is.na(same_colony) & same_colony ~ "confirmed_absent",
        TRUE ~ "unobserved"
      ),
      connected = dplyr::case_when(
        trail_status == "observed_present" ~ TRUE,
        trail_status == "confirmed_absent" ~ FALSE, TRUE ~ NA
      ),
      distGeo_m = geographic_distance(from_lon, from_lat, to_lon, to_lat,
        from_valid_gps & to_valid_gps
      ),
      mean_endpoint_size = dplyr::if_else(from_valid_size & to_valid_size,
        (from_nest_size + to_nest_size) / 2, NA_real_
      ),
      # Models 2 and 3: absolute difference in income per nest population.
      food_difference = abs(from_foraging_income - to_foraging_income),
      # Model 3: relative workforce allocation along an observed trail.
      trail_investment = dplyr::if_else(listed_trail & is.finite(strength) &
        is.finite(trail_length_cm) & !is.na(mean_endpoint_size),
      strength * trail_length_cm / mean_endpoint_size, NA_real_
      ),
      coordinate_proxy = from_coordinate_proxy | to_coordinate_proxy
    )
  stopifnot(
    nrow(dyads) == expected_pairs,
    !anyDuplicated(dyads$pair_id), all(dyads$from_label != dyads$to_label),
    all(dyads$same_colony[dyads$listed_trail]),
    identical(pair_key(dyads$to_label, dyads$from_label), dyads$pair_id)
  )
  # Every symmetric predictor is unchanged when endpoint order is exchanged.
  stopifnot(
    isTRUE(all.equal(
      dyads$food_difference,
      abs(dyads$to_foraging_income - dyads$from_foraging_income)
    )),
    isTRUE(all.equal(
      dyads$mean_endpoint_size,
      ifelse(dyads$from_valid_size & dyads$to_valid_size,
        (dyads$to_nest_size + dyads$from_nest_size) / 2, NA_real_
      )
    )),
    isTRUE(all.equal(dyads$distGeo_m,
      geographic_distance(
        dyads$to_lon, dyads$to_lat, dyads$from_lon, dyads$from_lat,
        dyads$from_valid_gps & dyads$to_valid_gps
      ),
      tolerance = 1e-8
    ))
  )

  # 4. Diploid genotypes and FST for all three models -----------------------------
  make_diploid <- function(swap = FALSE) {
    result <- lapply(loci, function(locus) {
      a <- microsat_raw[[paste0(locus, if (swap) "_2" else "_1")]]
      b <- microsat_raw[[paste0(locus, if (swap) "_1" else "_2")]]
      # Confirmed for this workbook: 000 denotes a failed allele call.
      a[a == "000"] <- NA_character_
      b[b == "000"] <- NA_character_
      dplyr::if_else(is.na(a) | is.na(b), NA_character_,
        paste0("A", a, "/A", b)
      )
    })
    names(result) <- loci
    as.data.frame(result, stringsAsFactors = FALSE)
  }
  diploid_calls <- make_diploid()
  typed_loci <- rowSums(!is.na(diploid_calls))
  keep_ant <- typed_loci > 0L
  genotype_exclusions <- tibble::tibble(
    microsat_sample_ID = microsat_raw$microsat_sample_ID,
    label = as.character(microsat_raw$label), typed_loci = typed_loci,
    included = keep_ant,
    reason = dplyr::if_else(keep_ant, "included", "no_typed_diploid_locus")
  )

  if (!any(keep_ant)) stop("No typed individuals remain.")
  to_genind <- function(calls, order = which(keep_ant)) {
    adegenet::df2genind(calls[order, , drop = FALSE],
      ploidy = 2, type = "codom", sep = "/",
      ind.names = microsat_raw$microsat_sample_ID[order],
      pop = factor(as.character(microsat_raw$label[order]), levels = genotyped_labels),
      NA.char = NA_character_
    )
  }
  genotype <- to_genind(diploid_calls)
  genotype_swapped <- to_genind(make_diploid(swap = TRUE))
  allele_table <- function(g) adegenet::tab(g)[, sort(colnames(adegenet::tab(g))), drop = FALSE]
  stopifnot(
    adegenet::nLoc(genotype) == 11L,
    identical(adegenet::locNames(genotype), loci),
    isTRUE(all.equal(allele_table(genotype), allele_table(genotype_swapped),
      check.attributes = FALSE
    ))
  )
  expected_allele_counts <- vapply(loci, function(locus) {
    calls <- diploid_calls[[locus]][keep_ant]
    length(unique(unlist(strsplit(calls[!is.na(calls)], "/", fixed = TRUE))))
  }, integer(1))
  stopifnot(identical(as.integer(adegenet::nAll(genotype)), as.integer(expected_allele_counts)))
  locus_audit$analysed_alleles <- as.integer(adegenet::nAll(genotype))
  locus_audit$typed_diploid_calls <- colSums(!is.na(diploid_calls))

  # Alphabetic allele labels force positive integer recoding in hierfstat.
  hf_data <- hierfstat::genind2hierfstat(genotype)
  hf_swap <- hierfstat::genind2hierfstat(genotype_swapped)
  # Integer encodings may reverse the two allele blocks without changing the genotype.
  stopifnot(isTRUE(all.equal(hierfstat::pop.freq(hf_data),
    hierfstat::pop.freq(hf_swap),
    check.attributes = FALSE
  )))
  expected_populations <- sort(unique(as.character(microsat_raw$label[keep_ant])))
  methods <- c(FST_Latter = "Fst", FST_WC84 = "WC84", FST_Nei = "Nei87")
  fst_matrices <- lapply(methods, function(method) {
    m <- as.matrix(hierfstat::genet.dist(hf_data, diploid = TRUE, method = method))
    stopifnot(
      identical(rownames(m), colnames(m)),
      setequal(rownames(m), expected_populations),
      isTRUE(all.equal(m, t(m), tolerance = 1e-12)),
      all(diag(m) == 0)
    )
    m
  })

  # Reversing individual order cannot change the population-labelled Latter matrix.
  hf_reversed <- hf_data[rev(seq_len(nrow(hf_data))), , drop = FALSE]
  reverse_latter <- as.matrix(hierfstat::genet.dist(hf_reversed, method = "Fst"))
  original_latter <- fst_matrices$FST_Latter
  stopifnot(isTRUE(all.equal(original_latter,
    reverse_latter[rownames(original_latter), colnames(original_latter)],
    tolerance = 1e-12
  )))
  for (name in names(fst_matrices)) {
    m <- fst_matrices[[name]]
    row_id <- match(dyads$from_label, rownames(m))
    col_id <- match(dyads$to_label, colnames(m))
    dyads[[name]] <- m[cbind(row_id, col_id)]
    stopifnot(isTRUE(all.equal(dyads[[name]], m[cbind(col_id, row_id)], tolerance = 1e-12)))
  }
  fst_support <- dplyr::bind_rows(lapply(names(methods), function(name) {
    x <- dyads[[name]]
    tibble::tibble(
      estimator = name, pairs = length(x), finite = sum(is.finite(x)),
      missing_or_nonfinite = sum(!is.finite(x)), negative = sum(x < 0, na.rm = TRUE),
      zero = sum(x == 0, na.rm = TRUE), one_or_above = sum(x >= 1, na.rm = TRUE),
      minimum = min(x, na.rm = TRUE), median = median(x, na.rm = TRUE), maximum = max(x, na.rm = TRUE)
    )
  }))
  fst_agreement <- dplyr::bind_rows(lapply(c("FST_WC84", "FST_Nei"), function(other) {
    dplyr::bind_rows(lapply(c("all", "within", "between"), function(group) {
      take <- switch(group,
        all = rep(TRUE, nrow(dyads)),
        within = !is.na(dyads$same_colony) & dyads$same_colony,
        between = !is.na(dyads$same_colony) & !dyads$same_colony
      )
      ok <- take & is.finite(dyads$FST_Latter) & is.finite(dyads[[other]])
      a <- dyads$FST_Latter[ok]
      b <- dyads[[other]][ok]
      tibble::tibble(
        estimator = other, pair_group = group, pairs = sum(ok),
        pearson = stats::cor(a, b), spearman = stats::cor(a, b, method = "spearman"),
        mean_difference = mean(b - a), median_absolute_difference = median(abs(b - a)),
        largest_absolute_difference = max(abs(b - a)),
        opposite_sign = sum(a * b < 0)
      )
    }))
  }))

  # 5. Model-specific samples and exclusions -------------------------------------
  # Reasons can overlap; the per-analysis sample table counts each pair once.
  analysis_masks <- list(
    # Model 1: genetic differentiation among pairs with known colonies.
    m1 = list(
      unknown_colony = !dyads$known_colony_pair,
      missing_genetic_estimate = !is.finite(dyads$FST_Latter),
      invalid_geography = !is.finite(dyads$distGeo_m)
    ),
    # Model 2: trail presence among eligible within-colony pairs.
    m2 = list(
      not_known_within_colony = is.na(dyads$same_colony) | !dyads$same_colony,
      unknown_trail_status = is.na(dyads$connected),
      missing_genetic_estimate = !is.finite(dyads$FST_Latter),
      invalid_geography = !is.finite(dyads$distGeo_m),
      unknown_food_difference = !is.finite(dyads$food_difference)
    ),
    # Model 3: positive investment among observed within-colony trails.
    m2b = list(
      not_known_within_colony = is.na(dyads$same_colony) | !dyads$same_colony,
      no_observed_trail = !dyads$listed_trail,
      missing_genetic_estimate = !is.finite(dyads$FST_Latter),
      invalid_geography = !is.finite(dyads$distGeo_m),
      unknown_food_difference = !is.finite(dyads$food_difference),
      missing_or_nonpositive_endpoint_size = !dyads$from_valid_size | !dyads$to_valid_size,
      invalid_positive_investment = !is.finite(dyads$trail_investment) | dyads$trail_investment <= 0
    )
  )
  analysis_exclusions <- dplyr::bind_rows(lapply(names(analysis_masks), function(analysis) {
    dplyr::bind_rows(lapply(names(analysis_masks[[analysis]]), function(reason) {
      mask <- analysis_masks[[analysis]][[reason]]
      tibble::tibble(analysis = analysis, pair_id = dyads$pair_id[which(mask)], reason = reason)
    }))
  }))
  analysis_data <- lapply(analysis_masks, function(reasons) {
    excluded <- Reduce(`|`, reasons)
    dyads[!excluded, , drop = FALSE]
  })
  sample_sizes <- dplyr::bind_rows(lapply(names(analysis_data), function(analysis) {
    x <- analysis_data[[analysis]]
    design_pairs <- switch(analysis,
      m1 = sum(dyads$known_colony_pair),
      m2 = sum(dyads$same_colony, na.rm = TRUE),
      m2b = sum(dyads$listed_trail & !is.na(dyads$same_colony) & dyads$same_colony)
    )
    tibble::tibble(
      analysis = analysis, candidate_pairs = nrow(dyads),
      design_eligible_pairs = design_pairs,
      included_pairs = nrow(x), excluded_pairs = nrow(dyads) - nrow(x),
      measurement_exclusions = design_pairs - nrow(x),
      nests = dplyr::n_distinct(c(x$from_label, x$to_label)),
      known_colonies = dplyr::n_distinct(c(x$from_colony, x$to_colony), na.rm = TRUE),
      observed_trails = sum(x$listed_trail),
      within_pairs = sum(x$same_colony, na.rm = TRUE),
      between_pairs = sum(!x$same_colony, na.rm = TRUE),
      unknown_colony_pairs = sum(is.na(x$same_colony))
    )
  }))

  coordinate_audit <- tibble::tibble(
    analysis = "m2b",
    retained_trails = nrow(analysis_data$m2b),
    trails_with_proxy_coordinates = sum(analysis_data$m2b$coordinate_proxy)
  )

  # 6. Analysis data and the audit tables used in the notebook --------------------
  list(
    analysis_data = analysis_data, dyads = dyads, nests = nests,
    sample_sizes = sample_sizes, sample_counts = sample_counts,
    exclusions = dplyr::count(analysis_exclusions, analysis, reason, name = "excluded_pairs"),
    exclusion_rows = analysis_exclusions, locus_audit = locus_audit,
    genotype_exclusions = genotype_exclusions, coordinate_audit = coordinate_audit,
    trail_mapping_issues = trail_mapping_issues, fst_support = fst_support,
    fst_agreement = fst_agreement, workbook_md5 = unname(tools::md5sum(input_file))
  )
}
