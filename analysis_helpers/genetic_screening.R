# Author: Daniel Franks

# Elva's screening rules, shared by AMOVA and all three Bayesian models.
# Missing diploid genotypes remain missing; no values are imputed.
screen_genotypes <- function(genotypes) {
  metadata <- adegenet::strata(genotypes)
  stopifnot(
    adegenet::nLoc(genotypes) == 11L,
    identical(as.character(metadata$sample_id), adegenet::indNames(genotypes)),
    identical(as.character(metadata$label), as.character(adegenet::pop(genotypes)))
  )
  count_sample <- function(screened_genotypes, stage) {
    tibble::tibble(stage = stage, ants = adegenet::nInd(screened_genotypes),
      nests = length(unique(adegenet::pop(screened_genotypes))),
      colonies = length(unique(stats::na.omit(adegenet::strata(screened_genotypes)$colony))),
      loci = adegenet::nLoc(screened_genotypes))
  }
  counts <- count_sample(genotypes, "Typed workbook genotypes")

  # 1. Keep ants from nests with known colony identity.
  known_colony <- !is.na(metadata$colony)
  unknown_colony <- metadata[!known_colony, , drop = FALSE]
  screened_genotypes <- genotypes[known_colony, ]
  counts <- dplyr::bind_rows(counts, count_sample(screened_genotypes, "Known colonies"))

  # 2. Remove ants missing more than 30% of the original 11 diploid loci.
  ant_missing_fraction <- 1 - adegenet::propTyped(screened_genotypes, by = "ind")
  high_missing <- adegenet::strata(screened_genotypes)[ant_missing_fraction > .30, , drop = FALSE]
  screened_genotypes <- screened_genotypes[ant_missing_fraction <= .30, ]
  counts <- dplyr::bind_rows(counts, count_sample(screened_genotypes, "Ant missingness at most 30%"))

  # 3. Remove loci exceeding 5% missingness after the ant screen.
  locus_missing_fraction <- 1 - adegenet::propTyped(screened_genotypes, by = "loc")
  locus_missingness <- tibble::tibble(locus = names(locus_missing_fraction),
    missing_fraction = unname(locus_missing_fraction), retained = locus_missing_fraction <= .05)
  screened_genotypes <- screened_genotypes[, loc = names(locus_missing_fraction)[locus_missing_fraction <= .05]]
  counts <- dplyr::bind_rows(counts, count_sample(screened_genotypes, "Locus missingness at most 5%"))

  # 4. Remove nests with fewer than four retained ants. Do not repeat the screens.
  nest_counts <- table(as.character(adegenet::pop(screened_genotypes)))
  removed_nests <- names(nest_counts)[nest_counts < 4L]
  keep_nest <- !as.character(adegenet::pop(screened_genotypes)) %in% removed_nests
  low_sample <- adegenet::strata(screened_genotypes)[!keep_nest, , drop = FALSE]
  screened_genotypes <- screened_genotypes[keep_nest, ]
  adegenet::strata(screened_genotypes) <- droplevels(adegenet::strata(screened_genotypes))
  counts <- dplyr::bind_rows(counts, count_sample(screened_genotypes, "At least four ants per nest"))
  counts <- counts |>
    dplyr::mutate(removed_ants = dplyr::lag(ants, default = dplyr::first(ants)) - ants,
      removed_nests = dplyr::lag(nests, default = dplyr::first(nests)) - nests,
      removed_loci = dplyr::lag(loci, default = dplyr::first(loci)) - loci)
  exclusions <- dplyr::bind_rows(
    dplyr::mutate(unknown_colony, reason = "unknown_colony"),
    dplyr::mutate(high_missing, reason = "more_than_30_percent_missing_loci"),
    dplyr::mutate(low_sample, reason = "fewer_than_four_retained_ants_in_nest")
  )
  stopifnot(identical(as.character(adegenet::strata(screened_genotypes)$sample_id), adegenet::indNames(screened_genotypes)),
    all(table(adegenet::pop(screened_genotypes)) >= 4L))
  list(genotypes = screened_genotypes, counts = counts, exclusions = exclusions,
    locus_missingness = locus_missingness,
    nest_counts = tibble::tibble(label = names(nest_counts), ants = as.integer(nest_counts),
      retained = as.logical(nest_counts >= 4L)))
}

rebuild_screened_data <- function(processed_data, genetic_screening) {
  labels <- as.character(unique(adegenet::pop(genetic_screening$genotypes)))
  # Keep the workbook's geographic, food and trail measurements unchanged.
  dyads <- processed_data$dyads |>
    dplyr::filter(from_label %in% labels, to_label %in% labels)
  stopifnot(nrow(dyads) == choose(length(labels), 2L), !anyDuplicated(dyads$pair_id))
  fst_matrices <- calculate_fst_matrices(genetic_screening$genotypes)
  genetics <- attach_fst_estimates(dyads, fst_matrices)
  samples <- select_analysis_samples(genetics$dyads)
  c(samples, list(dyads = genetics$dyads, fst_matrices = fst_matrices,
    fst_support = genetics$fst_support, fst_agreement = genetics$fst_agreement))
}
