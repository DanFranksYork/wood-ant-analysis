# Authors: Elva (AMOVA analysis)

# Input: the screened genind returned by screen_genotypes() in genetic_screening.R.
# Ant, locus and nest exclusions are applied there, once for AMOVA and FST.

run_wood_ant_amova <- function(genotypes, seed = 20261021, permutations = 999) {
  metadata <- adegenet::strata(genotypes)
  stopifnot(
    inherits(genotypes, "genind"),
    identical(as.character(metadata$sample_id), adegenet::indNames(genotypes)),
    !anyNA(metadata$colony), !anyNA(metadata$label),
    all(adegenet::ploidy(genotypes) == 2)
  )
  adegenet::strata(genotypes) <- droplevels(metadata)

  # Partition differences among diploid ants: colonies, nests and ants within nests.
  # Within-ant heterozygosity is not a separate component; missing calls are not imputed.
  # poppr looks up its correction by name in the calling environment.
  quasieuclid <- ade4::quasieuclid
  fit <- poppr::poppr.amova(
    genotypes, ~colony/label, within = FALSE, freq = TRUE, missing = "ignore",
    correction = "quasieuclid", method = "ade4", quiet = TRUE
  )
  set.seed(seed)
  permutation_test <- ade4::randtest(fit, nrepet = permutations)

  variation_levels <- c("Between colonies", "Between nests within colonies",
              "Between ants within nests", "Total")
  amova_summary <- data.frame(
    Source = variation_levels, df = fit$results$Df,
    sum_squares = fit$results$`Sum Sq`, mean_squares = fit$results$`Mean Sq`
  )
  variance_components <- data.frame(
    Source = variation_levels, variance = fit$componentsofcovariance$Sigma,
    percent = fit$componentsofcovariance$`%`
  )
  phi_statistics <- data.frame(
    Comparison = c("Nests relative to the total", "Nests within colonies",
                   "Colonies relative to the total"),
    Phi = fit$statphi$Phi
  )
  # ade4 orders permutation tests from the lowest to the highest level.
  permutation_results <- data.frame(
    Source = rev(variation_levels[1:3]), alternative = permutation_test$alter,
    p_value = permutation_test$pvalue, permutations = permutations
  )
  list(fit = fit, permutation_test = permutation_test, summary = amova_summary,
       variance = variance_components, phi = phi_statistics, permutation = permutation_results)
}
