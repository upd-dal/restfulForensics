project_path <- function(...) {
   testthat::test_path("..", "..", ...)
}
source(project_path("R", "helper_calcStatistics.R"), 
       local = FALSE)

test_that("standardize_names normalizes names correctly", {
   input <- c("North-East", "  SOUTH  ", "Pop_01", "Site.Three")
   
   result <- standardize_names(input)
   
   expect_equal(
      result,
      c("north east", "south", "pop 01", "site three")
   )
})

test_that("pop_breakdown counts unique samples per population", {
   input <- data.frame(
      Sample = c("S1", "S2", "S3", "S4", "S4"),
      Population = c("A", "A", "B", "B", "B")
   )
   
   result <- pop_breakdown(input, "Population")
   result <- result[order(result$Population), ]
   rownames(result) <- NULL
   
   expect_equal(
      result,
      data.frame(
         Population = c("A", "B"),
         Total = c(2L, 2L)
      )
   )
})


test_that("pop_breakdown errors for a missing column", {
   input <- data.frame(
      Sample = c("S1", "S2"),
      Population = c("A", "B")
   )
   
   expect_error(
      pop_breakdown(input, "MissingColumn"),
      "Ensure specified column is present"
   )
})

make_test_genind <- function() {
   genotypes <- data.frame(
      Locus1 = c(
         "A/A", "A/G", "G/G", "A/A", "A/G", "G/G",
         "A/A", "A/A", "A/G", "G/G", "A/G", "G/G"
      ),
      Locus2 = c(
         "C/C", "C/T", "T/T", "C/T", "C/C", "T/T",
         "C/C", "C/T", "T/T", "C/C", "C/T", "T/T"
      ),
      Locus3 = c(
         "G/G", "G/T", "T/T", "G/G", "G/T", "T/T",
         "G/T", "G/T", "T/T", "G/G", "G/G", "T/T"
      )
   )
   
   populations <- factor(rep(c("PopA", "PopB"), each = 6))
   
   adegenet::df2genind(
      genotypes,
      ploidy = 2,
      sep = "/",
      pop = populations,
      type = "codom"
   )
}


test_that("compute_pop_stats returns all expected components", {
   gen <- make_test_genind()
   
   result <- compute_pop_stats(gen)
   
   expect_named(
      result,
      c("mar_list", "heterozygosity", "ttest", "inbreeding_coeff")
   )
   
   expect_s3_class(result$mar_list, "data.frame")
   expect_s3_class(result$heterozygosity, "data.frame")
   expect_s3_class(result$ttest, "data.frame")
   expect_s3_class(result$inbreeding_coeff, "data.frame")
   
   expect_true(all(c("Population", "Variable", "Value") %in%
                      names(result$heterozygosity)))
   
   expect_true(all(c("Population", "Fis") %in%
                      names(result$inbreeding_coeff)))
   
   expect_true(all(result$heterozygosity$Variable %in% c("Ho", "He", "Ae")))
})


test_that("compute_af returns marker and allele frequency columns", {
   gen <- make_test_genind()
   
   result <- compute_af(gen)
   
   expect_s3_class(result, "data.frame")
   expect_equal(names(result)[1], "markers")
   expect_gt(nrow(result), 0)
   expect_gt(ncol(result), 1)
   
   expect_true(all(vapply(result[-1], is.numeric, logical(1))))
})


test_that("compute_hwe returns the expected result components", {
   gen <- make_test_genind()
   
   result <- compute_hwe(gen, correction = NULL, chains = 100)
   
   expect_named(
      result,
      c(
         "hw_summary",
         "hw_dataframe",
         "hw_mc",
         "loci_HWE_failure",
         "pops_out_of_HWE"
      )
   )
   
   expect_true(all(c("rsID") %in% names(result$hw_summary)))
   expect_true("rsID" %in% names(result$hw_mc))
   expect_true("Population" %in% names(result$hw_dataframe))
   expect_true("rsID" %in% names(result$loci_HWE_failure))
   expect_true("Populations" %in% names(result$pops_out_of_HWE))
   
   expect_gt(nrow(result$hw_summary), 0)
   expect_gt(nrow(result$hw_dataframe), 0)
})


test_that("compute_hwe adds adjusted summary columns when correction is used", {
   gen <- make_test_genind()
   
   result <- compute_hwe(
      gen,
      correction = "bonferroni",
      alpha = 0.05,
      chains = 100
   )
   
   expect_true("Chisq_bonferroni" %in%
                  names(result$loci_HWE_failure))
   expect_true("MC_bonferroni" %in%
                  names(result$loci_HWE_failure))
   
   expect_true("Chisq_bonferroni" %in%
                  names(result$pops_out_of_HWE))
   expect_true("MC_bonferroni" %in%
                  names(result$pops_out_of_HWE))
})


test_that("compute_hwe rejects invalid chain counts", {
   gen <- make_test_genind()
   
   expect_error(
      compute_hwe(gen, chains = -1),
      "non-negative integer"
   )
   
   expect_error(
      compute_hwe(gen, chains = 1.5),
      "non-negative integer"
   )
})


test_that("compute_hwe accepts supported p-value correction methods", {
   gen <- make_test_genind()
   
   for (method in c("holm", "BH", "bonferroni")) {
      result <- compute_hwe(gen, correction = method, chains = 100)
      
      expect_true(
         paste0("Chisq_", method) %in% names(result$loci_HWE_failure)
      )
   }
})


test_that("compute_fst returns matrix and long-format results", {
   gen <- make_test_genind()
   
   result <- compute_fst(gen)
   
   expect_named(result, c("fst_matrix", "fst_dataframe"))
   expect_s3_class(result$fst_dataframe, "data.frame")
   
   expect_true(
      all(c("Site1", "Site2", "Fst") %in% names(result$fst_dataframe))
   )
   
   expect_gt(nrow(result$fst_dataframe), 0)
   expect_true(all(is.finite(result$fst_dataframe$Fst) |
                      is.na(result$fst_dataframe$Fst)))
})



test_that("calc_observed_genotype_freq counts genotypes correctly", {
   input <- data.frame(
      Sample = paste0("S", 1:6),
      Population = c("PopA", "PopA", "PopA", "PopB", "PopB", "PopB"),
      rs1 = c("A/A", "A/G", "A/A", "G/G", "A/G", "G/G"),
      rs2 = c("C/C", "C/T", "N", "C/T", "T/T", NA_character_),
      stringsAsFactors = FALSE
   )
   
   result <- calc_observed_genotype_freq(input)
   
   expect_named(result, c("PopA", "PopB"))
   expect_s3_class(result$PopA, "data.frame")
   expect_true(
      all(c(
         "population", "marker", "genotype",
         "count", "n_genotyped", "observed_freq"
      ) %in% names(result$PopA))
   )
   
   rs1_aa <- result$PopA %>%
      dplyr::filter(marker == "rs1", genotype == "A/A")
   
   expect_equal(rs1_aa$count, 2L)
   expect_equal(rs1_aa$n_genotyped, 3L)
   expect_equal(rs1_aa$observed_freq, 2 / 3)
   
   rs2_pop_a <- result$PopA %>%
      dplyr::filter(marker == "rs2")
   
   expect_equal(sum(rs2_pop_a$count), 2L)
})


test_that("calc_observed_genotype_freq handles a marker with no genotypes", {
   input <- data.frame(
      Sample = c("S1", "S2"),
      Population = c("PopA", "PopA"),
      rs1 = c("N", NA_character_),
      stringsAsFactors = FALSE
   )
   
   result <- calc_observed_genotype_freq(input)
   
   expect_named(result, "PopA")
   expect_equal(nrow(result$PopA), 0L)
})


test_that("calc_expected_genotype_freq calculates diploid frequencies", {
   input <- data.frame(
      markers = c("rs1.A", "rs1.G"),
      PopA = c(0.7, 0.3),
      PopB = c(0.4, 0.6)
   )
   
   result <- calc_expected_genotype_freq(input)
   
   expect_s3_class(result, "data.frame")
   expect_true(all(c(
      "marker", "population", "p", "q",
      "homozygous1", "heterozygous", "homozygous2"
   ) %in% names(result)))
   
   pop_a <- result %>%
      dplyr::filter(marker == "rs1", population == "PopA")
   
   expect_equal(pop_a$p, 0.7)
   expect_equal(pop_a$q, 0.3)
   expect_equal(pop_a$homozygous1, 0.49)
   expect_equal(pop_a$heterozygous, 0.42)
   expect_equal(pop_a$homozygous2, 0.09)
   
   expect_equal(
      pop_a$homozygous1 + pop_a$heterozygous + pop_a$homozygous2,
      1
   )
})


test_that("calc_expected_genotype_freq handles a single allele", {
   input <- data.frame(
      markers = "rs2.C",
      PopA = 1
   )
   
   result <- calc_expected_genotype_freq(input)
   
   expect_equal(result$n_alleles, 1L)
   expect_equal(result$q, 0)
   expect_equal(result$homozygous1, 1)
   expect_equal(result$heterozygous, 0)
   expect_equal(result$homozygous2, 0)
})


test_that("calc_expected_genotype_freq applies the five-event floor", {
   input <- data.frame(
      markers = c("rs1.A", "rs1.G"),
      PopA = c(0.99, 0.01)
   )
   
   result <- calc_expected_genotype_freq(input, pop = 10)
   
   # floor = 5 / (2 * 10) = 0.25
   expect_equal(result$p, 0.99)
   expect_equal(result$q, 0.25)
   expect_equal(result$homozygous1, 0.99^2)
   expect_equal(result$heterozygous, 2 * 0.99 * 0.25)
})


test_that("af_to_long preserves population-specific allele frequencies", {
   input <- data.frame(
      markers = c("rs1.A", "rs1.G"),
      PopA = c(0.7, 0.3),
      PopB = c(0.4, 0.6)
   )
   
   result <- af_to_long(input)
   
   expect_equal(nrow(result), 4L)
   expect_true(all(c("marker", "allele", "population", "freq") %in%
                      names(result)))
   
   value <- result %>%
      dplyr::filter(
         marker == "rs1",
         allele == "G",
         population == "PopB"
      )
   
   expect_equal(value$freq, 0.6)
})

test_that("plot_heterozygosity creates a PNG file", {
   out_dir <- tempfile("het_test_")
   dir.create(out_dir)
   
   input <- data.frame(
      Population = c("PopA", "PopA", "PopA", "PopB", "PopB", "PopB"),
      Variable = c("Ho", "He", "Ae", "Ho", "He", "Ae"),
      Value = c(0.2, 0.3, 1.43, 0.25, 0.35, 1.54)
   )
   
   result <- plot_heterozygosity(input, out_dir)
   
   expect_true(file.exists(result))
   expect_equal(basename(result), "heterozygosity_plot.png")
   expect_gt(file.info(result)$size, 0)
   
   unlink(out_dir, recursive = TRUE)
})


test_that("plot_fst creates a heatmap PNG file", {
   out_dir <- tempfile("fst_test_")
   dir.create(out_dir)
   
   input <- data.frame(
      Site1 = c("PopA", "PopA", "PopB", "PopB"),
      Site2 = c("PopA", "PopB", "PopA", "PopB"),
      Fst = c(0, 0.12, 0.12, 0)
   )
   
   result <- plot_fst(input, out_dir)
   
   expect_true(file.exists(result))
   expect_equal(basename(result), "fst_heatmap.png")
   expect_gt(file.info(result)$size, 0)
   
   # test optional labels
   result_with_values <- plot_fst(
      input,
      out_dir,
      show_values = TRUE
   )
   
   expect_true(file.exists(result_with_values))
   expect_gt(file.info(result_with_values)$size, 0)
   
   unlink(out_dir, recursive = TRUE)
})


test_that("export_pop_results creates an Excel workbook", {
   skip_if_not_installed("openxlsx")
   
   out_dir <- tempfile("export_test_")
   dir.create(out_dir)
   
   stats_matrix <- list(
      mar_list = data.frame(
         MAR = c(1.5, 1.8),
         row.names = c("Locus1", "Locus2")
      ),
      heterozygosity = data.frame(
         Population = c("PopA", "PopA", "PopB", "PopB"),
         Variable = c("Ho", "He", "Ho", "He"),
         Value = c(0.2, 0.3, 0.25, 0.35)
      ),
      ttest = data.frame(
         Locus = c("Locus1", "Locus2"),
         Ho = c(0.2, 0.25),
         Hs = c(0.3, 0.35)
      ),
      inbreeding_coeff = data.frame(
         Population = c("PopA", "PopB"),
         Fis = c(0.1, -0.05)
      )
   )
   
   hw_matrix <- list(
      hw_summary = data.frame(rsID = "Locus1", P = 0.5),
      hw_dataframe = data.frame(Population = "PopA", Locus1 = 0.5),
      hw_mc = data.frame(rsID = "Locus1", PopA = 0.6),
      loci_HWE_failure = data.frame(
         rsID = "Locus1", Chisq = 0, MC = 0
      ),
      pops_out_of_HWE = data.frame(
         Populations = "PopA", Chisq = 0, MC = 0
      )
   )
   
   fst_matrix <- list(
      fst_matrix = list(
         PopA = c(0, 0.12),
         PopB = c(0.12, 0)
      )
   )
   
   allele_freq <- data.frame(
      markers = c("Locus1.A", "Locus1.G"),
      PopA = c(0.7, 0.3)
   )
   
   priv_alleles <- data.frame(
      Population = c("PopA", "PopB"),
      PrivateAlleles = c(1, 2)
   )
   
   result <- export_pop_results(
      allele_freq = allele_freq,
      priv_alleles = priv_alleles,
      stats_matrix = stats_matrix,
      hw_matrix = hw_matrix,
      fst_matrix = fst_matrix,
      dir = out_dir
   )
   
   expect_true(file.exists(result))
   expect_match(basename(result), "^population-statistics-results_.*\\.xlsx$")
   expect_gt(file.info(result)$size, 0)
   
   sheets <- openxlsx::getSheetNames(result)
   
   expect_true("Private Alleles" %in% sheets)
   expect_true("Mean Allelic Richness" %in% sheets)
   expect_true("Heterozygosities" %in% sheets)
   expect_true("Allele Frequencies" %in% sheets)
   expect_true("Hardy-Weinberg Equilibrium" %in% sheets)
   expect_true("Pairwise Fst Matrix" %in% sheets)
   
   unlink(out_dir, recursive = TRUE)
})

test_that("af_to_long reshapes AF table", {
   af <- data.frame(
      markers = c("rs-1.A", "rs1.G", "rs2.C"),
      Population_A = c(0.6, 0.4, 0.7),
      Population_B = c(0.5, 0.5, 0.8)
   )
   
   result <- af_to_long(af)
   expect_equal(nrow(result), 6)
   expect_named(result, c("allele", "marker", "population", "freq"))
   expect_setequal(unique(result$marker), c("rs1", "rs2"))
   expect_setequal(unique(result$allele), c("A", "G", "C"))
   expect_setequal(unique(result$population), c("Population_A", "Population_B"))
   expect_true(is.numeric(result$freq))
   expect_equal(result$freq[result$marker == "rs1" &
                               result$allele == "A" &
                               result$population == "Population_A"], 0.6)
})

test_that("prep_af_table assigns values as expected", {
   af <- data.frame(
      allele = c("A ", "G"),
      marker = c("rs1", "rs1"),
      population = c("Pop1", "Pop1"),
      freq = c("0.6", "0.4")
   )
   
   result <- prep_af_table(af)
   expect_named(result, c("allele", "marker", "population", "freq", "marker_std"))
   expect_equal(result$allele, c("A", "G"))
   expect_equal(result$marker_std, c("rs1", "rs1"))
   expect_type(result$freq, "double")
   expect_equal(result$freq, c(0.6, 0.4))
})

test_that("calc_rmp calculates gt frequencies and rmp", {
   af <- data.frame(
      allele = c("A ", "G", "C", "T"),
      marker = c("rs1", "rs1", "rs2", "rs2"),
      population = rep("Pop1", 4),
      freq = c(0.6, 0.4, 0.7, 0.3)
   )
   
   profile <- data.frame(
      marker = c("rs1", "rs2"),
      genotype = c("A/G", "C/C")
   )
   
   theta = 0.01
   
   result <- calc_rmp(profile = profile, af_table = af, n_ref = 100, theta = theta)
   locus_results <- result$locus_results_by_pop[["Pop1"]]
   summary <- result$summary
   expected_heterozygous <- 2 * 0.6 * 0.4 * (1 - theta)
   expected_homozygous <- 0.7^2 + theta * 0.7 * (1 - 0.7)
   expected_rmp <- expected_heterozygous * expected_homozygous
   
   expect_equal(locus_results$p, c(0.6, 0.7))
   expect_equal(locus_results$q, c(0.4, 0.7))
   expect_equal(
      locus_results$gen_prob,
      c(expected_heterozygous, expected_homozygous)
   )
   expect_equal(summary$population, "Pop1")
   expect_equal(summary$loci_used, 2)
   expect_equal(summary$RMP, expected_rmp)
   expect_equal(summary$reciprocal_RMP, 1/expected_rmp)
   expect_equal(nrow(result$missing), 0)
})

test_that("calc_rmp excludes incomplete genotypes", {
   af <- data.frame(
      allele = c("A", "G", "C", "T"),
      marker = c("rs1", "rs1", "rs2", "rs2"),
      population = rep("Pop1", 4),
      freq = c(0.6, 0.4, 0.7, 0.3)
   )
   
   profile <- data.frame(
      marker = c("rs1", "rs2", "rs3"),
      genotype = c("A/G", "C/", "A/NA")
   )
   
   result <- calc_rmp(profile, af)
   
   expect_equal(nrow(result$missing), 2)
   expect_equal(result$missing$marker, c("rs2", "rs3"))
   expect_equal(
      result$missing$reason,
      c("Second allele missing", "Second allele missing")
   )
   
   # Only the complete rs1 genotype should remain.
   expect_equal(
      nrow(result$locus_results_by_pop[["Pop1"]]),
      1
   )
   expect_equal(result$summary$loci_used, 1)
})

test_that("calc_rmp returns NA when an af is unavailable", {
   af <- data.frame(
      allele = c("A", "G"),
      marker = c("rs1", "rs1"),
      population = c("Pop1", "Pop1"),
      freq = c(0.6, 0.4)
   )
   
   profile <- data.frame(
      marker = "rs1",
      genotype = "A/T"
   )
   
   result <- calc_rmp(profile, af)
   
   expect_true(is.na(result$locus_results_by_pop[["Pop1"]]$gen_prob))
   expect_true(is.na(result$summary$RMP))
   expect_equal(result$summary$loci_used, 0)
   expect_equal(nrow(result$missing), 0)
})