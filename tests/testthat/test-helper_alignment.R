project_path <- function(...) {
   testthat::test_path("..", "..", ...)
}
source(project_path("R", "helper_alignment.R"), 
       local = FALSE)

test_that("calc_msa returns aligned DNA sequences", {
   sequences <- Biostrings::DNAStringSet(c(
      Sample1 = "ACGTACGT",
      Sample2 = "ACGTTCGT",
      Sample3 = "ACGTACGA"
   ))
   
   result <- calc_msa(files = sequences, algorithm = "ClustalW")
   expect_named(result, c("alignment", "scores", "adjusted", "staggered"))
   
   expect_s4_class(result$alignment, "MsaDNAMultipleAlignment")
   expect_s4_class(result$adjusted, "DNAStringSet")
   expect_s4_class(result$staggered, "DNAStringSet")
   expect_setequal(names(result$adjusted), names(sequences))
})

test_that("build_nj_tree returns plot", {
   sequences <- Biostrings::DNAStringSet(c(
      Sample1 = "ACGTACGT",
      Sample2 = "ACGTTCGT",
      Sample3 = "ACGTACGA",
      Sample4 = "TCGTACGA"
   ))
   
   expect_warning(
      result <- build_nj_tree(alignment = sequences,
                              outgroup = "Sample1",
                              model = "raw"),
      "Alignment has fewer than 10 sites. Skipping bootstrap."
   )
   
   expect_s3_class(result, "ggplot")
   expect_true("Sample3" %in% result$data$label)
})


test_that("build_upgma_tree returns plot", {
   sequences <- Biostrings::DNAStringSet(c(
      Sample1 = "ACGTACGT",
      Sample2 = "ACGTTCGT",
      Sample3 = "ACGTACGA",
      Sample4 = "TCGTACGA"
   ))
   
   expect_warning(
      result <- build_upgma_tree(alignment = sequences,
                              outgroup = "Sample1",
                              model = "raw"),
      "Alignment has fewer than 10 sites. Skipping bootstrap."
   )
   
   expect_s3_class(result, "ggplot")
   expect_true("Sample3" %in% result$data$label)
})

