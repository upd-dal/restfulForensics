project_path <- function(...) {
   testthat::test_path("..", "..", ...)
}
source(project_path("R", "helper_snpExtraction.R"), 
       local = FALSE)

expect_nonempty <- function(path) {
   testthat::expect_true(file.exists(path))
   testthat::expect_gt(file.info(path)$size, 0)
}

test_that("create_range_file creates a nonempty file", {
   pos_input <- data.frame(
      rsID = c("rs1", "rs2", "rs3"),
      CHR = c(1, 2, 3),
      POS = c(123, 456, 789)
   )
   
   result <- create_range_file(pos_input, addID = FALSE)
   expect_type(result, "list")
   expect_true("range_file" %in% names(result))
   expect_nonempty(result$range_file)
   
   range_data <- read.delim(result$range_file, header = FALSE)
   expect_equal(nrow(range_data), 3)
   expect_equal(ncol(range_data), 3)
   expect_equal(range_data[[1]], c(1, 2, 3))
   expect_equal(range_data[[2]], c(123, 456, 789))
   expect_equal(range_data[[3]], c(123, 456, 789))
   unlink(result$range_file)
})

test_that("create_range_file creates range and ID file", {
   pos_input <- data.frame(
      rsID = c("rs1", "rs2"),
      CHR = c(1, 2),
      POS = c(123, 456)
   )
   
   result <- create_range_file(pos_input, addID = TRUE)
   expect_true(all(c("range_file", "updated_file") %in% names(result)))
   expect_nonempty(result$range_file)
   expect_nonempty(result$updated_file)
   
   updated_data <- read.delim(
      result$updated_file,
      header = FALSE,
      col.names = c("new", "id")
   )
   
   expect_equal(updated_data$new, c("1:123", "2:456"))
   expect_equal(updated_data$id, c("rs1", "rs2"))
   unlink(c(result$range_file, result$updated_file))
})

test_that("plot_concordance returns result and plot", {
   concordance_input <- data.frame(
      ID = c("s1", "s2", "s3"),
      SNP1.x = c("A/A", "C/T", "N"),
      SNP1.y = c("A/A", "C/C", "G/G"),
      SNP2.x = c("G/G", "T/T", "A/G"),
      SNP2.y = c("G/G", "T/T", "A/A")
   )
   
   result <- plot_concordance(concordance_input)
   expect_type(result, "list")
   expect_true(all(c("results", "plot") %in% names(result)))
   expect_true(is.data.frame(result$results))
   expect_gt(nrow(result$results), 0)
   expect_s3_class(result$plot, "ggplot")
})