project_path <- function(...) {
   testthat::test_path("..", "..", ...)
}
source(project_path("R", "helper_performPCA.R"), 
       local = FALSE)

make_gen_test <- function() {
   geno <- data.frame(
      SNP1 = c("A/A", "A/A", "A/G", "G/G", "A/G", "G/G"),
      SNP2 = c("C/C", "C/T", "C/C", "T/T", "C/T", "T/T"),
      SNP3 = c("G/G", "G/T", "G/G", "T/T", "G/T", "T/T"),
      SNP4 = c("A/C", "A/A", "A/C", "C/C", "A/A", "C/C"),
      check.names = FALSE
   )
   
   gen <- adegenet::df2genind(
      geno,
      sep = "/",
      ploidy = 2,
      type = "codom"
   )
   
   adegenet::indNames(gen) <- paste0("sample", 1:6)
   adegenet::pop(gen) <- factor(
      c("PopA", "PopA", "PopB", "PopB", "PopC", "PopC"),
      levels = c("PopA", "PopB", "PopC")
   )
   gen
}

test_that("compute_pca to return list of dataframes", {
   gen <- make_gen_test()
   result <- compute_pca(gen)
   expect_type(result, "list")
   expect_true(all(c(
      "pca1", "percent", "ind_coords", "centroid"
   ) %in% names(result)))
   
   expect_s3_class(result$pca1, "pca")
   expect_type(result$percent, "double")
   expect_true(length(result$percent) > 0)
   expect_true(is.data.frame(result$ind_coords))
   expect_true(all(c("PC1", "PC2", "Ind", "Site") %in% names(result$ind_coords)))
   
   expect_equal(nrow(result$ind_coords), 6)
   expect_equal(result$ind_coords$Ind,
                paste0("sample", 1:6))
   expect_true(is.data.frame(result$centroid))
   expect_equal(nrow(result$centroid), 3)
   expect_true("Site" %in% names(result$centroid))
})

test_that("compute_pca works without popinfo", {
   gen <- make_gen_test()
   result <- compute_pca(gen, popinfo = FALSE)
   expect_true(is.data.frame(result$ind_coords))
   expect_true(all(c("PC1", "PC2", "Ind") %in% names(result$ind_coords)))
   expect_equal(nrow(result$ind_coords), 6)
   expect_null(result$centroid)
   expect_false("Site" %in% names(result$ind_coords))
})

test_that("compute_pca fails with insufficient SNPs", {
   gen <- make_gen_test()
   gen_loci <- gen[, 1, drop = FALSE]
   expect_error(compute_pca(gen_loci),
                "argument is of length zero")
})

test_that("get_labels returns default", {
   gen <- make_gen_test()
   result <- get_labels(gen)
   expect_type(result, "list")
   expect_true(all(c("labels", "colors", "shapes") %in% names(result)))
   expect_equal(result$labels, c("PopA", "PopB", "PopC"))
   expect_equal(length(result$colors), 3)
   expect_equal(length(result$shapes), 3)
   expect_true(all(result$shapes %in% 21:25))
   expect_true(all(!is.na(result$colors)))
})

test_that("get_labels for custom labels", {
   gen <- make_gen_test()
   
   expect_error(get_labels(gen, use_default = FALSE),
                "Please provide a label file.")
})

test_that("get_labels rejects file with few columns", {
   gen <- make_gen_test()
   label_file <- tempfile(fileext = ".csv")
   write.csv(
      data.frame(
         Population = c("PopA", "PopB", "PopC"),
         Color = c("red", "blue", "green")
      ),
      label_file,
      row.names = FALSE
   )
   
   expect_error(get_labels(
      gen,
      use_default = FALSE,
      label_file = label_file
   ), "Label file must contain at least three columns.")
   
   unlink(label_file)
})

test_that("get_labels reject duplication population names", {
   gen <- make_gen_test()
   label_file <- tempfile(fileext = ".csv")
   
   write.csv(
      data.frame(
         Population = c("PopA", "PopA", "PopC"),
         Color = c("red", "blue", "green"),
         Shape = c(21, 22, 23)
      ),
      label_file,
      row.names = FALSE
   )
   
   expect_error(get_labels(
      gen,
      use_default = FALSE,
      label_file = label_file
   ), "Duplicate population names found.")
   
   unlink(label_file)
})

test_that("get_labels reject mismatched population names", {
   gen <- make_gen_test()
   label_file <- tempfile(fileext = ".csv")
   
   write.csv(
      data.frame(
         Population = c("PopA", "PopB", "PopD"),
         Color = c("red", "blue", "green"),
         Shape = c(21, 22, 23)
      ),
      label_file,
      row.names = FALSE
   )
   
   expect_error(get_labels(
      gen,
      use_default = FALSE,
      label_file = label_file
   ), "Populations in the dataset and given labels do not match.")
   
   unlink(label_file)
})

test_that("get_labels that accept valid file", {
   gen <- make_gen_test()
   label_file <- tempfile(fileext = ".csv")
   
   write.csv(
      data.frame(
         Population = c("PopA", "PopB", "PopC"),
         Color = c("red", "blue", "green"),
         Shape = c(21, 22, 23)
      ),
      label_file,
      row.names = FALSE
   )
   
   result <- get_labels(gen, use_default = FALSE, label_file = label_file)
   expect_equal(result$labels, c("PopA", "PopB", "PopC"))
   expect_equal(unname(result$colors), c("red", "blue", "green"))
   expect_equal(unname(result$shapes), c(21, 22, 23))
   
   unlink(label_file)
})

test_that("plot_pca returns ggplot with pop data",{
   gen <- make_gen_test()
   pca <- compute_pca(gen)
   labels <- get_labels(gen)
   result <- plot_pca(
      ind_coords = pca$ind_coords,
      centroid = pca$centroid,
      percent = pca$percent,
      labels_colors = labels
   )
   
   expect_s3_class(result, "ggplot")
   expect_true(length(result$layers) >= 3)
})

test_that("plot_pca returns ggplot without pop data",{
   gen <- make_gen_test()
   pca <- compute_pca(gen, popinfo = FALSE)

   result <- plot_pca(
      ind_coords = pca$ind_coords,
      centroid = NULL,
      percent = pca$percent,
      labels_colors = NULL
   )
   
   expect_s3_class(result, "ggplot")
})

test_that("plot_pca returns highlight pops",{
   gen <- make_gen_test()
   pca <- compute_pca(gen)
   labels <- get_labels(gen)
   result <- plot_pca(
      ind_coords = pca$ind_coords,
      centroid = pca$centroid,
      percent = pca$percent,
      highlight_pop = "PopA",
      labels_colors = labels,
      centroid_label = FALSE
   )
   
   expect_s3_class(result, "ggplot")
   expect_equal(nrow(result$data), 6)
   expect_true("highlight" %in% names(result$data))
   expect_equal(result$data$highlight, c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE))
})