project_path <- function(...) {
   testthat::test_path("..", "..", ...)
}
source(project_path("R", "helper_fileConversions.R"), 
       local = FALSE)

test_that("validate_gt_df accepts biallelic genotypes", {
   valid <- data.frame(
      Sample = c("s1", "s2"),
      Population = c("P1", "P2"),
      SNP1 = c("A/G", "T/C"),
      check.names = FALSE
   )
   
   expect_no_error(validate_gt_df(valid))
})

test_that("validate_gt_df rejects non-biallelic genotypes", {
   invalid <- data.frame(
      Sample = c("s1", "s2"),
      Population = c("P1", "P2"),
      SNP1 = c("AString", "123")
   )
   
   expect_error(validate_gt_df(invalid),
                "Unsupported format: Genotype should be present by the third column. Only biallelic markers are accepted.")
})

test_that("clean_input_data to normalize values", {
   input <- data.frame(
      Sample = c("s1", "s2", "s3", "s4"),
      Population = c("P1", "P2", "P2", "P3"),
      SNP1 = c("AG", "A/G", NA, "N/A"),
      SNP2 = c("C|C", ".", "NA", "TT"),
      check.names = FALSE
   ) 
   
   result <- clean_input_data(input)
   expect_true(is.data.frame(result))
   expect_true(all(vapply(result, is.character, logical(1))))
   expect_equal(result$SNP1[3], "N")
   expect_equal(result$SNP1[4], "N")
   expect_equal(result$SNP2[2:3], c("N", "N"))
   expect_equal(result$SNP2[4], "T/T")
   expect_equal(result$SNP1[1], "A/G")
   expect_equal(result$SNP2[1], "C/C")
})

test_that("convert_to_genind to create a genind object", {
   input <- data.frame(
      Sample = c("s1", "s2", "s3", "s4"),
      Population = c("P1", "P2", "P2", "P3"),
      SNP1 = c("A/G", "A/G", "C/G", "G/G"),
      SNP2 = c("C/C", "T/C", "A/A", "T/T"),
      check.names = FALSE
   ) 
   
   result <- convert_to_genind(input, popinfo = TRUE)
   expect_s4_class(result, "genind")
   expect_equal(as.character(adegenet::pop(result)),
                   c("P1", "P2", "P2", "P3"))
   expect_equal(adegenet::nInd(result), 4)
   expect_equal(adegenet::nLoc(result), 2)
})

test_that("convert_to_genind requiring pop info for .str file", {
   input <- data.frame(
      Sample = c("s1", "s2", "s3", "s4"),
      SNP1 = c("A/G", "A/G", "C/G", "G/G"),
      SNP2 = c("C/C", "T/C", "A/A", "T/T"),
      check.names = FALSE
   ) 
   
   expect_error(convert_to_genind(input, to_str = TRUE, popinfo = FALSE),
                "Population data are required when converting to STRUCTURE format.")
})

test_that("subset_genind_pop to keep selected populations only", {
   input <- data.frame(
      Sample = c("s1", "s2", "s3", "s4"),
      Population = c("P1", "P2", "P2", "P3"),
      SNP1 = c("A/G", "A/G", "C/G", "G/G"),
      SNP2 = c("C/C", "T/C", "A/A", "T/T"),
      check.names = FALSE
   ) 
   
   gen <- convert_to_genind(input, popinfo = TRUE)
   res <- subset_genind_pop(gen, "P1")
   expect_s4_class(res, "genind")
   expect_equal(adegenet::nInd(res), 1)
   expect_true(all(as.character(adegenet::pop(res)) == "P1"))
})

test_that("convert_files to reject missing input", {
   temp_dir <- tempfile("missing_input")
   dir.create(temp_dir)
   on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)
   
   expect_error(
      convert_files(file.path(temp_dir, "nonexistent.vcf"), output.dir = temp_dir),
      "File not found."
   )
})

test_that("vcf_to_csv to reject non-VCF file",{
   expect_error(
      vcf_to_csv("input.csv"),
      "Input is not a VCF file."
   )
})

test_that("add_metadata to support one population",{
   input <- data.frame(
      Sample = c("s1", "s2", "s3", "s4"),
      SNP1 = c("A/G", "A/G", "C/G", "G/G"),
      SNP2 = c("C/C", "T/C", "A/A", "T/T"),
      check.names = FALSE
   ) 
   
   result <- add_metadata(input, "Test")
   expect_type(result, "list")
   expect_named(result, c("with_meta", "missing"))
   expect_null(result$missing)
   expect_equal(result$with_meta$Pop, c("Test", "Test", "Test", "Test"))
})

test_that("add_metadata to support multiple population",{
   input <- data.frame(
      Sample = c("s1", "s 2", "s3", "unmatched"),
      SNP1 = c("A/G", "A/G", "C/G", "G/G"),
      SNP2 = c("C/C", "T/C", "A/A", "T/T"),
      check.names = FALSE
   ) 
   
   metadata <- data.frame(
      ID = c("s-1", "S2", "s3"),
      Population = c("P1", "P2", "P2")
   )
   
   result <- add_metadata(input, metadata)
   expect_equal(nrow(result$with_meta), 3)
   expect_equal(nrow(result$missing), 1)
   expect_equal(result$with_meta$Pop, c("P1", "P2", "P2"))
})

test_that("add_metadata errors no matching IDs", {
   input <- data.frame(
      Sample = c("s1"),
      SNP1 = c("A/G"),
      SNP2 = c("C/C"),
      check.names = FALSE
   )
   
   metadata <- data.frame(
      ID = c("random"),
      Population = c("P1")
   )
   
   expect_error(
      add_metadata(input, metadata),
      "Sample IDS do not match between VCF and metadata."
   )
})

test_that("alignment_to_dnabin to reject unsupported extension",{
   expect_error(
      alignment_to_dnabin("alignment.csv"),
      "Unsupported file format: csv"
   )
})

test_that("build_arp_per_population to write Arlequin file", {
   input <- data.frame(
      Sample = c("s1", "s2"),
      Population = c("P1&1", "P2"),
      SNP1 = c("A/G", NA),
      SNP2 = c("C/C", "?"),
      check.names = FALSE
   )
   
   result <- build_arp_per_population(input, output.prefix = "test_arp", data_type = "STANDARD")
   on.exit(unlink(dirname(result), recursive = TRUE), add = TRUE)
   
   expect_true(file.exists(result))
   lines <- readLines(result)
   expect_true(any(lines == "[Profile]"))
   expect_true(any(lines == "NbSamples=2"))
   expect_true(any(lines == 'SampleName="P1and1"'))
   expect_true(any(lines == 'SampleName="P2"'))
   expect_true(any(grepl('MissingData="\\?"', lines)))
})

test_that("csv_to_gtype_format convert missing gts to N/N", {
   input_file <- tempfile(fileext = ".csv")
   on.exit(unlink(input_file), add = TRUE)
   
   utils::write.csv(data.frame(
      Sample = c("s1", "s2"),
      Population = c("P1", "P2"),
      SNP1 = c("N", "A/G"),
      check.names = FALSE
   ),
   input_file,
   row.names = FALSE,
   na = "")
   
   res <- csv_to_gtype_format(input_file)
   expect_true(is.data.frame(res))
   expect_equal(nrow(res), 2)
   expect_equal(res[[3]][1], "N")
   expect_equal(res[[4]][1], "N")
})

