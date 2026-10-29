project_path <- function(...) {
   testthat::test_path("..", "..", ...)
}
source(project_path("R", "helper_utils_io.R"), 
       local = FALSE)

test_that("Unpack input files", {
   sample_zipped <- project_path("www", "endemic_birds.zip")
   output_dir <- tempfile("io_test")
   dir.create(output_dir)
   
   on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)
   expect_true(file.exists(sample_zipped))
   
   res <- unpack_input_file(sample_zipped, output.dir = output_dir)
   
   expect_type(res, "list")
   expect_true(all(c("data_path", "data_files") %in% names(res)))
   expect_true(dir.exists(res$data_path))
   expect_gt(length(res$data_files), 0)
   expect_true(all(file.exists(res$data_files)))
})

test_that("Read input files - csv", {
   sample_file <- project_path("www", "sample.csv")
   expect_true(file.exists(sample_file))
   res <- load_csv_xlsx_files(sample_file)
   expect_s3_class(res, "data.frame")
   expect_gt(nrow(res), 0)
   expect_gt(ncol(res), 0)
})

test_that("Read input files - excel", {
   sample_file <- project_path("www", "marker_info2.xlsx")
   expect_true(file.exists(sample_file))
   res <- load_csv_xlsx_files(sample_file)
   expect_s3_class(res, "data.frame")
   expect_gt(nrow(res), 0)
   expect_gt(ncol(res), 0)
})

test_that("Read input files - txt", {
   sample_file <- project_path("www", "markers.txt")
   expect_true(file.exists(sample_file))
   res <- load_csv_xlsx_files(sample_file)
   expect_s3_class(res, "data.frame")
   expect_gt(nrow(res), 0)
   expect_gt(ncol(res), 0)
})

test_that("Read input files - vcf", {
   sample_file <- project_path("www", "sample_hgdp.vcf")
   expect_true(file.exists(sample_file))
   res <- load_vcf_files(sample_file)
   expect_s3_class(res, "data.frame")
   expect_gt(nrow(res), 0)
   expect_gt(ncol(res), 0)
})

test_that("Read input files - fasta", {
   sample_file <- project_path("www", "endemic_birds.zip")
   output_dir <- tempfile("fasta_test")
   dir.create(output_dir)
   
   on.exit(unlink(output_dir, recursive = TRUE),
           add = TRUE)
   res <- read_fasta(sample_file, output_dir)
   expect_s4_class(res, "DNAStringSet")
   expect_true(all(width(res)>0))
})

test_that("Read input files - not zipped", {
   expect_error(unpack_input_file(project_path("www", "custom.css")),
                "Not a zipped file.")
})

test_that("Read input files - not valid", {
   expect_error(unpack_input_file(project_path("www", "nothere.zip")),
                "File does not exist in the working directory")
})

test_that("Read input files - invalid for csv/excel", {
   expect_error(load_csv_xlsx_files(project_path("www", "custom.css")),
                "Input file should be in CSV, text, or Excel format.")
})

test_that("Read input files - invalid for vcf", {
   expect_error(load_vcf_files(project_path("www", "sample.csv")),
                "Unsupported file type for VCF input.")
})

