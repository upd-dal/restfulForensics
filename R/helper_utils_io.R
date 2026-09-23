#' Get PLINK2 executable path
#'
#' @returns The file path for PLINK 2.0.
get_plink2_path <- function() {
  return("./plink/plink2.exe")
}

#' Get PLINK1.9 executable path
#'
#' @returns The file path for PLINK 1.9.
get_plink_path <- function() {
  return("./plink/plink.exe")
}

#' Get Arlecore path
#'
#' @returns The file path for the 64-bit arlecore executable
get_arlecore_path <- function() {
  normalizePath("./arlequin/arlecore64.exe",
    winslash = "\\",
    mustWork = TRUE
  )
}


#' Unpack compressed files
#'
#' @param files The file path of the zipped file.
#' @param output.dir The directory to save the unpacked files. Default is the working directory.
#'
#' @returns File path of the directory where files are unpacked and the list of actual files.
unpack_input_file <- function(files, output.dir = ".") {
  if (!file.exists(files)) {
    stop("File does not exist in the working directory")
  } else {
    archive_ext <- tolower(tools::file_ext(files))
    data_path <- tempfile(pattern = "unpacked_", tmpdir = output.dir)
    dir.create(data_path, recursive = TRUE, showWarnings = FALSE)
    
    if (tools::file_ext(files) == "zip") {
      utils::unzip(files,
        files = NULL,
        list = FALSE,
        overwrite = TRUE,
        exdir = data_path
      )

    } else if (tools::file_ext(files) == "tar") {
      untar(files, 
            files = NULL, 
            list = FALSE, 
            exdir = data_path
            )

    } else {
      stop("Not a zipped file. Accepted are zipped (.zip) and tar (.tar) files")
    }

    data_files <- list.files(
      path = data_path, recursive = TRUE, full.names = TRUE, include.dirs = FALSE
    )
    data_files <- data_files[!file.info(data_files)$isdir]
    return(list(data_path = data_path, data_files = data_files))
  }
}

#' Read CSV/XLSX/TXT files
#'
#' @param input The file path of input.
#'
#' @returns A dataframe containing information from the input file.
load_csv_xlsx_files <- function(input) {
  if (!file.exists(input)) {
    stop("File provided does not exist.")
  }

  if (tools::file_ext(input) == "csv") {
    return(read.csv(input, check.names = FALSE))
  } else if (tools::file_ext(input) %in% c("xlsx", "xlsm", "xlsb", "xls")) {
    return(readxl::read_excel(input))
  } else if (tools::file_ext(input) == "txt") {
    return(read.table(input, quote = "\"", comment.char = ""))
  } else {
    stop("Input file should be in CSV, text, or Excel format.")
  }
}

#' Read VCF file and convert to DF object
#'
#' @param vcf The VCF file path.
#' @param output.dir The directory to save the unpacked and merged files if files are zipped. Default is the working directory.
#'
#' @returns The genotype dataframe with sample ID.
load_vcf_files <- function(vcf, output.dir = ".") {
  vcf_ext <- tools::file_ext(vcf)
  if (vcf_ext == "vcf") {
    vcf_object <- vcfR::read.vcfR(vcf, verbose = FALSE)
    genotypes <- vcfR::extract.gt(vcf_object, return.alleles = TRUE)
    columns <- as.data.frame(vcfR::getFIX(vcf_object))
    ID <- columns$ID
    raw_df <- data.frame(ID, genotypes)

    final_df <- data.frame(t(raw_df)) %>%
      janitor::row_to_names(row_number = 1) %>%
      tibble::rownames_to_column(var = "Sample")
  } else if (vcf_ext %in% c("tar", "zip")) {
    vcf_zipped <- unpack_input_file(vcf, output.dir = output.dir)
    wb <- vcf_zipped$data_files
#    wb <- list.files(path = file.path(output.dir), pattern = ".vcf$", full.names = TRUE)
    dflist <- lapply(wb, function(x) {
      vcf_obj <- vcfR::read.vcfR(x, verbose = FALSE)
      genotypes <- vcfR::extract.gt(vcf_obj, return.alleles = TRUE)
      columns <- as.data.frame(vcfR::getFIX(vcf_obj))
      ID <- columns$ID
      temp_df <- data.frame(ID, genotypes)

      data.frame(t(temp_df)) %>%
        janitor::row_to_names(row_number = 1) %>%
        tibble::rownames_to_column(var = "Sample")
    })

    final_df <- dplyr::bind_rows(dflist, .id = "source", .fill = TRUE)
  } else {
    stop("Unsupported file type for VCF input.")
  }
  return(final_df)
}

#' Read FASTA files
#'
#' @param zipped The file path of the zipped file.
#' @param directory The directory to save the unpacked files. Default is the working directory.
#'
#' @returns A DNA bin.
read_fasta <- function(zipped, directory) {
  utils::unzip(zipped,
    files = NULL,
    list = FALSE,
    overwrite = TRUE,
    exdir = file.path(directory, "fasta_files")
  )

  data_path <- file.path(directory, "fasta_files")
  fasta_patterns <- paste("\\.fasta$", "\\.fa$", "\\.fna$", "\\fas$", sep = "|")
  fasta_files <- list.files(path = data_path, pattern = fasta_patterns, full.names = TRUE)
  dna_sequences <- Biostrings::readDNAStringSet(fasta_files)

  return(dna_sequences)
}
