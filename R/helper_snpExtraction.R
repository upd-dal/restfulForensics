#' Create POS range file for merging
#' @description
#'
#' Creates a range file necessary when extracting SNPs by base pair position.
#' The rsIDs are assumed to be absent from input file. The function also creates a file for adding rsID if addID set to TRUE.
#' Default is FALSE.
#'
#' @param pos_input The path to file containing marker information. Required columns are 1. rsID, 2. chromosome, 3. position if addID is TRUE. Else, only the last two columns are expected.
#' @param addID Indicate whether rsIDs will be added to the extracted SNPs.
#' @param output_dir The directory to save the range files.
#'
#' @returns A list of file paths of range files.
create_range_file <- function(pos_input, addID = FALSE) {
  if (is.character(pos_input)) {
    pos_df <- load_csv_xlsx_files(pos_input)
  } else {
    pos_df <- pos_input
  }

  pos_df <- as.data.frame(pos_df)

  if (isTRUE(addID)) {
    chr <- as.character(trimws(pos_df[[2]]))
    pos <- as.numeric(pos_df[[3]])
    label <- as.character(trimws(pos_df[[1]]))

    range_df <- data.frame(
      CHR = chr,
      START = pos,
      END = pos,
      LABEL = label,
      stringsAsFactors = FALSE
    )

    range_file <- tempfile(pattern = "range", fileext = ".txt")

    write.table(
      range_df,
      file = range_file,
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE,
      sep = "\t"
    )

    update_name <- data.frame(
      new = paste0(chr, ":", pos, sep = ""),
      id = label,
      stringsAsFactors = FALSE
    )

    updated_file <- tempfile(pattern = "update_name", fileext = ".txt")

    write.table(
      update_name,
      file = updated_file,
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE,
      sep = "\t"
    )

    return(list(range_file = range_file, updated_file = updated_file))
  } else {
    chr <- as.character(pos_df[[2]])
    pos <- as.numeric(pos_df[[3]])

    range_df <- data.frame(
      CHR = chr,
      START = pos,
      END = pos,
      stringsAsFactors = FALSE
    )

    range_file <- tempfile(pattern = "range", fileext = ".txt")

    write.table(
      range_df,
      file = range_file,
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE,
      sep = "\t"
    )
    return(list(range_file = range_file))
  }
}


#' Extract SNPs by rsID
#'
#' @param pgen_prefix The prefix of the PLINK2.0 files.
#' @param snps_list The dataframe of the list of SNPs to be extracted.
#' @param output_dir The directory to save output file. Default is working directory.
#' @param merged_file The name of the output file.
#'
#' @returns The file path of result.
extract_by_ID_pgen <- function(pgen_prefix,
                               snps_list,
                               output_dir = ".",
                               merged_file = "extracted_file") {
  plink_path <- get_plink2_path()
  out_prefix <- file.path(output_dir, merged_file)

  cmd <- paste(
    shQuote(plink_path),
    "--pfile", shQuote(pgen_prefix),
    "--extract", shQuote(snps_list),
    "--export vcf",
    "--out", shQuote(out_prefix)
  )

  # system(cmd)
  status <- system(cmd, intern = FALSE, ignore.stdout = FALSE, ignore.stderr = FALSE)
  vcf_file <- paste0(out_prefix, ".vcf")

  if (status != 0) {
    stop("PLINK failed: ", status)
  }

  return(vcf_file)
}


#' Extract SNPs by GRCh38/37 positions
#'
#' @param pos_list The file path of the list of SNPs to be extracted.
#' @param pgen_prefix The prefix of the PLINK2.0 files.
#' @param output_dir The directory to save output file. Default is working directory.
#' @param merged_file The name of the output file.
#'
#' @returns The file path of extracted file.
extract_by_pos_pgen <- function(pos_list,
                                pgen_prefix,
                                output_dir = ".",
                                merged_file = "extracted_file") {
  plink_path <- get_plink2_path()
  range_file <- create_range_file(pos_list, addID = FALSE)

  out_prefix <- file.path(output_dir, merged_file)

  print(read.delim(range_file$range_file, header = FALSE))

  cmd <- paste(
    shQuote(plink_path),
    "--pfile", shQuote(pgen_prefix),
    "--extract", "range", shQuote(range_file$range_file),
    "--export vcf",
    "--out", shQuote(out_prefix)
  )

  status <- system(cmd, intern = FALSE, ignore.stdout = FALSE, ignore.stderr = FALSE)
  vcf_file <- paste0(out_prefix, ".vcf")

  if (status != 0) {
    stop("PLINK failed: ", status)
  }

  return(vcf_file)
}


#' Extract SNPs by GRCh38/37 positions and add rsID
#'
#' @param pos_list The file path of the list of SNPs to be extracted.
#' @param pgen_prefix The prefix of the PLINK2.0 files.
#' @param output_dir The directory to save output file. Default is working directory.
#'
#' @returns The file path of result.
extract_POStoID_pgen <- function(pos_list,
                                 pgen_prefix,
                                 output_dir = ".") {
  plink_path <- get_plink2_path()
  range_file <- create_range_file(pos_list, addID = TRUE)

  extracted_prefix <- file.path(output_dir, "pos_extract")
  cmd_extract <- paste(
    shQuote(plink_path),
    "--pfile", shQuote(pgen_prefix),
    "--extract", "range", shQuote(range_file$range_file),
    "--make-pgen",
    "--out", shQuote(extracted_prefix)
  )
  system(cmd_extract)

  renamed_prefix <- file.path(output_dir, "updated_file")
  cmd_rename <- paste(
    shQuote(plink_path),
    "--pfile", shQuote(extracted_prefix),
    "--set-all-var-ids @:#",
    "--make-pgen",
    "--out", shQuote(renamed_prefix)
  )
  system(cmd_rename)

  updated_prefix <- file.path(output_dir, "extracted_file")
  cmd_updated <- paste(
    shQuote(plink_path),
    "--pfile", shQuote(renamed_prefix),
    "--update-name", shQuote(range_file$updated_file),
    "--recode vcf",
    "--out", shQuote(updated_prefix)
  )
  system(cmd_updated)

  vcf_file <- paste0(updated_prefix, ".vcf")

  if (!file.exists(vcf_file)) {
    stop("PLINK extraction failed: no VCF generated.")
  }

  return(vcf_file)
}


#' Match samples for concordance analysis
#'
#' @param file1 The file path of sequencing results using technique 1. Expected sample and genotype information.
#' @param file2 The file path of sequencing results using technique 2. Expected sample and genotype information.
#' @param phased Indicates if genotypes are phased. Default is FALSE.
#'
#' @returns The data frame containing the genotype of samples present in both file1 and file2.
calc_concordance <- function(file1, file2, phased = FALSE) {
  file1 <- clean_input_data(file1)
  file2 <- clean_input_data(file1)
  
  file1 <- dplyr::rename(file1, Ind = 1)
  file2 <- dplyr::rename(file2, Ind = 1)
  file_list <- list(file1, file2)
  overlaps <- as.list(intersect(file1$Ind, file2$Ind))

  if (is.null(overlaps)) {
    stop("Samples names from the two files are different.")
  }

  file_list2 <- lapply(
    file_list,
    function(x) {
      x[x$Ind %in% overlaps, ]
    }
  )

  file_list3 <- lapply(
    file_list2,
    function(x) {
      library(dplyr)
      data.frame(t(x)) %>%
        janitor::row_to_names(row_number = 1) %>%
        tibble::rownames_to_column(., var = "markers")
    }
  )

  overlaps <- as.character(overlaps)
  markers1 <- file_list3[[1]]$markers
  markers2 <- file_list3[[2]]$markers

  file_list4 <- lapply(
    file_list3,
    function(x) {
      relocate(x, any_of(overlaps))
    }
  )

  file_list4[[1]]$markers <- markers1
  file_list4[[2]]$markers <- markers2
  merged <- file_list4 %>% purrr::reduce(full_join, by = "markers")
  ID <- merged$markers
  merged <- clean_input_data(merged)

  if (phased == TRUE) {
    message("Assuming the data are haplotypes.")
    merged <- merged %>% dplyr::mutate(across(tidyselect::everything(), ~ case_when(
      . == "A" ~ "A/A",
      . == "T" ~ "T/T",
      . == "C" ~ "C/C",
      . == "G" ~ "G/G",
      TRUE ~ .x
    )))
  } else if (phased == FALSE) {
    merged <- merged %>% dplyr::mutate(across(tidyselect::everything(), ~ case_when(
      . == "A" ~ "A/A",
      . == "T" ~ "T/T",
      . == "C" ~ "C/C",
      . == "G" ~ "G/G",
      . == "T/C" ~ "C/T",
      . == "T/G" ~ "G/T",
      . == "T/A" ~ "A/T",
      . == "G/A" ~ "A/G",
      . == "C/G" ~ "G/C",
      . == "C/A" ~ "A/C",
      TRUE ~ .x
    )))
  } else {
    stop("Parameter haplotype is required.")
  }

  overlap1 <- merged %>% dplyr::select(dplyr::ends_with(".x"))
  overlap2 <- merged %>% dplyr::select(dplyr::ends_with(".y"))
  for_conc <- QurvE::zipFastener(overlap1, overlap2, along = 2)
  for_conc2 <- data.frame(ID, for_conc)
  for_conc2[is.na(for_conc2)] <- "N"
  names(for_conc2) <- sub("^X", "", names(for_conc2))

  return(for_conc2)
}


#' Perform concordance analysis
#' @description
#' Compares genetic data information generated by different sequencing techniques.
#' It compares the genetic data of the same individual from technique 1 and technique 2.
#'
#' @param dataframe The dataframe containing genotype information of samples sequenced using two different techniques.
#'
#' @returns A list containing the dataframe of the tally of concordant and discordant calls and the concordance plot.
plot_concordance <- function(dataframe) {
  dataframe <- dplyr::rename(dataframe, ID = 1)
  x_cols <- dataframe %>% dplyr::select(dplyr::ends_with(".x"))
  y_cols <- dataframe %>% dplyr::select(dplyr::ends_with(".y"))
  Total <- ncol(x_cols)
  Incomparable <- rowSums(x_cols == "N" | y_cols == "N")
  Concordant <- rowSums(x_cols == y_cols & x_cols != "N" & y_cols != "N")
  Discordant <- Total - (Concordant + Incomparable)

  concordance <- dataframe %>%
    dplyr::mutate(
      Total = Total,
      Incomparable = Incomparable,
      Concordant = Concordant,
      Discordant = Discordant
    ) %>%
    dplyr::relocate(Total, Incomparable, Concordant, Discordant, .after = ID)

  ID <- concordance[, 1]
  pivot2 <- concordance[, 3:5]
  pivot <- data.frame(ID, pivot2)
  pivot <- pivot %>%
    tidyr::pivot_longer(!ID,
      names_to = "Condition",
      values_to = "Count"
    )

  Count <- as.integer(pivot$Count)
  rsID <- pivot$ID
  Condition <- pivot$Condition
  visual <- data.frame(rsID, Count, Condition)

  plot_conc <- visual %>%
    arrange(Count) %>%
    arrange(Condition) %>%
    mutate(rsID = forcats::fct_inorder(rsID)) %>%
    ggplot(aes(fill = Condition, x = rsID, y = Count)) +
    geom_bar(position = "stack", stat = "identity") +
    theme(
      axis.text.x = element_text(
        angle = 90,
        vjust = .3
      ),
      panel.background = element_blank()
    ) +
    scale_fill_manual(values = c(
      "#1ca7ec",
      "#fb7a8e",
      "#1f2f98"
    ))

  return(list(
    results = concordance,
    plot = plot_conc
  ))
}
