#' Homogenize data format
#'
#' @param input The input as a dataframe.
#'
#' @returns A dataframe containing information from the input file.
clean_input_data <- function(file) {
  # Get first value of the third column, assumption is genotype
  val <- file[1, 3]
  # check if it contains two letters
  is_a_char <- stringr::str_count(val, "[A-Za-z]") == 2

  if (isFALSE(is_a_char)) {
    stop("Unsupported file format: Genotype should be present by the third column. Only biallelic markers are accepted.")
  }

  if (isTRUE(grepl("/", val))) {
    file <- as.data.frame(file)
    file[is.na(file)] <- "N"
  } else if (isTRUE(grepl("|", val))) {
    file <- lapply(file, function(x) gsub("|", "/", x, fixed = TRUE))
    file <- as.data.frame(file)
    file[is.na(file)] <- "N"
  } else {
    file[, 3:ncol(file)] <- lapply(file[, 3:ncol(file), drop = FALSE], function(x) {
      x <- as.character(x)
      x <- sub("^([A-Za-z])([A-Za-z])$", "\\1/\\2", x)
      x
    })
    file[is.na(file)] <- "N"
  }

  file <- file %>%
    mutate(across(everything(), as.character)) %>%
    mutate(across(everything(), ~ case_when(
      .x %in% c("N/A", "NA", ".") ~ "N",
      TRUE ~ .x
    )))
  file <- as.data.frame(file)

  return(file)
}

#' Convert dataframes to genind object
#'
#' @param file A dataframe containing sample, population, and genotype information.
#' @param to_str To indicate whether input file will be used to generate .str files. Default is FALSE.
#' @param popinfo To indicate if population metadata is present. Default is TRUE.
#'
#' @returns A genind object.
convert_to_genind <- function(file, to_str = FALSE, popinfo = TRUE) {
  if (isFALSE(popinfo) && isTRUE(to_str)) {
    stop("Population data are required when converting to STRUCTURE format.")
  }

  file <- clean_input_data(file)

  if (isTRUE(popinfo)) {
    file <- dplyr::rename(file, Ind = 1, Pop = 2)
  } else {
    file <- dplyr::rename(file, Ind = 1)
  }

  if (isTRUE(to_str)) {
    # For Plotting
    populations_df <- file[, 1:2]
    colnames(populations_df) <- c("Label", "Population")
    populations_df$Label <- rownames(populations_df)

    ### Change pop to numeric - for STRUCTURE
    pop_df <- as.data.frame(file$Pop) %>%
      dplyr::rename(pops = 1)

    # Get total no. of pops
    pop_df_unique <- as.data.frame(pop_df[!duplicated(pop_df), ]) %>%
      dplyr::rename(pops = 1)
    # add row names as numbers
    pop_df_unique$num <- seq_len(nrow(pop_df_unique))

    # replace the pops in the original df (pop_df) with the numbers
    pop_df_corr <- dplyr::left_join(pop_df, pop_df_unique, by = "pops")
    pops <- pop_df_corr$num

    # Change Ind to numeric
    ind_only <- as.data.frame(file[, 1])
    ind_only$num <- rownames(ind_only)
    
    # New file
    Label <- file$Ind
    STR_label <- ind_only$num
    Pop <- file$Pop
    STR_pop <- pop_df_corr$num
    csv_revised <- data.frame(Label, STR_label, Pop, STR_pop, file[,-c(1,2)])

    ind <- as.character(ind_only$num)
    pop <- as.character(pops)
    fsnps_geno <- file[, 3:ncol(file)]
  } else {
    if (isFALSE(popinfo)) {
      ind <- as.character(file$Ind)
      pop <- NULL
      fsnps_geno <- file[, 2:ncol(file)]
    } else {
      ind <- as.character(file$Ind)
      pop <- as.character(file$Pop)
      fsnps_geno <- file[, 3:ncol(file)]
    }
  }

  fsnps_gen <- adegenet::df2genind(fsnps_geno,
    ind.names = ind,
    pop = pop,
    sep = "/",
    NA.char = "N",
    ploidy = 2,
    type = "codom"
  )

  if (isTRUE(to_str)) {
    fsnps_gen@pop <- as.factor(pop)
    return(list(
      fsnps_gen = fsnps_gen,
      populations = pop_df_unique,
      csv_revised = csv_revised
    ))
  }

  if (isTRUE(popinfo)) {
    fsnps_gen@pop <- as.factor(file$Pop)
  }

  return(fsnps_gen)
}

#' Subset populations for PCA recalculation
#'
#' @param fsnps_gen The genind object containing sample, population, and genotype data.
#' @param pops The populations to be used for PCA recalculation.
#'
#' @return Genind object.
subset_genind_pop <- function(fsnps_gen, pops) {
  keep <- adegenet::pop(fsnps_gen) %in% pops
  fsnps_gen[keep, ]
}

#' Convert files to PLINK2.0 files
#'
#' @param input.file The filepath of the dataset (either VCF, VCF.GZ, BCF, or PLINK1.9 files).
#' @param original_name The filename of the input file, NULL if PLINK files.
#' @param isplink Indicates whether input file/s are PLINK files. Default is FALSE.
#' @param name The prefix of the output file.
#' @param output_chr The chromosome formatting for output files based on PLINK2.0. Default is '26' indicating numeric codes.
#'
#' @returns The prefix of PLINK2.0 files.
convert_to_plink2 <- function(input.file,
                              original_name = NULL,
                              isplink = FALSE,
                              name = "convert_to_plink2",
                              output_chr = "26") {
  plink_path <- get_plink2_path()
  check_name <- if (is.null(original_name)) {
    input.file
  } else {
    original_name
  }

  is_vcf <- grepl("\\.vcf(\\.gz)?$", check_name, ignore.case = TRUE)
  is_bcf <- grepl("\\.bcf$", check_name, ignore.case = TRUE)

  if (isplink == TRUE) {
    args <- c(
      "--bfile", input.file,
      "--make-pgen",
      "--sort-vars",
      "--output-chr", output_chr,
      "--out", name
    )
  } else if (is_vcf || is_bcf) {
    input_flag <- if (is_bcf) "--bcf" else "--vcf"
    args <- c(
      input_flag, input.file,
      "--make-pgen",
      "--sort-vars",
      "--output-chr", output_chr,
      "--out", name
    )
  } else {
    stop("Unsupported file type. Please provide a VCF, VCF.GZ, or BCF.")
  }

  res <- system2(plink_path, args = args, stdout = TRUE, stderr = TRUE)
  return(name)
}

#' Convert files to PLINK1.9 files
#'
#' @param input.file The filepath of the dataset (either VCF, VCF.GZ, BCF, or PLINK1.9 files).
#' @param output.dir The filename of the input file, NULL if PLINK files.
#' @param isplink Indicates whether input file/s are PLINK files. Default is FALSE.
#' @param name The prefix of the output file.
#' @param output_chr The chromosome formatting for output files based on PLINK2.0. Default is '26' indicating numeric codes.
#'
#' @returns The prefix of PLINK2.0 files.
convert_to_plink <- function(input.file, name = "converted_to_plink") {
  plink_path <- get_plink_path()

  file_extension <- tools::file_ext(input.file)

  if (file_extension == "bcf") {
    system2(plink_path, args = c(
      "--bcf", input.file,
      "--make-bed",
      "--out", name
    ))
  } else if (file_extension %in% c("vcf", "vcf.gz", ".gz")) {
    system2(plink_path, args = c(
      "--vcf", input.file,
      "--make-bed",
      "--out", name
    ))
  } else {
    stop("Unsupported file type. Please provide a VCF, VCF.GZ, or BCF.")
  }

  return(name)
}


#' Convert zipped files to PLINK 2.0 files
#'
#' @param input_file The filepath of the zipped file.
#' @param output.dir The directory to save the unpacked files. Default is the working directory.
#'
#' @returns The prefix of the merged dataset.
#prepare_input_dataset_archive <- function(df_list, ext, output.dir = ".") {
convert_merge_to_plink <- function(df_list, plink_files = FALSE, output.dir = ".") {
  
  work_dir <- file.path(output.dir, "for_processing")
  dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
  merge_list_path <- file.path(work_dir, "merge_list.txt")
  prefixes <- c()

  if (isFALSE(plink_files)){  
    for (f in df_list) {
      base <- sub("\\.vcf(\\.gz)?$|\\.bcf$", "", basename(f), ignore.case = TRUE)
      out_pref <- file.path(work_dir, paste0(base, "_p2"))
      converted <- convert_to_plink(
        input.file = f,
        name = out_pref
      )
      prefixes <- c(prefixes, converted)
    }
  } else {
      prefixes <- df_list
    }

  if (length(prefixes) == 0) {
    stop("No valid files in the directory.")
  }

  writeLines(prefixes[-1], con = merge_list_path)
  merged_prefix <- file.path(work_dir, "merged_dataset")

  merge_plink2_files(
    base_prefix = prefixes[1],
    merge_list = merge_list_path,
    output_prefix = merged_prefix
  )

  return(list(pgen_prefix = merged_prefix))
}

#' Evaluate and convert data files to PLINK
#' @return Plink2 file prefix
convert_files <- function(input_file, output.dir = ".") {
  
  if (!file.exists(input_file)) {
    stop("File not found.")
  }
  
  # check if zipped or not
  file_extension <- tools::file_ext(input_file)
  if (file_extension %in% c("zip", "tar")) {
    unpacked <- unpack_input_file(input_file, output.dir)
    files <- unpacked$data_files
    ext <- tools::file_ext(files[[1]])
    
    if (ext %in% c("bed", "bim", "fam")) {
      unique_files <- unique(unlist(lapply(files, tools::file_path_sans_ext), use.names = FALSE))
      unique_files <- as.list(unique_files)
      merged_plink <- convert_merge_to_plink(unique_files, output.dir = output.dir, plink_files = TRUE)
    } else {
      merged_plink <- convert_merge_to_plink(files, output.dir = output.dir, plink_files = FALSE)
    }
    return(merged_plink$pgen_prefix)
    
  } else if (file_extension == "bed") {
    input_prefix <- tools::file_path_sans_ext(input_file)
    pgen_prefix <- file.path(output.dir, "input_pgen")
    convert_to_plink2(input_file, original_name = input_prefix, isplink = TRUE, name = pgen_prefix, output_chr = "26")
    return(pgen_prefix)
    
  } else {
    pgen_prefix <- file.path(output.dir, "input_pgen")
    convert_to_plink2(input_file, original_name = NULL, isplink = FALSE, name = pgen_prefix, output_chr = "26")
    return(pgen_prefix)
    
  }
}

#' Convert PLINK 1.9 files to other formats
#'
#' @param prefix The prefix of the PLINK files.
#' @param output_type The desired output file type.
#' @param output.dir The directory to save file. Default is the working directory.
#' @param ref The population metadata of samples. Required only if output_type is CSV.
#'
#' @returns The file path of the output.
convert_from_plink2 <- function(prefix,
                                output_type,
                                output.dir = ".",
                                ref = NULL) {
  plink2_path <- get_plink2_path()
  out_prefix <- file.path(output.dir, "converted")

  if (output_type == "vcf2") {
    system2(plink2_path, args = c(
      "--pfile", prefix,
      "--recode", "vcf",
      "--out", out_prefix
    ))
    return(paste0(out_prefix, ".vcf"))
  }

  if (output_type == "plink2") {
    files <- list.files(
      path = dirname(out_prefix),
      pattern = paste0("^", basename(out_prefix), "\\."),
      full.names = TRUE
    )
    zip_path <- paste0(out_prefix, ".zip")
    zip::zipr(zipfile = zip_path, files = files)

    return(zip_path)
  }

  if (output_type == "plink1") {
    system2(plink2_path, args = c(
      "--pfile", prefix,
      "--make-bed",
      "--out", out_prefix
    ))

    files <- list.files(
      path = dirname(out_prefix),
      pattern = paste0("^", basename(out_prefix), "\\."),
      full.names = TRUE
    )

    zip_path <- paste0(out_prefix, "_plink.zip")
    zip::zipr(zipfile = zip_path, files = files)
    return(zip_path)
  }

  if (output_type == "csv2") {
    vcf_file <- convert_from_plink2(prefix, "vcf2", output.dir)
    csv <- vcf_to_csv(
      vcf_file,
      ref = ref,
      output.dir = output.dir
    )
    return(csv)
  }
}

#' Merge PLINK files
#'
#' @param merge_list The list of PLINK 2.0 files.
#' @param output_prefix The prefix of the output (expected to be PLINK2.0).
#'
#' @returns The file path of the output.
merge_plink2_files <- function(base_prefix, merge_list, output_prefix) {
  plink_path <- get_plink_path()
  plink2_path <- get_plink2_path()

  args <- c(
    "--bfile", base_prefix,
    "--merge-list", merge_list,
    "--make-bed",
    "--out", output_prefix
  )
  system2(plink_path, args = args)

  args2 <- c(
    "--bfile", output_prefix,
    "--make-pfile",
    "--out",
    output_prefix
  )
  system2(plink2_path, args = args2)

  return(output_prefix)
}


#' Convert VCF to CSV
#'
#' @param files The file path of input file.
#' @param ref (optional) The dataframe containing metadata of samples.
#' @param output.dir (optional) The output directory to save input files if zipped. Default is working directory.
#'
#' @returns The dataframe of merged metadata with genotype information.
vcf_to_csv <- function(files, ref = NULL, output.dir = ".") {
  extension <- tools::file_ext(files)

  if (extension %in% c("vcf", "gz")) {
    raw_file <- load_vcf_files(files, output.dir = output.dir)
    raw_file <- dplyr::rename(raw_file, Sample = 1)
    #raw_file$Sample <- gsub("HGDP([0-9]+)_HGDP\\1$", "HG\\1", raw_file$Sample)
    raw_file$Sample <- trimws(raw_file$Sample)
  } else {
    stop("Input is not a VCF file.")
  }

  final_df <- raw_file

  if (is.null(ref)) {
     
    return(list(
       with_meta = as.data.frame(final_df),
       missing = NULL
    ))
  } else {
    merged <- add_metadata(final_df, ref)
    return(merged)
  }
}


#' Add metadata to dataframes
#' 
#' @param df The dataframe containing sample and genotype information.
#' @param metadata The dataframe containing the sample and population information. Sample name should match the samples with genotype information.
add_metadata <- function(df, metadata) {
  # check if single population
  if (is.character(metadata) && length(metadata) == 1) {
    df$Pop <- metadata
    df <- df %>% dplyr::relocate(tidyselect::last_col(), .after = 1)
    return(list(
      with_meta = df,
      missing = NULL
    ))
  } else {
    ref_data <- data.frame(metadata)
    ref_data <- dplyr::rename(ref_data, Sample = 1)
    final_df <- dplyr::rename(df, Sample = 1)
    
    samples_vcf <- final_df$Sample
    samples_ref <- ref_data$Sample
    
    if (sum(samples_vcf %in% samples_ref) == 0) {
      stop("Sample IDS do not match between VCF and metadata.")
    }
    
    cols <- colnames(ref_data)
    
    with_meta_data <- final_df %>%
      dplyr::inner_join(ref_data, by = "Sample") %>%
      relocate(cols, .after = 1) 
    
    missing_meta <- final_df %>% dplyr::anti_join(ref_data, by = "Sample")
    
    return(list(
      with_meta = with_meta_data,
      missing = missing_meta
    ))
  }
}

#' Convert SNP genotypes to dosages
#'
#' @param df The dataframe of SNP calls.
#' @param markers The reference file containing metadata of samples with information on 1. rsID 2. REF 3. ALT.
#'
#' @returns The dataframe of dosages per marker and samples.
to_binary <- function(df, markers = marker.file) {
  df <- as.data.frame(df)
  rownames(df) <- paste(df[, 1], "id", sep = "_")
  data <- df[, -c(1, 2)]
  data <- data.frame(t(data))
  row_name <- data.frame(rownames(data))
  data <- lapply(data, function(x) gsub("/", "", x, fixed = TRUE))
  data <- as.data.frame(data)
  revised_data <- data.frame(row_name, data)
  revised_data <- dplyr::rename(revised_data, ID = 1)

  marker <- as.data.frame(markers)
  marker <- marker[, -c(2, 3, 4)]
  marker <- marker[, 1:3]
  marker <- marker %>% mutate(across(tidyselect::everything(), ~ case_when(
    . == "A" ~ "AA",
    . == "T" ~ "TT",
    . == "C" ~ "CC",
    . == "G" ~ "GG",
    TRUE ~ .x
  )))
  marker <- dplyr::rename(marker, ID = 1, REF = 2, ALT = 3)

  df_marker <- merge(marker, revised_data, by = "ID", all.x = TRUE)
  df_marker <- df_marker %>%
    mutate(across(
      ends_with("_id"),
      ~ case_when(
        .x == df_marker$REF ~ 0,
        .x == df_marker$ALT ~ 2
      )
    ))

  df_marker[is.na(df_marker)] <- 1
  final_df <- df_marker[, -c(1, 2, 3)]
  final_df <- data.frame(t(final_df))
  return(final_df)
}


#' Convert file to gen_tibble object
#'
#' @param file The dataframe of SNP calls.
#' @param loci.meta The file path to the reference file containing metadata of markers with rsID, chr, pos, genetic distance, ref allele, and alt allele information.
#'
#' @returns A gentibble object.
csv_to_gentibble <- function(file, loci.meta = loci.meta) {
  df <- load_csv_xlsx_files(file)
  meta <- df[, 1:2]
  meta <- dplyr::rename(meta, id = 1, population = 2)

  loci_meta <- load_csv_xlsx_files(loci.meta)
  loci_meta <- dplyr::rename(loci_meta,
    name = 1,
    chromosome = 2,
    position = 3,
    genetic_dist = 4,
    allele_ref = 5,
    allele_alt = 6
  )

  geno <- to_binary(df, markers = loci_meta)
  geno <- as.matrix(geno)

  gentibble <- tidypopgen::gen_tibble(
    x = geno,
    loci = loci_meta,
    indiv_meta = meta,
    valid_alleles = c("A", "T", "C", "G"),
    quiet = TRUE
  )

  return(gentibble)
}

#' Widen long genotype file
#'
#' @param files A zipped file containing long genotype files.
#' @param population (optional) The metadata of samples.
#' @param output.dir The directory to save unpacked files. Default is working directory.
#'
#' @returns The dataframe of widened and merged genotype files.
widen_genotype_file <- function(files = files,
                                #population = NULL,
                                output.dir = ".") {
  if (!file.exists(files)) {
    stop("File does not exist in the working directory")
  } else {
    files_raw <- unpack_input_file(files, output.dir = output.dir)
}

  data_list <- files_raw$data_files
  all.list <- list()

  # check file extension
  f_ext <- tools::file_ext(data_list[1])

  if (f_ext == "xlsx") {
    for (x in data_list) {
      all.list[[x]] <- readxl::read_excel(x,
        sheet = 1,
        col_names = TRUE
      )
    }
  } else if (f_ext == "csv") {
    for (x in data_list) {
      all.list[[x]] <- read.csv(x, check.names = FALSE, row.names = 1)
    }
  } else {
    stop("Zippes files in unsupported format. Ensure they are CSV/XLSX files.")
  }

  new_colnames <- c("Sample", "ID", "Allele")
  dflist_new <- lapply(all.list, setNames, new_colnames)

  dflist_corrected <- lapply(
    dflist_new,
    function(x) {
      stats::aggregate(Allele ~ Sample + ID, x, paste, collapse = "/")
    }
  )

  df_list <- purrr::map(dflist_corrected, ~ (tidyr::pivot_wider(.x,
    names_from = Sample,
    values_from = Allele
  )))

  merged <- df_list %>% purrr::reduce(full_join, by = "ID")
  id <- merged$ID
  correct_alleles <- merged
  correct_alleles <- data.frame(t(correct_alleles))
  names(correct_alleles) <- correct_alleles[1, ]
  correct_alleles <- correct_alleles[-1, ]

  samples <- data.frame(colnames(merged[, -1]))
  corrected <- correct_alleles

  Samples <- rownames(corrected)
  corrected <- data.frame(Samples, corrected)
  rownames(corrected) <- NULL
  
  return(corrected)

#  if (is.null(population)) {
#    return(corrected)
#  } else {
#    pop_data <- load_csv_xlsx_files(population)
#    pop_data <- dplyr::rename(pop_data, Sample = 1, Population = 2)
#    matched <- corrected %>% dplyr::left_join(pop_data, by = "Sample")
#    data_length <- as.integer(ncol(corrected) - 1)
#    data_matched <- matched[, 2:data_length]
#    meta_begin <- as.integer(ncol(corrected) + 1)
#    meta <- matched[, meta_begin:ncol(matched)]

#    final_df <- dplyr::bind_cols(matched$Sample, meta, data_matched)
#    names(final_df)[names(final_df) == "matched$Sample"] <- "Sample"
#    return(final_df)
#  }
}

#' Convert dataframe to structure file
#'
#' @param file The dataframe containing sample, population, and genotype information.
#' @param to_str whether input file will be used to generate .str files.
#'
#' @returns The file path to the .str file
#'
#' @seealso [convert_to_genind()] to convert sample, population, and genotype file to STRUCTURE v2.3.4-compatible format. Set 'to_str' to TRUE.
revise_structure_file <- function(file, output.dir = ".", system = "Windows") {
  fsnps_gen_sub <- poppr::popsub(file)
  path <- file.path(output.dir, "structure_file.str")

  if (system == "Windows") {
    genind2structure2(fsnps_gen_sub, file = path, pops = TRUE, unix = FALSE)
  } else if (system == "Linux") {
    genind2structure2(fsnps_gen_sub, file = path, pops = TRUE, markers = TRUE, unix = TRUE)
    system(paste("tr '\t' ' '", shQuote(path), ">", shQuote(path)))
    system(paste("sed -e 's/ /\t/2' -e 's/ /\t/1'", shQuote(path), ">", shQuote(path)))
  }
  return(path)
}


#' Convert dataframe to SNIPPER-analysis-compatible file
#'
#' @param input The dataframe with genotype and sample information.
#' @param references The dataframe containing population metadata of input. Sample names should match.
#' @param target.pop This indicates whether test data will be used. Otherwise, all data will be used included as a training set. Default is TRUE.
#' @param population.name The population to be used if test data will be used.
#' @param markers The number of markers in the dataset.
#'
#' @returns A dataframe formatted for SNIPPER compatibility.
to_snipper <- function(input,
                       references,
                       target.pop = TRUE,
                       population.name = NULL,
                       markers = snps) {
  if (is.data.frame(input)) {
    input.file <- input
  } else {
    stop("Not a dataframe.")
  }

  tosnipper <- lapply(
    input.file,
    function(x) {
      gsub(pattern = "/", replacement = "", x = x, fixed = TRUE)
    }
  )

  tosnipper <- as.data.frame(tosnipper)

  if (class(tosnipper$Sample) != "character") {
    tosnipper$Sample <- as.character(tosnipper$Sample)
  }

  if (is.data.frame(references)) {
    reference <- references
  } else {
    stop("Not a dataframe.")
  }

  if (class(reference$Sample) != "character") {
    reference$Sample <- as.character(reference$Sample)
  }

  matched <- tosnipper %>% dplyr::left_join(reference, by = "Sample")
  last.col <- as.integer(ncol(matched))
  sec.last <- last.col - 1
  Superpop <- as.data.frame(matched[, last.col])
  Population <- as.data.frame(matched[, sec.last])
  Sample <- as.data.frame(matched$Sample)
  data <- as.data.frame(matched[, 2:ncol(matched) - 1])
  drops <- "Sample"
  data <- data[, !(names(data) %in% drops)]

  to_excel <- dplyr::bind_cols(Population, Superpop, Sample, data)
  to_excel <- dplyr::rename(to_excel, Population = 1, Superpop = 2, Sample = 3)

  tosnpr_split <- split(to_excel, to_excel$Population)
  tosnpr_split <- tosnpr_split %>% map(`rownames<-`, NULL)
  tosnpr_split <- lapply(
    tosnpr_split,
    function(x) {
      x$No <- rownames(x)
      as.data.frame(x)
      data <- as.data.frame(x[, 2:ncol(x) - 2])
      Sample <- x[, 1]
      x <- dplyr::bind_cols(x$No, Sample, data)
    }
  )

  merged <- plyr::ldply(tosnpr_split, data.frame)
  merged <- merged[, -c(1, 3)]

  if (target.pop == TRUE) {
    target <- merged[merged$Population == population.name, ]
    target$snpr <- "0"
    non.target <- merged[merged$Population != population.name, ]
    non.target$snpr <- "1"
    merged2 <- dplyr::bind_rows(target, non.target)
  } else if (target.pop == FALSE) {
    merged2 <- merged
    merged2$snpr <- "1"
  } else {
    stop("Parameter target.pop is not stated.")
  }

  sample.count <- as.integer(nrow(merged2))
  pop.only <- as.data.frame(merged2$Superpop)
  pop.only <- pop.only[!duplicated(pop.only), ]
  pop.only <- as.data.frame(pop.only)
  pop.count <- as.integer(nrow(pop.only))
  merged3 <- merged2[, -2]
  merged3 <- as.data.frame(merged3)

  names(merged3)[names(merged3) == "...1"] <- sample.count
  names(merged3)[names(merged3) == "Superpop"] <- markers
  names(merged3)[names(merged3) == "Sample"] <- pop.count
  names(merged3)[names(merged3) == "snpr"] <- ""

  merged3 <- rbind(NA, merged3)
  merged3 <- rbind(NA, merged3)
  merged3 <- rbind(NA, merged3)
  merged3 <- rbind(NA, merged3)

  return(merged3)
}

#' Convert alignment to DNA bin
#'
#' @param path The alignment object.
#'
#' @returns DNA bin.
alignment_to_dnabin <- function(path) {
  ext <- tolower(tools::file_ext(path))

  if (ext %in% c("fasta", "fa", "fas")) {
    fmt <- "fasta"
  } else if (ext == "aln") {
    fmt <- "clustal"
  } else if (ext == "msf") {
    fmt <- "msf"
  } else if (ext == "mase") {
    fmt <- "mase"
  } else if (ext == "phylip") {
    fmt <- "phylip"
  } else {
    stop("Unsupported file format: ", ext)
  }

  alignment <- seqinr::read.alignment(
    path,
    format = fmt
  )
  mat <- do.call(rbind, lapply(alignment$seq, function(x) {
    strsplit(toupper(x), "")[[1]]
  }))
  rownames(mat) <- alignment$nam
  return(ape::as.DNAbin(mat))
}

#' Generate Arlequin-compatible file (.arp)
#'
#' @param df The dataframe containing marker and population information.
#' @param genotypic_data To indicate if the data to be analyzed is haplotypic or genotypic.
#' @param recessive_data Indicates if genotypic data present are recessive.
#' @param locus_sep The character used to separate the alleles of different loci (locus separator).
#' @param data_type The type of data to be analyzed. Standard is DNA.
#' @param output.prefix The prefix of output (.arp) file.
#'
#' @returns The file path of the output .arp file.
build_arp_per_population <- function(df,
                                     genotypic_data = 1,
                                     gametic_phase = 0,
                                     recessive_data = 0,
                                     locus_sep = "WHITESPACE",
                                     output.prefix = "arp_file",
                                     data_type = "STANDARD") {
  workdir <- tempfile("arlequin_")
  dir.create(workdir)
  # workdir <- "."
  out_path <- file.path(workdir, paste0(output.prefix, ".arp"))

  if (file.exists(out_path)) {
    file.remove(out_path)
  }

  # get snps
  names(df) <- tolower(trimws(names(df)))
  names(df)[2] <- "population"

  # need to ensure no XML reserved characters sa population (<, >, or &)
  df$population <- gsub("&", "and", df$population)
  # df <- as.data.frame(df)

  subgroups <- unique(df$population)
  loci <- colnames(df)[-c(1, 2)]

  # Base Arlequin header
  arp <- c(
    "[Profile]",
    paste0('Title="Structure Analysis per Population"'),
    paste0("NbSamples=", length(subgroups)),
    paste0("NbLoci=", length(loci)),
    paste0("GenotypicData=", genotypic_data),
    paste0("GameticPhase=", gametic_phase),
    paste0("RecessiveData=", recessive_data),
    paste0("DataType=", data_type), # Dynamic: MICROSAT for LB, STANDARD for SB
    paste0("LocusSeparator=", locus_sep),
    'MissingData="?"',
    "",
    "[Data]",
    "[[Samples]]"
  )

  # Append each individual population as its own Sample section
  for (sub in subgroups) {
    sub_df <- df[df$population == sub, ]
    arp <- c(
      arp,
      paste0('SampleName="', sub, '"'),
      paste0("SampleSize=", nrow(sub_df)),
      "SampleData={"
    )

    for (i in 1:nrow(sub_df)) {
      allele1_row <- c()
      allele2_row <- c()

      for (loc in loci) {
        val <- sub_df[[loc]][i]

        if (is.na(val) || val %in% c("", "?")) {
          allele1_row <- c(allele1_row, "?")
          allele2_row <- c(allele2_row, "?")
          next
        }

        # Split alleles by slash
        alleles <- unlist(strsplit(as.character(val), "/"))
        alleles <- trimws(alleles)

        if (length(alleles) == 1) {
          # If only 1 allele is detected, female is homozygous. Duplicate the allele.
          allele1_row <- c(allele1_row, alleles[1])
          allele2_row <- c(allele2_row, alleles[1])
        } else {
          # Heterozygous
          allele1_row <- c(allele1_row, alleles[1])
          allele2_row <- c(allele2_row, alleles[2])
        }
      }

      # Build the two lines per sample (Arlequin format for GenotypicData=1)
      arp <- c(arp, paste(sub_df$sample[i], "1", paste(allele1_row, collapse = " ")))
      arp <- c(arp, paste("    ", paste(allele2_row, collapse = " ")))
    }
    arp <- c(arp, "}", "")
  }

  #   arp <- c(arp,
  #            "",
  #            "[[Structure]]",
  #            "",
  #            'StructureName="Population Structure"',
  #            paste0("NbGroups=", length(subgroups)))

  #   for (sub in subgroups) {
  #      arp <- c(
  #         arp,
  #         "Group={",
  #         paste0('"', sub, '"'),
  #         "}"
  #      )
  #   }

  writeLines(arp, out_path)
  message("SUCCESS: Generated Arlequin project file:", out_path)
  return(out_path)
}

#' Convert CSV file to gtype class (strataG-specific object)
#'
#' @param file The file path of the .csv/.xlsx file.
#'
#' @returns A dataframe in gtype format.
csv_to_gtype_format <- function(file) {
  df <- load_csv_xlsx_files(file)
  df <- clean_input_data(df)
  Sample <- df[, 1]
  Pop <- df[, 2]
  data <- df[3:ncol(df)]

  data <- data %>%
    mutate(across(everything(), as.character)) %>%
    mutate(across(everything(), ~ case_when(
      .x == "N" ~ "N/N",
      TRUE ~ .x
    )))

  data_clean <- data %>% tidyr::separate_wider_delim(
    cols = everything(),
    delim = "/",
    names_sep = "."
  )

  data_file <- dplyr::bind_cols(Sample, Pop, data_clean)

  return(data_file)
}
