#' Plot Depth of Coverage
#'
#' @param vcf The file path of the VCF input file.
#' @param output.dir The directory to save the plots. Default is the working directory.
#' @param reference The file path to the population data of the samples.
#' @param palette The color to differentiate the samples.
#' @param width The width of the png output.
#' @param height The height of the png output.
#' @param dpi The quality of the png output.
#'
#' @returns Depth plot.
depth_from_vcf <- function(vcf,
                           output.dir = ".",
                           reference,
                           palette = NULL,
                           width = 10,
                           height = 8,
                           dpi = 300) {
  vcf.file <- vcfR::read.vcfR(vcf)
  depth <- vcfR::extract.gt(vcf.file, element = "DP", as.numeric = TRUE)
  depth <- as.data.frame(t(depth))
  depth$Sample <- rownames(depth)

  depth_long <- depth %>% tidyr::pivot_longer(
    !Sample,
    names_to = "rsID",
    values_to = "Depth"
  )

  if (!is.null(reference)) {
    ref <- load_csv_xlsx_files(reference)
    ref <- dplyr::rename(ref, Sample = 1, highlight = 2)
    depth_long <- depth_long %>% dplyr::left_join(ref, by = "Sample")
    fill2 <- depth_long$highlight
  } else {
    fill2 <- NULL
  }

  p_rsid <- ggplot2::ggplot(depth_long, ggplot2::aes(x = rsID, y = Depth, fill = fill2)) +
    ggplot2::geom_boxplot() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = .4)) +
    ggplot2::scale_fill_brewer(palette = palette)

  out_dp_rsid <- file.path(output.dir, "Depth_marker.png")
  ggsave(out_dp_rsid, plot = p_rsid, width = width, height = height, dpi = dpi)

  # plot of depth per sample
  p_sample <- ggplot2::ggplot(depth_long, ggplot2::aes(x = Sample, y = Depth, fill = fill2)) +
    ggplot2::geom_boxplot() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = .4)) +
    ggplot2::scale_fill_brewer(palette = palette)

  out_dp_sample <- file.path(output.dir, "Depth_samples.png")
  ggsave(out_dp_sample, plot = p_sample, width = width, height = height, dpi = dpi)

  return(list(
    plot_marker = out_dp_rsid,
    plot_sample = out_dp_sample
  ))
}

#' Compute heterozygosities, MAR, and inbreeding
#'
#' @param fsnps_gen The genind data containing genotype and population information.
#'
#' @returns A list containing the metrices.
compute_pop_stats <- function(fsnps_gen) {
  mar_matrix <- hierfstat::allelic.richness(hierfstat::genind2hierfstat(fsnps_gen))$Ar %>%
    apply(MARGIN = 2, FUN = mean) %>%
    round(digits = 3)
  mar_list <- as.data.frame(mar_matrix)

  basic_fsnps <- hierfstat::basic.stats(fsnps_gen, diploid = TRUE)
  Ho <- apply(basic_fsnps$Ho, 2, mean, na.rm = TRUE) %>% round(2)
  He <- apply(basic_fsnps$Hs, 2, mean, na.rm = TRUE) %>% round(2)
  inv_het <- 1 / (1 - He) %>% round(3)
  heterozygosity_df <- data.frame(Population = names(Ho), Ho = Ho, He = He, Ae = inv_het) %>%
    tidyr::pivot_longer(cols = c("Ho", "He", "Ae"), names_to = "Variable", values_to = "Value")

  # T-Test between Heterozygosities
  ttest_df <- data.frame(
    Locus = rownames(basic_fsnps$perloc),
    basic_fsnps$perloc,
    row.names = NULL
  )

  # Inbreeding Coefficient
  fis_values <- apply(basic_fsnps$Fis, 2, mean, na.rm = TRUE) %>%
    round(3)
  fis_df <- data.frame(Population = names(fis_values), Fis = fis_values)

  return(list(
    mar_list = mar_list,
    heterozygosity = heterozygosity_df,
    ttest = ttest_df,
    inbreeding_coeff = fis_df
  ))
}

#' Compute heterozygosities, MAR, and inbreeding
#'
#' @inheritParams compute_pop_stats
#'
#' @returns An allele frequency dataframe.
compute_af <- function(fsnps_gen) {
  fsnps_gpop <- adegenet::genind2genpop(fsnps_gen)
  allele_freqs <- t(adegenet::makefreq(fsnps_gpop, quiet = FALSE, missing = NA)) %>%
    as.data.frame()

  allele_freqs <- data.frame(rownames(allele_freqs), allele_freqs)
  allele_freqs <- dplyr::rename(allele_freqs, markers = 1)
  rownames(allele_freqs) <- NULL

  return(allele_freqs)
}

#' Perform HWE test
#'
#' @inheritParams compute_pop_stats
#' @param correction The method of correction to be applied in calculating the alpha to identify out of HWE. Default is Bonferroni.
#' @param alpha The alpha value for the threshold. Default is 0.05.
#'
#' @returns A list of HWE stats results.
compute_hwe <- function(fsnps_gen, correction = "Bonferroni", alpha = 0.05, chains = 100000) {
  # Hardy-Weinberg Equilibrium (List for export)
  fsnps_hwe <- as.data.frame(round(pegas::hw.test(fsnps_gen, B = chains), 6))
  fsnps_hwe <- data.frame(rownames(fsnps_hwe), fsnps_hwe)
  fsnps_hwe <- dplyr::rename(fsnps_hwe, rsID = 1)
  rownames(fsnps_hwe) <- NULL

  # Chi-square test (Matrix for export, Data Frame for ggplot)
  fsnps_hwe_test <- data.frame(sapply(
    adegenet::seppop(fsnps_gen),
    function(ls) pegas::hw.test(ls, B = 0)[, 3]
  ))

  fsnps_hwe_chisq_matrix <- t(data.matrix(fsnps_hwe_test))
  fsnps_hwe_chisq_df <- as.data.frame(fsnps_hwe_chisq_matrix) %>% tibble::rownames_to_column("Population")

  # Monte Carlo: p value
  fsnps_hwe_mc <- data.frame(sapply(
    adegenet::seppop(fsnps_gen),
    function(ls) pegas::hw.test(ls, B = chains)[, 4]
  ))
  fsnps_hwe_mc <- t(data.matrix(fsnps_hwe_mc))
  fsnps_hwe_mc_df <- data.frame(t(fsnps_hwe_mc))
  fsnps_hwe_mc_df <- tibble::rownames_to_column(fsnps_hwe_mc_df, var = "rsID")

  loc_total <- adegenet::nLoc(fsnps_gen)
  alpha_value <- alpha / loc_total
  # loci out of HWE
  if (is.null(correction)) {
    loci_HWE_failure <- data.frame(
      rsID = colnames(fsnps_hwe_chisq_matrix),
      Chisq = apply(fsnps_hwe_chisq_matrix < alpha, 2, mean),
      MC = apply(fsnps_hwe_mc < alpha_value, 2, mean)
    )

    pops_out_of_HWE <- data.frame(
      Populations = rownames(fsnps_hwe_chisq_matrix),
      Chisq = apply(fsnps_hwe_chisq_matrix < alpha, 1, mean),
      MC = apply(fsnps_hwe_mc < alpha_value, 1, mean)
    )
  } else if (correction == "Bonferroni") {
    loci_HWE_failure <- data.frame(
      rsID = colnames(fsnps_hwe_chisq_matrix),
      Chisq = apply(fsnps_hwe_chisq_matrix < alpha_value, 2, mean),
      MC = apply(fsnps_hwe_mc < alpha_value, 2, mean)
    )

    pops_out_of_HWE <- data.frame(
      Populations = rownames(fsnps_hwe_chisq_matrix),
      Chisq = apply(fsnps_hwe_chisq_matrix < alpha_value, 1, mean),
      MC = apply(fsnps_hwe_mc < alpha_value, 1, mean)
    )
  } else if (correction == "FDR") {
    Chisq_FDR <- apply(fsnps_hwe_chisq_matrix, 1, p.adjust, method = "fdr")
    Chisq_FDR <- t(Chisq_FDR)

    MC_FDR <- apply(fsnps_hwe_mc, 1, p.adjust, method = "fdr")
    MC_FDR <- t(MC_FDR)

    loci_HWE_failure <- data.frame(
      rsID = colnames(fsnps_hwe_chisq_matrix),
      Chisq = apply(fsnps_hwe_chisq_matrix < alpha_value, 2, mean),
      MC = apply(fsnps_hwe_mc < alpha_value, 2, mean),
      Chisq_FDR = apply(Chisq_FDR < alpha_value, 2, mean),
      MC_FDR = apply(MC_FDR < alpha_value, 2, mean)
    )

    pops_out_of_HWE <- data.frame(
      Populations = rownames(fsnps_hwe_chisq_matrix),
      Chisq = apply(fsnps_hwe_chisq_matrix < alpha_value, 1, mean),
      MC = apply(fsnps_hwe_mc < alpha_value, 1, mean),
      Chisq_FDR = apply(Chisq_FDR < alpha_value, 1, mean),
      MC_FDR = apply(MC_FDR < alpha_value, 1, mean)
    )
  } else {
    stop("Correction method not supported")
  }

  rownames(loci_HWE_failure) <- NULL
  rownames(pops_out_of_HWE) <- NULL

  return(list(
    hw_summary = fsnps_hwe,
    hw_dataframe = fsnps_hwe_chisq_df,
    hw_mc = fsnps_hwe_mc_df,
    loci_HWE_failure = loci_HWE_failure,
    pops_out_of_HWE = pops_out_of_HWE
  ))
}

#' Compute FST
#'
#' @inheritParams compute_pop_stats
#'
#' @returns A list of FST stats results.
compute_fst <- function(fsnps_gen) {
  # Pairwise Fst matrix using Weir & Cockerham 1984 method
  fst_matrix_raw <- hierfstat::genet.dist(fsnps_gen, method = "WC84") %>%
    round(3)

  fst_list <- if (length(fst_matrix_raw) == 0) {
    list(message = "No Fst values calculated")
  } else {
    as.list(as.matrix(fst_matrix_raw))
  }

  fst_df <- as.data.frame(as.matrix(fst_matrix_raw)) %>%
    tibble::rownames_to_column("Site1") %>%
    tidyr::pivot_longer(cols = -Site1, names_to = "Site2", values_to = "Fst")

  return(list(
    fst_matrix = fst_list,
    fst_dataframe = fst_df
  ))
}

#' Plot heterozygosity tables
#'
#' @param Het_fsnps_df A dataframe containing the observed and expected heterozygosities ("Ho" and "He").
#' @param out_dir The directory to save the plot.
#'
#' @returns Plot of observed vs expected heterozygosity per population.
plot_heterozygosity <- function(Het_fsnps_df, out_dir) {
  out_path <- file.path(out_dir, "heterozygosity_plot.png")

  Het_fsnps_df <- Het_fsnps_df %>%
    dplyr::filter(Variable %in% c("Ho", "He"))

  Het_fsnps_df$Variable <- as.factor(Het_fsnps_df$Variable)

  p <- ggplot(Het_fsnps_df, aes(x = Population, y = Value, fill = Variable)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.6), colour = "black") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 0.5)) +
    scale_fill_manual(
      values = c("royalblue", "#bdbdbd"),
      labels = c(expression(italic("H")[o]), expression(italic("H")[e]))
    ) +
    labs(y = "Heterozygosity") +
    theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, face = "bold"))

  ggsave(out_path, plot = p, width = 9, dpi = 300)
  return(out_path)
}

#' Plot FST across populations
#'
#' @param Het_fsnps_df A dataframe containing the calculated FST values.
#' @param out_dir The directory to save the plot.
#'
#' @returns PNG file of the heatmap.
plot_fst <- function(fst_df, out_dir) {
  out_path <- file.path(out_dir, "fst_heatmap.png")

  p <- ggplot(fst_df, aes(x = Site1, y = Site2, fill = Fst, label = round(Fst, 3))) +
    geom_tile(color = "black") +
    geom_text(aes(label = round(Fst, 3)), size = 3, color = "black") +
    scale_fill_gradient2(
      low = "blue", mid = "pink", high = "red",
      midpoint = max(fst_df$Fst, na.rm = TRUE) / 2
    ) +
    labs(x = "Site 1", y = "Site 2", fill = "Fst") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  ggsave(out_path, plot = p, width = 8, height = 6, dpi = 300)
  return(out_path)
}


#' Compile and export calculated population summary statistics
#'
#' @param allele_freq The dataframe of calculated allele frequencies.
#' @param priv_alleles A list/string/dataframe of the total number of calculated private alleles per population.
#' @param stats_matrix A list of the calculated heterozygosities, MAR, and inbreeding coefficient.
#' @param hw_matrix A list of the HWE test results.
#' @param fst_matrix A list of the FST metrices.
#' @param dir The directory to save the result.
#'
#' @returns An excel (.xlsx) file containing the results in individual sheets.
export_pop_results <- function(allele_freq, priv_alleles, stats_matrix, hw_matrix, fst_matrix, dir = tempdir()) {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
  out_file <- file.path(dir, paste0("population-statistics-results_", timestamp, ".xlsx"))

  # Revise formats
  #------ MAR
  mar <- data.frame(rownames(stats_matrix$mar_list), stats_matrix$mar_list)
  mar <- dplyr::rename(mar, Pop = 1)
  rownames(mar) <- NULL

  #------ Heterozygosity
  het <- stats_matrix$heterozygosity
  het <- tidyr::pivot_wider(
    data = het,
    names_from = Variable,
    values_from = Value
  )

  #------- IC
  ic <- stats_matrix$inbreeding_coeff
  rownames(ic) <- NULL

  datasets <- list(
    "Private Alleles" = as.data.frame(priv_alleles),
    "Mean Allelic Richness" = mar,
    "Heterozygosities" = het,
    "T-test per Locus" = as.data.frame(stats_matrix$ttest),
    "Inbreeding Coefficient" = as.data.frame(stats_matrix$inbreeding_coeff),
    "Allele Frequencies" = allele_freq,
    "Hardy-Weinberg Equilibrium" = as.data.frame(hw_matrix$hw_summary),
    "Chi-square test HWE" = as.data.frame(hw_matrix$hw_dataframe),
    "Monte Carlo test HWE" = as.data.frame(hw_matrix$hw_mc),
    "Loci out of HWE" = as.data.frame(hw_matrix$loci_HWE_failure),
    "Pops out of HWE" = as.data.frame(hw_matrix$pops_out_of_HWE),
    "Pairwise Fst Matrix" = as.data.frame(fst_matrix$fst_matrix)
  )

  openxlsx::write.xlsx(datasets, file = out_file)
  return(out_file)
}


#' Evaluate if input file is an allele frequency or gt frequency table
#'
#' @param df The genotype input as a dataframe.
#' @param sample_size The total number of samples in the file.
#' @param genotype The format of genotype data.
#'
#' @returns A string indicating the data format.
evaluate_file <- function(df, sample_size = 50, genotype = "^[A-Z]/[A-Z]$") {
  all_vals <- unlist(df, use.names = FALSE)
  all_vals <- all_vals[!is.na(all_vals) & all_vals != "N"]
  sample_vals <- if (length(all_vals) > sample_size) {
    sample(all_vals, sample_size)
  } else {
    all_vals
  }

  num_gts <- sum(stringr::str_detect(sample_vals, genotype))
  num_vals <- suppressWarnings(as.numeric(sample_vals))
  num_freqs <- sum(!is.na(num_vals) & num_vals >= 0 & num_vals <= 1)

  if (num_gts > num_freqs) {
    return("gts")
  } else if (num_freqs > num_gts) {
    return("freqs")
  } else {
    return("unknown data format")
  }
}


#' Calculate the genotype frequencies per population
#'
#' @param df The dataframe containing sample, population, and genotype information.
#' @param pop The total number of samples in the data. Required if employing five-event minimum allele frequency.
#'
#' @returns A list containing genotype frequency (overall and by population).
calc_observed_genotype_freq <- function(df) {
  df <- dplyr::rename(df, pop = 2)
  pops_df <- split(df, df$pop)
  
  result <- lapply(names(pops_df), function(x) {
    pop_sub_df <- pops_df[[x]]
    pop_sub_df <- pop_sub_df[, -c(1,2), drop = FALSE]
    
    # Loop through each marker
    marker_results <- lapply(names(pop_sub_df), function(marker) {
      gt <- trimws(as.character(pop_sub_df[[marker]]))
      gt <- gt[!is.na(gt) & gt != "N" & gt != ""]
      
      if (length(gt) == 0) {
        return(NULL)
      }
      
      gt_counts <- table(gt)
      
      data.frame(
        population = x,
        marker = marker,
        genotype = names(gt_counts),
        count = as.integer(gt_counts),
        n_genotyped = length(gt),
        observed_freq = as.numeric(gt_counts)/length(gt),
        row.names = NULL
      )
    })
    bind_rows(marker_results)
  })
  names(result) <- names(pops_df)
  return(result)
}

calc_expected_genotype_freq <- function(df, pop = NULL) {
  df <- dplyr::rename(df, markers = 1)

  df <- df %>%
    mutate(
      marker = stringr::str_remove(markers, "\\..*"),
      allele = stringr::str_remove(markers, ".*\\.")
    ) %>%
    dplyr::select(-markers)

  df_long <- df %>% tidyr::pivot_longer(
    cols = -c(marker, allele),
    names_to = "population",
    values_to = "freq"
  )

  # double check single alleles
  geno_freqs1 <- df_long %>%
    dplyr::group_by(marker, population) %>%
    dplyr::summarise(
      n_alleles = dplyr::n(),
      allele1 = dplyr::first(allele),
      allele2 = dplyr::last(allele),
      p = dplyr::first(freq),
      q = dplyr::last(freq),
      .groups = "drop"
    )

  if (is.null(pop)) {
    geno_freqs2 <- geno_freqs1 %>%
      dplyr::mutate(
        q = dplyr::if_else(n_alleles == 1, 0, q),
        allele2 = dplyr::if_else(n_alleles == 1, NA_character_, allele2),
        homozygous1 = p^2,
        heterozygous = dplyr::if_else(n_alleles == 1, 0, 2 * p * q),
        homozygous2 = dplyr::if_else(n_alleles == 1, 0, q^2)
      )
  } else if (!is.null(pop)) {
    floor <- 5 / (2 * pop)

    geno_freqs2 <- geno_freqs1 %>%
      dplyr::mutate(
        p = pmax(p, floor),
        q = pmax(q, floor),
        homozygous1 = p^2,
        heterozygous = dplyr::if_else(n_alleles == 1, 0, 2 * p * q),
        homozygous2 = dplyr::if_else(n_alleles == 1, 0, q^2)
      )
  }
  return(geno_freqs2)

#  by_pop <- split(geno_freqs2, geno_freqs2$population)
#  by_pop <- lapply(by_pop, function(x) {
#    x <- x[, -2]
#  })

#  clean_names <- names(by_pop) %>%
#    stringr::str_replace_all("[._]", " ")
#  names(by_pop) <- clean_names

#  return(list(
#    gt_complete = geno_freqs2,
#    gt_by_pop = by_pop
#  ))
}

standardize_names <- function(x) {
  x <- gsub("[^[:alnum:]]+", " ", x)
  x <- trimws(x)
  x <- tolower(x)
  x
}

af_to_long <- function(df) {
  df %>%
    dplyr::rename(markers = 1) %>%
    dplyr::mutate(
      allele = sub("^.*\\.", "", markers),
      marker = sub("\\.[^.]*$", "", markers)
    ) %>%
    dplyr::select(-markers) %>%
    tidyr::pivot_longer(
      cols = -c(marker, allele),
      names_to = "population",
      values_to = "freq"
    ) %>%
    dplyr::mutate(
      marker = standardize_names(marker),
      allele = trimws(allele),
      freq = as.numeric(freq)
    )
}

#' Calculate forensic parameters for iisnps
#'
#' @param geno_freqs The dataframe of genotype frequencies per population.
#' @param profile (optional) If target profile will be provided for identification.
#' @param theta The coefficient of ancestry or inbreeding to adjust calculation of match probabilities. Required if profile is not NULL.
#'
#' @returns A list of dataframes containing results.
calc_iisnps_params <- function(geno_freqs, af) {
  
  res <- lapply(geno_freqs, function(x) {
    # split by marker
    pops_marker <- split(x, x$marker)
    
    # calculate Gsqr
    marker_res <- lapply(names(pops_marker), function(y) {
      curr <- pops_marker[[y]]
      Gsqr <- sum(curr$observed_freq^2, na.rm = TRUE)
      alleles <- stringr::str_split_fixed(curr$genotype, "/", 2)
      
      curr <- curr %>% dplyr::mutate(
        gt_type = ifelse(alleles[,1] == alleles[,2], "Homozygous", "Heterozygous")
      )
      
      total <- unique(na.omit(curr$n_genotyped))[1]
      
      gt_res <- curr %>%
        dplyr::group_by(gt_type) %>%
        dplyr::summarise(
          Count = sum(count, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        dplyr::mutate(Frequency = Count/total) %>%
        dplyr::select(gt_type, Frequency) %>%
        tidyr::pivot_wider(
          names_from = gt_type,
          values_from = Frequency,
          values_fill = 0
        )
      
      gt_res <-  gt_res %>% dplyr::mutate(Gsqr = Gsqr)
      data.frame(y, gt_res) %>%
        dplyr::rename(marker = 1)
      
    })
    return(marker_res)
  })
  
  marker_metrics <- lapply(res, function(x) {
    x <- x %>%
      purrr::list_flatten() %>%
      purrr::list_flatten() %>%
      bind_rows()

    x <- dplyr::rename(x, RMP = 4)
    x <- x %>% dplyr::mutate(
      PD = 1 - RMP,
      PE = (Heterozygous^2)*(1 - (2*(Heterozygous)*(Homozygous^2))),
      TPI = 1/(2*Homozygous)
      )
  })
  
  pic <- af %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      PIC = 1 - (homozygous1 + homozygous2) - (2 * (homozygous1 * homozygous2)), # 2 * homozygous1 * homozygous2 * (1-2*homozygous1*homozygous2),
    ) %>%
    dplyr::ungroup()
  
  pic_pop <- split(pic, pic$population)
  
  # match 
  metrics1 <- standardize_names(names(marker_metrics))
  metrics2 <- standardize_names(names(pic_pop))
  match_df <- match(metrics1, metrics2)
  
  metrics_updated <- lapply(seq_along(marker_metrics), function(x){
    if (is.na(match_df[x])) {
      return(marker_metrics[[x]])
    }
    
    df1 <- marker_metrics[[x]]
    df2 <- pic_pop[[match_df[x]]]
    
    df1 <- df1 %>% left_join(
      df2 %>% dplyr::select(marker, PIC),
      by = "marker"
    )
    return(df1)
  })
  names(metrics_updated) <- names(marker_metrics)
  return(metrics_updated)
}

calc_rmp <- function(profile, af_table, n_ref, theta = 0.01) {
  profile <- profile %>% dplyr::mutate(
    marker_std = standardize_names(marker),
    genotype_std = clean_input_data(genotype)
  )
  freq_floor <- 5/(2*n_ref)
  af_table <- af_table %>% 
    dplyr::group_by(marker) %>%
    dplyr::mutate(
      freq_adj = pmax(freq, freq_floor),
      freq_adj = freq_adj / sum(freq_adj)
    ) %>%
    dplyr::ungroup
  
  locus_res <- profile %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      allele1 = strsplit(genotype, "/", fixed = TRUE)[[1]][1],
      allele2 = strsplit(genotype, "/", fixed = TRUE)[[1]][2],
    ) %>%
    dplyr::ungroup() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      
      p = af_table$freq_adj[
        match(paste(marker, allele1), paste(af_table$marker, af_table$allele))
      ],
      
      q = af_table$freq_adj[
        match(paste(marker, allele2), paste(af_table$marker, af_table$allele))
      ],
      
      gen_prob = dplyr::case_when(
        is.na(p) | is.na(q) ~ NA_real_,
        allele1 == allele2 ~ p^2 + theta*p*(1-p),
        TRUE ~ 2*p*q*(1-theta)
      )
    ) %>%
    dplyr::ungroup()
    rmp <- prod(locus_res$gen_prob)
    list(locus_results = locus_res,
         rmp = rmp)
}

#' Calculate the population breakdown of samples
#'
#' @param file The input file as dataframe.
#' @param column (string) The basis for grouping samples, it should exist as a column in the dataframe.
#'
#' @returns The dataframe of tally count.
pop_breakdown <- function(file, column) {
  col_name <- as.character(column)

  if (!(col_name %in% colnames(file))) {
    stop("Ensure specified column is present in the dataset.")
  }

  file <- file %>%
    dplyr::rename(Sample = 1, Pop = col_name)

  for_breakdown <- data.frame(file$Sample, file$Pop)
  for_breakdown <- for_breakdown %>%
    dplyr::rename(Sample = 1, Population = 2)

  total <- for_breakdown %>%
    unique() %>%
    dplyr::group_by(Population) %>%
    dplyr::summarize(Total = n())

  return(as.data.frame(total))
}
