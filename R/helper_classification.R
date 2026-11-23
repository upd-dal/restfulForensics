#' Perform Naive Bayes classification on forensic SNPs
#'
#' @param file The file path to file with sample, population, and genotype information.
#'
#' @returns A list containing the prediction rate and other model metrics.
calculate_loo <- function(file) {
  data_fsnps <- load_csv_xlsx_files(file)
  data_fsnps <- clean_input_data(data_fsnps)
  # validate_file(file)
  data_fsnps <- dplyr::rename(data_fsnps, Sample = 1, Pop = 2)
  data_fsnps$Sample <- as.character(data_fsnps$Sample)
  data_fsnps$Pop <- factor(data_fsnps$Pop)

  if (any(is.na(data_fsnps$Pop))) {
    stop("Missing population labels detected.")
  }

  samples <- data_fsnps$Sample
  populations <- data_fsnps$Pop
  population_levels <- levels(populations)

  snp_cols <- setdiff(colnames(data_fsnps), c("Sample", "Pop"))

  if (length(snp_cols) == 0) {
    stop("No SNP predictors detected.")
  }

  observed_snps <- rowSums(!is.na(data_fsnps[, snp_cols, drop = FALSE]))

  # training the naive bayes classifier using a leave-one-out cross validation method
  res <- lapply(seq_len(nrow(data_fsnps)), function(i) {
    train <- data_fsnps[-i, , drop = FALSE]
    test <- data_fsnps[i, , drop = FALSE]

    usable_snps <- snp_cols[colSums(!is.na(
      train[, snp_cols, drop = FALSE]
    )) > 0]

    # if there are no available snps
    if (length(usable_snps) == 0 || observed_snps[i] == 0) {
      return(data.frame(
        label = as.character(populations[i]),
        pred = NA_character_,
        n_SNPs = 0
      ))
    }

    for (snp in usable_snps) {
      train[[snp]] <- droplevels(factor(train[[snp]]))
      test[[snp]] <- factor(
        as.character(test[[snp]]),
        levels = levels(train[[snp]])
      )
    }

    available_snps <- usable_snps[!is.na(test[1, usable_snps, drop = TRUE])]

    if (length(available_snps) == 0) {
      return(data.frame(
        label = as.character(populations[i]),
        pred = NA_character_,
        n_SNPs = 0
      ))
    }

    # train classifier
    model_formula <- stats::reformulate(available_snps, response = "Pop")
    fit <- e1071::naiveBayes(formula = model_formula, data = train, laplace = 1, na.action = na.pass)
    prediction <- predict(fit, newdata = test[, available_snps, drop = FALSE])
    data.frame(
      label = as.character(populations[i]),
      pred = as.character(prediction),
      n_SNPs = length(available_snps)
    )
  })


  # summarise predictions
  res <- do.call(rbind, res)
  res$label <- factor(res$label, levels = population_levels)
  res$pred <- factor(res$pred, levels = population_levels)
  valid <- !is.na(res$pred)

  if (!any(valid)) {
    stop("No samples were classified using the available data")
  }

  # prepare the confusion matrix from the results
  confMatrix <- caret::confusionMatrix(res$pred[valid], data_fsnps$Pop[valid], mode = "everything")

  # convert to table, get pred and ref values
  pred <- as.data.frame(confMatrix$table)

  predWide <- tidyr::pivot_wider(
    data = pred,
    names_from = Reference,
    values_from = Freq
  )

  # access accuracy values
  predStat <- as.data.frame(confMatrix$overall)

  # get other stats
  otherStat <- as.data.frame(confMatrix$byClass)

  Sample <- data_fsnps$Sample
  predictions <- data.frame(Sample = samples, res, stringsAsFactors = FALSE)
  missing_summary <- data.frame(
    Sample = samples,
    Population = populations,
    Observed_SNPs = observed_snps,
    Missing_SNPs = length(snp_cols) - observed_snps,
    Classified = valid
  )

  return(list(
    predTable = predWide,
    predStat = predStat,
    otherStat = otherStat,
    preds = predictions,
    missingData = missing_summary
  ))
}

calculate_crossval <- function(file, optimize = FALSE, k = 5, reps = 3) {
  data_fsnps <- dplyr::rename(file, Sample = 1, Pop = 2)
  data_fsnps$Pop <- factor(data_fsnps$Pop)
  predictor_cols <- setdiff(colnames(data_fsnps), c("Sample", "Pop"))
  data_fsnps[predictor_cols] <- lapply(data_fsnps[predictor_cols], factor)

  pop_counts <- table(data_fsnps$Pop)
  lowest_pop <- names(pop_counts)[which.min(pop_counts)]
  lowest_count <- min(pop_counts)

  if (any(pop_counts < k)) {
    stop(paste(
      "Each population must contain at least k samples.",
      lowest_pop, "only has", lowest_count, "values."
    ))
  }

  set.seed(123)

  tr_control <- caret::trainControl(
    method = "repeatedcv",
    number = k,
    repeats = reps,
    savePredictions = "final"
  )

  if (isTRUE(optimize)) {
    tune_grid <- expand.grid(
      laplace = c(0, 0.5, 1, 2),
      usekernel = c(TRUE, FALSE),
      adjust = c(0.5, 1.0, 1.5)
    )
  } else {
    tune_grid <- expand.grid(
      laplace = 0,
      usekernel = FALSE,
      adjust = 1
    )
  }

  tune_data <- caret::train(
    Pop ~ . - Sample,
    data = data_fsnps,
    method = "naive_bayes",
    trControl = tr_control,
    tuneGrid = tune_grid
  )

  cv_preds <- tune_data$pred

  pred_class <- factor(as.character(cv_preds$pred), levels = levels(data_fsnps$Pop))
  ref_class <- factor(as.character(cv_preds$obs), levels = levels(data_fsnps$Pop))

  # prepare the confusion matrix from the results
  confMatrix <- caret::confusionMatrix(data = pred_class, reference = ref_class, mode = "everything")

  # convert to table, get pred and ref values
  pred <- as.data.frame(confMatrix$table)

  predWide <- tidyr::pivot_wider(
    data = pred,
    names_from = Reference,
    values_from = Freq
  )

  # access accuracy values
  predStat <- as.data.frame(confMatrix$overall)

  # get other stats
  otherStat <- as.data.frame(confMatrix$byClass)

  predictions_df <- data.frame(
    Sample = data_fsnps$Sample[cv_preds$rowIndex],
    label = ref_class,
    pred = pred_class
  )

  return(list(
    predTable = predWide,
    predStat = predStat,
    otherStat = otherStat,
    preds = predictions_df,
    model = tune_data,
    bestTune = tune_data$bestTune
  ))
}
