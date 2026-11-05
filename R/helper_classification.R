#' Perform Naive Bayes classification on forensic SNPs
#'
#' @param file The file path to file with sample, population, and genotype information.
#'
#' @returns A list containing the prediction rate and other model metrics.
calculate_loo <- function(file) {
  #data_fsnps <- load_csv_xlsx_files(file)
  #data_fsnps <- clean_input_data(data_fsnps)
  validate_file(file)
  data_fsnps <- dplyr::rename(data_fsnps, Sample = 1, Pop = 2)
  data_fsnps[] <- lapply(data_fsnps, factor)
  predictors <- !grepl("Pop", colnames(data_fsnps))
  label <- "Pop"
  data_fsnps <- as.data.frame(data_fsnps)

  # training the naive bayes classifier using a leave-one-out cross validation method
  res <- lapply(1:nrow(data_fsnps), function(i) {
    fit <- e1071::naiveBayes(
      y = factor(data_fsnps[-i, label]),
      x = as.matrix(data_fsnps[-i, predictors])
    )
    data.frame(
      label = data_fsnps[i, label],
      pred = predict(fit, as.matrix(data_fsnps[i, predictors], nrow = 1))
    )
  })

  # summarise predictions
  res <- do.call(rbind, res)

  # prepare the confusion matrix from the results
  confMatrix <- caret::confusionMatrix(res$pred, data_fsnps$Pop, mode = "everything")

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
  predictions <- data.frame(Sample, res)

  return(list(
    predTable = predWide,
    predStat = predStat,
    otherStat = otherStat,
    preds = predictions
  ))
}

calculate_crossval <- function(file, optimize = FALSE, k = 5) {
  validate_file(file)
  data_fsnps <- dplyr::rename(file, Sample = 1, Pop = 2)
  data_fsnps$Pop <- factor(data_fsnps$Pop)
  predictor_cols <- setdiff(colnames(data_fsnps), c("Sample", "Pop"))
  data_fsnps[predictor_cols] <- lapply(data_fsnps[predictor_cols], factor)
  
  pop_counts <- table(data_fsnps$Pop)
  
  if (any(pop_counts < k)) {
    stop("Each population must contain at least k samples.")
  }
  
  set.seed(123)
  
  tr_control <- caret::trainControl(
    method = "repeatedcv",
    number = k,
    repeats = 3,
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
  
  predictions_df <- data.frame(Sample = data_fsnps$Sample[cv_preds$rowIndex], 
                               label = ref_class, 
                               pred = pred_class)
  
  return(list(
    predTable = predWide,
    predStat = predStat,
    otherStat = otherStat,
    preds = predictions_df,
    model = tune_data,
    bestTune = tune_data$bestTune
  ))
}

