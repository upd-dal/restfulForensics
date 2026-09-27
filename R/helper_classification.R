#' Perform Naive Bayes classification on forensic SNPs
#'
#' @param file The file path to file with sample, population, and genotype information.
#'
#' @returns A list containing the prediction rate and other model metrics.
calculate_naive_bayes <- function(file) {
  data_fsnps <- load_csv_xlsx_files(file)
  data_fsnps <- clean_input_data(data_fsnps)
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
