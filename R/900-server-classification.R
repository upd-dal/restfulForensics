classification_server <- function(input, output, session, rv) {
  # ================== LOO =====================#
  classificationRef <- data.frame(
    Sample = c("Sample1", "Sample2", "Sample3", "Sample4", "..."),
    Population = c("POP1", "POP2", "POP3", "POP4", "..."),
    rs101 = c("A/A", "A/T", "A/A", "T/T", "..."),
    rs102 = c("G/G", "C/C", "G/C", "G/G", "..."),
    rs_n = c("...", "...", "...", "...", "...")
  )

  output$classificationRef_UI <- DT::renderDataTable(
    {
      req(classificationRef)
      classificationRef
    },
    options = list(
      scrollX = TRUE,
      pageLength = 5
    )
  )

  observe({
    file_ready <- !is.null(input$forPredFile)
    shinyjs::toggleState("runNaiveBayes", condition = file_ready)
  })

  predTable <- reactiveVal(NULL)
  predStat <- reactiveVal(NULL)
  predModel <- reactiveVal(NULL)
  PredictionList <- reactiveVal(NULL)
  
  observeEvent(input$runNaiveBayes, {
    disable("runNaiveBayes")
    req(input$forPredFile)

    withProgress(message = "Ongoing: ", {
      tryCatch(
        {
      incProgress(0.2, detail = "Loading and calculating...")
          file_ext <- tools::file_ext(input$forPredFile$datapath)
          if (!file_ext %in% c(".csv", "xlsx", "xlsm", "xlsb", "xls")) {
            stop("Unsupported file format: Should be a CSV/Excel file")
          }
          
        validate_file(load_csv_xlsx_files(input$forPredFile$datapath))
        calculate_loo(input$forPredFile$datapath)
    
    incProgress(0.4, detail = "Processing results...")
      stats <- as.data.frame(result$predStat)
      stats <- data.frame(row.names(stats), stats)
      other_stat <- t(as.data.frame(result$otherStat))
      other_stat <- data.frame(row.names(other_stat), other_stat)
      
      predTable(as.data.frame(result$predTable))
      predStat(stats)
      predModel(other_stat)
      PredictionList(result$preds)

      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
        return(NULL)
      })
      
    })
    enable("runNaiveBayes")
  })

  output$predictionTableResult <- DT::renderDataTable({
     req(predTable())
     DT::datatable(
        predTable(),
        options = list(
           pageLength = 10,
           autoWidth = TRUE,
           searchHighlight = TRUE
        ),
        filter = "top",
        selection = "multiple"
     )
  })

  output$statbyClassResult <- DT::renderDataTable({
    req(predModel())
     DT::datatable(
        predModel(),
        options = list(
           pageLength = 10,
           autoWidth = TRUE,
           searchHighlight = TRUE
        ),
        filter = "top",
        selection = "multiple"
     )
  })

  output$overallStatResult <- DT::renderDataTable({
    req(predStat())
     DT::datatable(
        predStat(),
        options = list(
           pageLength = 10,
           autoWidth = TRUE,
           searchHighlight = TRUE
        ),
        filter = "top",
        selection = "multiple"
     )
  })
  
  output$predictionList <- DT::renderDataTable({
     req(PredictionList())
     DT::datatable(
        PredictionList(),
        options = list(
           pageLength = 10,
           autoWidth = TRUE,
           searchHighlight = TRUE
        ),
        filter = "top",
        selection = "multiple"
     )
  })

  output$downloadClassification <- downloadHandler(
    filename = function() {
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
      paste0("classification-results_", timestamp, ".xlsx")
    },
    content = function(file) {
      dataset <- list(
        "Table" = predTable(),
        "Stats per Class" = predModel(),
        "Overall Stats" = predStat(),
        "Predictions by Individual" = PredictionList()
      )
      openxlsx::write.xlsx(dataset, file = file)
    }
  )

  output$downloadClassification_UI <- renderUI({
    req(predTable(), predModel(), predStat(), PredictionList())
    downloadButton("downloadClassification", "Download Results")
  })
  
  #=================CROSS VALIDATION
  classificationRefCV <- data.frame(
    Sample = c("Sample1", "Sample2", "Sample3", "Sample4", "..."),
    Population = c("POP1", "POP2", "POP3", "POP4", "..."),
    rs101 = c("A/A", "A/T", "A/A", "T/T", "..."),
    rs102 = c("G/G", "C/C", "G/C", "G/G", "..."),
    rs_n = c("...", "...", "...", "...", "...")
  )
  
  output$classificationRefCV_UI <- DT::renderDataTable(
    {
      req(classificationRefCV)
      classificationRefCV
    },
    options = list(
      scrollX = TRUE,
      pageLength = 5
    )
  )
  
  predTableCV <- reactiveVal(NULL)
  predStatCV <- reactiveVal(NULL)
  predModelCV <- reactiveVal(NULL)
  PredictionListCV <- reactiveVal(NULL)
  modelInfoCV <- reactiveVal(NULL)
  
  observe({
    file_ready <- !is.null(input$crossValFile)
    shinyjs::toggleState("performCrossval", file_ready)
  })
  
  observeEvent(input$performCrossval, {
    disable("performCrossval")
    req(input$crossValFile)
    withProgress("Performing cross-validation", {
      
  tryCatch({
    incProgress(0.2, detail = "Loading and validating input data...")
    file_ext <- tools::file_ext(input$crossValFile$datapath)
    if (!file_ext %in% c(".csv", "xlsx", "xlsm", "xlsb", "xls")) {
      stop("Unsupported file format: Should be a CSV/Excel file")
    }
        input.file <- load_csv_xlsx_files(input$crossValFile$datapath)
        input.file <- clean_input_data(input.file)
        validate_file(input.file)
        calculate_crossval(file = input.file,
                           optimize = input$optimizeVals,
                           k = as.integer(input$kValue))
      
        incProgress(0.5, detail = "Processing results...")
        stats <- as.data.frame(result$predStat)
        stats <- data.frame(
          Statistic = rownames(stats), stats, row.names = NULL
        )
        
        other_stat <- as.data.frame(result$otherStat)
        other_stat <- data.frame(Population = rownames(other_stat), other_stat, row.names = NULL)
        
        predTableCV(result$predTable)
        predStatCV(stats)
        predModelCV(other_stat)
        PredictionListCV(result$preds)
        
        model_info <- data.frame(
          Setting = c("Validation method", "Number of folds", "Repeats", "Optimization", "Best tuning parameters"),
          Value = c("Repeated K-fold cross-validation", input$kValue, 3, if(isTRUE(input$optimizeVals)) "Yes" else "No",
                    paste(names(result$bestTune), result$bestTune[1, ], sep = " = ", collapse = "; "))
        )
        
        modelInfoCV(model_info)
    
      }, error = function(e) {
        showNotification(paste("Cross-validation error:", e$message), type = "error")
        return(NULL)
      })
      
    }
      
    )
    enable("performCrossval")
  })
  
  output$predictionTableResultCV <- DT::renderDataTable({
    req(predTableCV())
    
    DT::datatable(predTableCV(),
                  options = list(
                    pageLength = 10, 
                    autoWidth = TRUE,
                    searchHighlight = TRUE
                  ), filter = "top", selection = "multiple")
  })
  
  output$statbyClassResultCV <- DT::renderDataTable({
    req(predModelCV())
    
    DT::datatable(predModelCV(),
                  options = list(
                    pageLength = 10, 
                    autoWidth = TRUE,
                    searchHighlight = TRUE
                  ), filter = "top", selection = "multiple")
  })
  
  output$statbyStatResultCV <- DT::renderDataTable({
    req(predStatCV())
    
    DT::datatable(predStatCV(),
                  options = list(
                    pageLength = 10, 
                    autoWidth = TRUE,
                    searchHighlight = TRUE
                  ), filter = "top", selection = "multiple")
  })
  
  output$predictionListCV <- DT::renderDataTable({
    req(PredictionListCV())
    
    DT::datatable(PredictionListCV(),
                  options = list(
                    pageLength = 10, 
                    autoWidth = TRUE,
                    searchHighlight = TRUE
                  ), filter = "top", selection = "multiple")
  })
  
  output$modelInfoCV <- DT::renderDataTable({
    req(modelInfoCV())
    
    DT::datatable(modelInfoCV(),
                  options = list(
                    pageLength = 10, 
                    autoWidth = TRUE,
                    searching = FALSE,
                    ordering = FALSE
                  ))
  })
  
  output$downloadPredictionTableCV <- downloadHandler(
    filename = function() {
      paste0("CV_confusion_matrix_", format(Sys.time(), "%Y%m%d"), ".xlsx")
    }, content = function(file) {
      openxlsx::write.xlsx(predTableCV(), file = file)
    }
  )
  
  output$downloadPredictionTableCV_UI <- renderUI({
    req(predTableCV())
    
    downloadButton("downloadPredictionTableCV", "Download Confusion matrix")
  })
  
  output$downloadstatbyClassResultCV <- downloadHandler(
    filename = function() {
      paste0("CV_statistics_bypop_", format(Sys.time(), "%Y%m%d"), ".xlsx")
    }, content = function(file) {
      openxlsx::write.xlsx(predModelCV(), file = file)
    }
  )
  
  output$downloadstatbyClassResultCV_UI <- renderUI({
    req(predModelCV())
    
    downloadButton("downloadstatbyClassResultCV", "Download Population Statistics")
  })
  
  output$downloadoverallStatResultCV <- downloadHandler(
    filename = function() {
      paste0("CV_overall_statistics_", format(Sys.time(), "%Y%m%d"), ".xlsx")
    }, content = function(file) {
      openxlsx::write.xlsx(predStatCV(), file = file)
    }
  )
  
  output$downloadstatbyClassResultCV_UI <- renderUI({
    req(predStatCV())
    
    downloadButton("downloadoverallStatResultCV", "Download Overall Statistics")
  })
  
  output$downloadpredictionListCV <- downloadHandler(
    filename = function() {
      paste0("CV_individual_predictions_", format(Sys.time(), "%Y%m%d"), ".xlsx")
    }, content = function(file) {
      openxlsx::write.xlsx(PredictionListCV(), file = file)
    }
  )
  
  output$downloadpredictionListCV_UI <- renderUI({
    req(PredictionListCV())
    
    downloadButton("downloadpredictionListCV", "Download Individual Predictions")
  })
  
  
}
