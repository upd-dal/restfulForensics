classification_server <- function(input, output, session, rv) {
  # ================== CLASSIFICATION =====================#

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
      
      incProgress(0.2, detail = "Loading and calculating...")
      tryCatch(
      {
        result <- calculate_naive_bayes(input$forPredFile$datapath)
        showNotification("Calculation complete!", type = "message")
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
    
    incProgress(0.4, detail = "Processing results...")
    stats <- as.data.frame(result$predStat)
    stats <- data.frame(rowname(stats), stats)
    other_stat <- t(as.data.frame(result$otherStat))
    other_stat <- data.frame(rowname(other_stat), other_stat)
    
    predTable(as.data.frame(result$predTable))
    predStat(stats)
    predModel(other_stat)
    PredictionList(result$preds)
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
    downloadButton("downloadClassification", "Download Results")
  })
}
