plot_structure_server <- function(input, output, session, rv) {
  
  populationNames <- reactiveVal(NULL)
  structureResults <- reactiveVal(NULL)
  clumppResults <- reactiveVal(NULL)
  strPlot <- reactiveVal(NULL)
  
  observe({
    inhouse_res <- !is.null(rv$structureRes) && (isTRUE(input$plotSR))
    new_files <- (!isTRUE(input$plotSR)) && !is.null(input$structureFilesZipped)
    file_ready <- inhouse_res || new_files
    shinyjs::toggleState("plotStructureResults", condition = file_ready)
  })
  
  output.dir <- tempdir()
  
  observeEvent(input$plotStructureResults, {
    disable("plotStructureResults")
    
    withProgress(message = "Analysis ongoing...", {
      incProgress(0.2, detail = "Loading input file...")
      
      if (isTRUE(input$plotSR)) {
        sr <- rv$structureRes
      } else if (!is.null(input$structureFilesZipped)) {
        str_results <- unpack_input_file(input$structureFilesZipped$datapath, output.dir = output.dir)
        matrices_directory <- str_results$data_path
        
        # list files
        log_files <- list.files(matrices_directory, pattern = "_out_f", full.names = TRUE)
        
        if (length(log_files) == 0) {
          stop("No STRUCTURE _f files found.")
        }
        
        # get the pops 
        incProgress(0.4, detail = "Generating indfiles...")
        pops <- structure_get_populations(log_files[1])
        populationNames(pops)
        sr <- structureImport(log_files, pops = pops)
      }
      
      structureResults(sr)
      incProgress(0.6, detail = "Running clumpp...")
      # run clumpp
      out_clumpp <- file.path(output.dir, "clumpp_res")
      clumpp_res <- strataG::clumpp(
        sr,
        k = input$kChoice,
        sim.stat = input$simStat,
        align.algorithm = input$greedyAlgo,
        greedy.option = input$greedyOption,
        repeats = input$repeatsClumpp,
        order.by.run = input$orderRun,
        label = out_clumpp,
        delete.files = FALSE
      )
      
      clumppResults(clumpp_res)
      
    })
    
    enable("plotStructureResults")
  })
  
  output$structure_result_plots <- renderPlot({
    req(clumppResults())
    p <- strataG::structurePlot(clumppResults(), plot = FALSE)
    strPlot(p)
    p
  }, res = 120)
  
  output$downloadStructurePlot <- downloadHandler(
    filename = function() {
      paste0("STRUCTURE_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      req(strPlot())
      
      ggplot2::ggsave(
        filename = file,
        plot = strPlot(),
        width = 10,
        height = 6,
        dpi = 300
      )
    },
    contentType = "image/png"
  )
  
  output$downloadStructurePDF <- downloadHandler(
    filename = function() {
      paste0("STRUCTURE_plot_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      req(strPlot())
      
      ggplot2::ggsave(
        filename = file,
        plot = strPlot(),
        width = 10,
        height = 6,
        device = "pdf"
      )
    },
    contentType = "application/pdf"
  )
  
  output$downloadStructure_UI <- renderUI({
    req(strPlot())
    tagList(
      downloadButton("downloadStructurePlot", "Download STRUCTURE Plot (PNG)"),
      downloadButton("downloadStructurePDF", "Download STRUCTURE Plot (PDF)")
    )
  })
  
}