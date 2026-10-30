plot_structure_server <- function(input, output, session, rv) {
  callModule(profvis_server, "profiler_strEvanno")
  output.dir <- tempdir()
  evannoRes <- reactiveVal(NULL)
  
  observe({
    inhouse_res <- !is.null(rv$structureRes) && (isTRUE(input$useSR))
    new_files <- (!isTRUE(input$useSR)) && !is.null(input$structureFilesEvanno)
    file_ready <- inhouse_res || new_files
    shinyjs::toggleState("runEvanno", condition = file_ready)
  })
  
  observeEvent(input$runEvanno, {
    disable("runEvanno")
    
    withProgress(message = "Analysis ongoing...", {
      incProgress(0.2, detail = "Loading input file...")
      if (isTRUE(input$useSR)) {
        sr <- rv$structureRes
      } else if (!is.null(input$structureFilesEvanno)) {
        str_results <- unpack_input_file(input$structureFilesEvanno$datapath, output.dir = output.dir)
        matrices_directory <- str_results$data_path
        
        # list files
        log_files <- list.files(matrices_directory, pattern = "_out_f", full.names = TRUE)
        
        if (length(log_files) == 0) {
          stop("No STRUCTURE _f files found.")
        }
        
        # get the pops 
        incProgress(0.4, detail = "Generating indfiles...")
        pops <- structure_get_populations(log_files[1])
        sr <- structureImport(log_files, pops = pops)
      }
      
      evanno_res <- strataG::evanno(sr, plot = TRUE)
      evannoRes(evanno_res)
    })
  })
  
  output$evannoTable <- DT::renderDataTable(
    {
      as.data.frame(evannoRes()$df)
    },
    options = list(pageLength = 10, scrollX = TRUE)
  )
  
  output$downloadEvannoTable <- downloadHandler(
    filename = function() {
      "evanno_results.csv"
    },
    content = function(file) {
      req(evannoRes())
      write.csv(as.data.frame(evannoRes()$df),
                file, row.names = FALSE, na = "")
    }
  )
  
  output$evannoMeanK <- renderPlot(
    {
      req(evannoRes())
      print(evannoRes()$plot$mean.ln.k)
    }
  )
  
  output$evannoLnPK <- renderPlot(
    {
      req(evannoRes())
      print(evannoRes()$plots$ln.pk)
    }
  )
  
  output$evannoPPK <- renderPlot(
    {
      req(evannoRes())
      print(evannoRes()$plots$ln.ppk)
    }
  )
  
  output$evannoDelta <- renderPlot(
    {
      req(evannoRes())
      print(evannoRes()$plots$delta.k)
    }
  )
  
  download_evanno_plots <- function(plot_name, file_name) {
    downloadHandler(
      filename = function() {
        file_name
      }, content = function(file) {
        req(evannoRes())
        p <- evannoRes()$plots[[plot_name]]
        req(!is.null(p))
        
        ggplot2::ggsave(
          filename = file,
          plot = p,
          device = "png",
          width = 8,
          height = 6,
          units = "in",
          dpi = 300
        )
      }
    )
  }
  
  output$downloadevannoMeanK <- download_evanno_plots(
    "mean.ln.k", "evanno_mean_ln_k.png"
  )
  
  output$downloadevannoLnPK <- download_evanno_plots(
    "ln.pk", "evanno_ln_pk.png"
  )
  
  output$downloadevannoPPK <- download_evanno_plots(
    "ln.ppk", "evanno_ln_ppk.png"
  )
  
  output$downloadevannoDelta <- download_evanno_plots(
    "delta.k", "evanno_delta_k.png"
  )
  
  output$downloadEvannoTable_UI <- renderUI({
    req(evannoRes())
    downloadButton("downloadEvannoTable", "Download Table (CSV)")
  })
  
  output$downloadevannoDelta_UI <- renderUI({
    req(evannoRes())
    downloadButton("downloadevannoDelta", "Download Delta K")
  })
  
  output$downloadevannoMeanK_UI <- renderUI({
    req(evannoRes())
    downloadButton("downloadevannoMeanK", "Download Mean Ln(K)")
  })
  
  output$downloadevannoLnPK_UI <- renderUI({
    req(evannoRes())
    downloadButton("downloadevannoLnPK", "Download Mean LnP(K)")
  })
  
  output$downloadevannoPPK_UI <- renderUI({
    req(evannoRes())
    downloadButton("downloadevannoPPK", "Download Mean LnP''(K)")
  })
  
  #==================== Plot STR Results ========================#
  callModule(profvis_server, "profiler_strClumpp")
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
    p <- structurePlot2(clumppResults(), plot = FALSE, type = "area")
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