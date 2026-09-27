run_structure_analysis <- function(input, output, session, rv) {
  
  structure_exec <- normalizePath(
    "./structure/structure.exe",
    mustWork = TRUE
  )
  
  examplePop_STR2 <- data.frame(
    Sample = c("Sample1", "Sample2", "Sample3", "Sample4", "..."),
    Population = c("POP1", "POP2", "POP3", "POP4", "..."),
    rs101 = c("A/A", "A/T", "A/A", "T/T", "..."),
    rs102 = c("G/G", "C/C", "G/C", "G/G", "..."),
    rs_n = c("...", "...", "...", "...", "...")
  )
  
  output$examplePop_STR2UI <- DT::renderDataTable(
    {
      req(examplePop_STR2)
      examplePop_STR2
    },
    options = list(
      scrollX = TRUE,
      pageLength = 5
    )
  )
  
  observe({
    file_ready <- !is.null(input$structureFile)
    shinyjs::toggleState("runStructure", condition = file_ready)
  })
  
  analysis_done <- reactiveVal(FALSE)
  structureLog <- reactiveVal("")
  qmatrices_result <- reactiveVal(NULL)
  output_dir <- tempfile("structure_")
  dir.create(output_dir)
  out_path <- output_dir
  
  observeEvent(input$runStructure, {
    old_wd <- getwd()
    setwd(output_dir)
    on.exit(setwd(old_wd), add = TRUE)
    disable("runStructure")
    req(input$structureFile)
    Sys.sleep(1.5)
    
    withProgress(message = "Analysis ongoing...", {
      incProgress(0.2, detail = "Loading input file...")
      df <- csv_to_gtype_format(input$structureFile$datapath)
      df_gtype <- strataG::df2gtypes(df, ploidy = 2, id.col = 1, strata.col = 2, loc.col = 3)
      
      incProgress(0.4, detail = "Running STRUCTURE analysis...")
      
      if (isTRUE(input$inferAlpha)){
        sr <- strataG::structureRun(df_gtype,
                                    k.range = input$kMin:input$kMax,
                                    num.k.rep = input$numKRep,
                                    burnin = input$burnin,
                                    numreps = input$numreps,
                                    noadmix = input$noadmix,
                                    freqscorr = input$freqScore,
                                    inferalpha = input$inferAlpha,
                                    alpha = input$alphaValStructure,
                                    exec = structure_exec,
                                    delete.files = FALSE,
                                    label = "structureRun",
                                    pop.prior = "usepopinfo"
        )
      } else {
        sr <- strataG::structureRun(df_gtype,
                                    k.range = input$kMin:input$kMax,
                                    num.k.rep = input$numKRep,
                                    burnin = input$burnin,
                                    numreps = input$numreps,
                                    noadmix = input$noadmix,
                                    freqscorr = input$freqScore,
                                    exec = structure_exec,
                                    delete.files = FALSE,
                                    label = "structureRun",
                                    pop.prior = "usepopinfo"
        )
      }
      stray_dir <- file.path(getwd(), "structureRun.structureRun")
      
      if (dir.exists(stray_dir)) {
        
        dest <- file.path(out_path, "structure_files")
        
        dir.create(dest, recursive = TRUE, showWarnings = FALSE)
        
        file.copy(
          list.files(stray_dir, full.names = TRUE),
          dest,
          recursive = TRUE,
          overwrite = TRUE
        )
        
        unlink(stray_dir, recursive = TRUE, force = TRUE)
      }
      
      rv$structureRes <- sr
      
      # get the q matrices
      incProgress(0.8, detail = "Extracting q matrices...")
      qmatrices_result(lapply(sr, function(x) x$q.mat))
      
      enable("runStructure")
    })
    analysis_done(TRUE)
  })
  
  #   output$structure_log <- renderText({
  #      structureLog()
  #   })
  #   
  #   output$downloadLogs <- downloadHandler(
  #      filename = function() {
  #         paste0("structure_logs_", Sys.Date(), ".zip")
  #      },
  #      content = function(file) {
  #         log_files <- list.files(out_path, pattern = "_log.*$", full.names = TRUE)
  #         if (length(log_files) == 0) {
  #            return(NULL)
  #         }
  #         
  #         log_files_renamed <- sapply(log_files, function(f) {
  #            ext <- tools::file_ext(f)
  #            if (ext == "") {
  #               new_f <- paste0(f, ".txt")
  #               file.rename(f, new_f)
  #               return(new_f)
  #            }
  #            return(f)
  #         })
  #         zip::zipr(zipfile = file, files = log_files_renamed, recurse = FALSE)
  #      },
  #      contentType = "application/zip"
  #   )
  
  output$downloadFOutputs <- downloadHandler(
    filename = function() {
      paste0("structure_outputs_", Sys.Date(), ".zip")
    },
    content = function(file) {
      
      f_files <- list.files(
        output_dir,
        pattern = "_out_f$",
        recursive = TRUE,
        full.names = TRUE
      )
      
      if (length(f_files) == 0) {
        stop("No STRUCTURE output files found.")
      }
      
      zip::zipr(zipfile = file, files = f_files)
    },
    contentType = "application/zip"
  )
  
  output$downloadQMatrixTxtZip <- downloadHandler(
    filename = function() {
      paste0("q_matrices_", Sys.Date(), ".zip")
    },
    content = function(file) {
      req(qmatrices_result())
      q_files <- tempfile()
      dir.create(q_files)
      lapply(names(qmatrices_result()), function(name) {
        mat <- qmatrices_result()[[name]]
        if (!is.null(mat)) {
          write.table(mat, file.path(q_files, paste0(name, ".indfile")),
                      row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t"
          )
        }
      })
      zip::zipr(zipfile = file, files = list.files(q_files, full.names = TRUE))
    },
    contentType = "application/zip"
  )
  
  output$downloadButtons <- renderUI({
    req(analysis_done())
    tagList(
      #downloadButton("downloadLogs", "Download Log Files (.zip)"),
      downloadButton("downloadFOutputs", "Download STRUCTURE _f Files (.zip)"),
      downloadButton("downloadQMatrixTxtZip", "Download Q Matrices (.zip)")
    )
  })
  
}