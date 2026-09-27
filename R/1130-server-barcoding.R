barcoding_server <- function(input, output, session, rv) {
  # ===================== BARCODING =======================#
  #------------------ Species identification
  observe({
    refReady <- !is.null(input$refBarcoding$datapath)
    queReady <- !is.null(input$queBarcoding$datapath)
    toggleState("identifySpecies", refReady && queReady)
  })

  resultIdentity <- reactiveVal(NULL)
  refseq <- reactiveVal(NULL)
  queseq <- reactiveVal(NULL)

  observeEvent(input$identifySpecies, {
    disable("identifySpecies")

    req(input$refBarcoding)
    req(input$queBarcoding)

    tryCatch({
      ref_seq <- alignment_to_dnabin(input$refBarcoding$datapath)
    que_seq <- alignment_to_dnabin(input$queBarcoding$datapath)

    refseq(ref_seq)
    queseq(que_seq)

    # If not using kmer method
    if (!isTRUE(input$kmerSelect)) {
      result_identity <- BarcodingR::barcoding.spe.identify(
        refseq(),
        queseq(),
        method = input$barcodingMethod
      )

      if (input$barcodingMethod == "Bayesian") {
        result <- result_identity[["output_identified"]]
      } else {
        result <- list(
          success_rate = result_identity[["success.rates.ref"]],
          output = result_identity[["output_identified"]]
        )
      }
    } else {
      if (input$kmerType == "Fuzzy-set Method and kmer") {
        result_identity <- BarcodingR::barcoding.spe.identify2(
          refseq(),
          queseq(),
          kmer = input$kmerValueFuzzy,
          optimization = input$optimizationKMER
        )

        result <- list(
          model_success = result_identity[["model.success"]],
          output = result_identity[["output_identified"]]
        )
      } else if (input$kmerType == "BP-based Method and kmer") {
        result_identity <- BarcodingR::bbsik(
          refseq(),
          queseq(),
          kmer = input$kmerValueBP,
          UseBuiltModel = input$builtModel,
          lr = input$lrValue,
          maxit = input$maxitValue
        )

        result <- list(
          success_rate = result_identity[["success.rates.ref"]],
          output = result_identity[["output_identified"]]
        )
      }
    }

    resultIdentity(result)
    showNotification("Calculation complete!", type = "message")
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    }
      )
    
    enable("identifySpecies")
  })

  output$identificationResult <- renderPrint({
    req(resultIdentity())
    if (!isTRUE(input$kmerSelect)) {
      if (input$barcodingMethod == "Bayesian") {
        print("Species identification success rate unavailable when using 'Bayesian' method.")
      }
      rate <- resultIdentity()[[1]]
      print(rate * 100)
    } else {
      rate <- resultIdentity()[[1]]
      print(rate * 100)
    }
  })

  output$identified_samples <- DT::renderDataTable({
    req(resultIdentity())
    if (!isTRUE(input$kmerSelect)) {
      if (input$barcodingMethod == "Bayesian") {
        DT::datatable(resultIdentity()[[1]],
          options = list(pageLength = 5),
          rownames = FALSE
        )
      }
    }

    DT::datatable(resultIdentity()[[2]],
      options = list(pageLength = 5),
      rownames = FALSE
    )
  })

  #----------------- Optimize kmer values
  observe({
    fileReady <- !is.null(input$optimizeKmerRef)
    toggleState("calOptimumKmer", fileReady)
  })

  kmerFile <- reactiveVal(NULL)
  optimalKmer <- reactiveVal(NULL)

  observeEvent(input$calOptimumKmer, {
    disable("calOptimumKmer")
    req(input$optimizeKmerRef)

    Sys.sleep(1.5)

    kmer_File <- alignment_to_dnabin(input$optimizeKmerRef$datapath)

    tmp_file <- tempfile(fileext = ".png")
    png(tmp_file, width = 1200, height = 800)

    optimal_Kmer <- BarcodingR::optimize.kmer(
      kmer_File,
      max.kmer = as.numeric(input$maxKmer)
    )

    dev.off()
    kmerFile(tmp_file)
    optimalKmer(optimal_Kmer)

    enable("calOptimumKmer")
  })

  output$kmerResult <- renderPrint({
    req(optimalKmer())
    as.data.frame(optimalKmer())
  })

  output$kmerPlot <- renderImage(
    {
      req(kmerFile())

      list(
        src = kmerFile(),
        contentType = "image/png",
        alt = "Optimum Kmer Plot"
      )
    },
    deleteFile = FALSE
  )

  output$downloadKmerPlot <- downloadHandler(
    filename = function() {
      paste0("optimum_kmer_", Sys.Date(), ".png")
    },
    content = function(file) {
      req(kmerFile())
      req(input$maxKmer)
      BarcodingR::optimize.kmer(kmerFile(), max.kmer = input$maxKmer)
    }, contentType = "image/png"
  )

  output$downloadKmerPlot_UI <- renderUI({
    req(kmerFile())
    downloadButton("downloadKmerPlot", "Download Kmer Plot")
  })

  #---------------------- barcoding gap
  observe({
    gapReady <- !is.null(input$barcodeRef)
    toggleState("gapBarcodes", gapReady)
  })

  refBarcode <- reactiveVal(NULL)
  barcodeGap <- reactiveVal(NULL)
  gapPlotFile <- reactiveVal(NULL)

  observeEvent(input$gapBarcodes, {
    disable("gapBarcodes")
    req(input$barcodeRef)
    req(input$gapModel)

    ref_Barcode <- alignment_to_dnabin(input$barcodeRef$datapath)
    refBarcode(ref_Barcode)

    tmp_plot <- tempfile(fileext = ".png")
    png(tmp_plot, width = 1200, height = 800)
    gap <- tryCatch(
      BarcodingR::barcoding.gap(ref_Barcode, dist = input$gapModel),
      error = function(e) {
        message(e$message)
        NULL
      }
    )

    dev.off()

    barcodeGap(gap)
    gapPlotFile(tmp_plot)
    enable("gapBarcodes")
  })

  output$barcodingResult <- renderPrint({
    req(barcodeGap())
    as.data.frame(barcodeGap())
  })

  output$BarcodingGapPlot <- renderImage(
    {
      req(gapPlotFile())

      list(
        src = gapPlotFile(),
        contentType = "image/png",
        alt = "Barcoding Gap Plot"
      )
    },
    deleteFile = FALSE
  )

  output$downloadGapPlot <- downloadHandler(
    filename = function() {
      paste0("barcoding_gap_", Sys.Date(), ".png")
    },
    content = function(file) {
      req(gapPlotFile())
      file.copy(gapPlotFile(), file, overwrite = TRUE)
    }, contentType = "image/png"
  )

  output$downloadGapPlot_UI <- renderUI({
    req(gapPlotFile())
    downloadButton("downloadGapPlot", "Download Gap Plot")
  })

  #--------------------------- barcodes eval
  observe({
    barcode1ready <- !is.null(input$barcode1)
    barcode2ready <- !is.null(input$barcode2)
    toggleState("evalBarcodes", barcode1ready && barcode2ready)
  })

  barcode1RV <- reactiveVal(NULL)
  barcode2RV <- reactiveVal(NULL)
  barcodeEvalRV <- reactiveVal(NULL)

  observeEvent(input$evalBarcodes, {
    disable("evalBarcodes")
    req(input$barcode1)
    req(input$barcode2)
    req(input$kmer1)
    req(input$kmer2)

    b1 <- alignment_to_dnabin(input$barcode1$datapath)
    b2 <- alignment_to_dnabin(input$barcode2$datapath)
    barcode1RV(b1)
    barcode2RV(b2)

    # convert to dataframe to download
    result <- tryCatch(
      {
        BarcodingR::barcodes.eval(b1,
          b2,
          kmer1 = as.numeric(input$kmer1),
          kmer2 = as.numeric(input$kmer2)
        )
      },
      error = function(e) {
        message(e$message)
        NULL
      }
    )
    barcodeEvalRV(result)
  })

  output$evalBarcodesResult <- renderTable(
    {
      req(barcodeEvalRV())
      as.data.frame(barcodeEvalRV())
    },
    rownames = TRUE
  )

  #---------------------------- tdr2
  observe({
    file1Ready <- !is.null(input$oneSpe)
    file2Ready <- !is.null(input$queSpe)
    toggleState("calculateTDR2", file1Ready && file2Ready)
  })

  oneSpeRV <- reactiveVal(NULL)
  queSpeRV <- reactiveVal(NULL)
  tdrRV <- reactiveVal(NULL)

  observeEvent(input$calculateTDR2, {
    disable("calculateTDR2")
    req(input$oneSpe)
    req(input$queSpe)
    req(input$bootValue1)
    req(input$bootValue2)

    oneSpe <- alignment_to_dnabin(input$oneSpe$datapath)
    queSpe <- alignment_to_dnabin(input$queSpe$datapath)

    oneSpeRV(oneSpe)
    queSpeRV(queSpe)

    result <- tryCatch(
      {
        BarcodingR::TDR2(oneSpe, queSpe, boot = as.numeric(input$bootValue1), boot2 = as.numeric(input$bootValue2))
      },
      error = function(e) {
        message(e$message)
        NULL
      }
    )

    tdrRV(result)
    enable("calculateTDR2")
  })

  # issue with results, it prints and not stores
  output$tdrValues <- renderTable(
    {
      req(tdrRV())
      data.frame(TDR = tdrRV())
    },
    rownames = TRUE
  )
}
