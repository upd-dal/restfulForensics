filtering_server <- function(input, output, session, rv) {
  # ===================== FILTERING =======================#

  observe({
    hasFile <- !is.null(input$markerFileFilter)
    anyFilter <- input$filterIndiv || input$filterVariant || input$filterAllele || input$filterQuality || input$filterHWE || input$filterLD || input$cutoffKing
    shinyjs::toggleState("calcDP", condition = hasFile && anyFilter)
  })

  output$filterWarning <- renderText({
    if (!(input$filterIndiv || input$filterVariant || input$filterAllele || input$filterQuality || input$filterHWE || input$filterLD || input$cutoffKing)) {
      "Please select at least one filtering option."
    } else {
      ""
    }
  })

  plink2_path <- get_plink2_path()

  # reset buttons with new file
  observeEvent(input$markerFileFilter, {
    ext <- tools::file_ext(input$markerFileFilter$name)
    if (tolower(ext) == "vcf") {
      shinyjs::enable("enableDP")
    } else {
      updateCheckboxInput(inputId = "enableDP", value = FALSE)
      shinyjs::disable("enableDP")
    }
  })

  temp_dir <- tempdir()
  depth_outputs <- reactiveVal(NULL)
  filtered_plink_file <- reactiveVal(NULL)

  observeEvent(input$calcDP, {
    req(input$markerFileFilter)

    vcf_path <- input$markerFileFilter$datapath
    ref_path <- if (!is.null(input$highlightRef)) input$highlightRef$datapath else NULL
    palette <- if (!is.null(input$colorPalette)) input$colorPalette else NULL
    ext <- tools::file_ext(input$markerFileFilter$name)

  withProgress(message = "Ongoing...", {
    incProgress(0.2, detail = "Loading input file...")
  if (tolower(ext) == "vcf" && input$enableDP) {
      dp <- depth_from_vcf(
        vcf = vcf_path,
        output.dir = temp_dir,
        reference = ref_path,
        palette = palette
      )
      depth_outputs(dp)
    }

    # CHECK THE FILE EXTENSIONS
    input_type <- if (!is.null(input$markerFileFilter)) {
      if (grepl("\\.bcf$", input$markerFileFilter$name, ignore.case = TRUE)) {
        "bcf"
      } else {
        "vcf"
      }
    } else {
      "plink"
    }
    incProgress(0.4, detail = "Preparing file for filtering...")
    # convert
    pgen_prefix <- file.path(temp_dir, "convert_to_plink2")
    if (input_type %in% c("vcf", "bcf")) {
      convert_to_plink2(input$markerFileFilter$datapath, original_name = NULL, isplink = FALSE, name = pgen_prefix)
    } else {
      bed_prefix <- tools::file_path_sans_ext(input$bedFileFilter$datapath)
      convert_to_plink2(bed_prefix, original_name = NULL, isplink = TRUE, name = pgen_prefix)
    }

    incProgress(0.6, detail = "Filtering...")
    # for plink filtering
    plink_cmds <- c(shQuote(plink2_path), "--pfile", shQuote(pgen_prefix), "--recode", "vcf", "--out", file.path(temp_dir, "filtered"))

    if (input$filterIndiv) {
      plink_cmds <- c(plink_cmds, "--mind", input$mindThresh)
    }
    if (input$filterVariant) {
      plink_cmds <- c(plink_cmds, "--geno", input$genoThresh)
    }
    if (input$filterAllele) {
      plink_cmds <- c(plink_cmds, "--maf", input$mafThresh)
    }
    if (input$filterQuality) {
      plink_cmds <- c(plink_cmds, "--qual-threshold", input$qualThresh)
    }
    if (input$filterHWE) {
      plink_cmds <- c(plink_cmds, "--hwe", input$qualHWE, input$kval)
    }
    if (input$filterLD) {
      plink_cmds <- c(plink_cmds, "--indep-pairwise", input$ldWindow, input$ldStep, input$ldR2)
    }
    if (input$cutoffKing) {
      plink_cmds <- c(plink_cmds, "--king-cutoff", input$kingThresh)
    }

    # for custom flags
    custom_flag <- strsplit(trimws(input$customFilter), "\\s+")[[1]]
    if (length(custom_flag) == 1 && custom_flag[1] == "") custom_flag <- NULL

    if (!is.null(custom_flag) && length(custom_flag) >= 2) {
      if (!is.null(input$extraFile1) && length(custom_flag) >= 2) {
        custom_flag[2] <- shQuote(input$extraFile1$datapath)
      }

      if (!is.null(input$extraFile2) && length(custom_flag) >= 4) {
        custom_flag[4] <- shQuote(input$extraFile2$datapath)
      }
      plink_cmds <- c(plink_cmds, custom_flag)
    }

    system(paste(plink_cmds, collapse = " "))

    filtered_path <- file.path(temp_dir, "filtered.vcf")
    if (file.exists(filtered_path)) {
      filtered_plink_file(filtered_path)
    }
})
    output$plinkCommandPreview <- renderText({
      paste("plink", paste(plink_cmds, collapse = " "))
    })
    
  }) # end of observe events

  output$depthMarkerPlot <- renderImage(
    {
      req(depth_outputs())
      list(src = depth_outputs()$plot_marker, contentType = "image/png", width = "100%")
    },
    deleteFile = FALSE
  )

  output$depthSamplePlot <- renderImage(
    {
      req(depth_outputs())
      list(src = depth_outputs()$plot_sample, contentType = "image/png", width = "100%")
    },
    deleteFile = FALSE
  )

  output$downloadFilteredFile <- downloadHandler(
    filename = function() {
      "filtered.vcf"
    },
    content = function(file) {
      req(filtered_plink_file())
      file.copy(filtered_plink_file(), file)
    }
  )

  output$downloadDepthMarkerPlot <- downloadHandler(
    filename = function() {
      "Depth_marker.png"
    },
    content = function(file) {
      req(depth_outputs())
      file.copy(depth_outputs()$plot_marker, file)
    }
  )

  output$downloadDepthSamplePlot <- downloadHandler(
    filename = function() {
      "Depth_samples.png"
    },
    content = function(file) {
      req(depth_outputs())
      file.copy(depth_outputs()$plot_sample, file)
    }
  )

  output$depthMarkerPlot_UI <- renderUI({
    req(depth_outputs()$plot_marker)
    downloadButton("downloadDepthMarkerPlot", "Download Marker Depth Plot")
  })

  output$depthSamplePlot_UI <- renderUI({
    req(depth_outputs()$plot_sample)
    downloadButton("downloadDepthSamplePlot", "Download Sample Depth Plot")
  })

  output$downloadFilteredFile_UI <- renderUI({
    req(filtered_plink_file())
    downloadButton("downloadFilteredFile", "Download Filtered File")
  })
}
