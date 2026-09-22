snp_extraction_server <- function(input, output, session, rv) {
  # ================ MARKER EXTRACTION ====================#

  output$examplersID <- renderTable({
    data.frame(
      rsID = c("rs101", "rs102", "rs103", "rs104", "...")
    )
  })

  output$examplePOS <- renderTable({
    data.frame(
      rsID = c("rs01", "rs04", "..."),
      chromosome = c("1", "2", "..."),
      position = c("104500", "205300", "...")
    )
  })

  output$exampleTable <- renderTable({
    data.frame(
      Ind = c("sample1", "sample2", "sample3", "..."),
      rs101 = c("A/A", "A/T", "T/T", "..."),
      rs102 = c("G/C", "G/C", "G/G", "..."),
      rs103 = c("C/C", "C/G", "G/G", "..."),
      rs_n = c("...", "...", "...", "...")
    )
  })

  variant_ids <- reactiveVal(NULL)
  rsID_available <- reactiveVal(FALSE)
  extracted_file <- reactiveVal(NULL)

  observe({
    toggleState("validateBtn", (!is.null(input$markerFile) || !is.null(input$zippedPLINK) ||
      (!is.null(input$bedFile) && !is.null(input$bimFile) && !is.null(input$famFile))))
  })

  plink2_path <- get_plink2_path()

  observeEvent(input$validateBtn, {
    if (!is.null(input$zippedPLINK)) {
      plink_prefixes <- prepare_input_dataset(
        input_file = input$zippedPLINK$datapath,
        output.dir = tempdir()
      )
      plink_prefixes <- tools::file_path_sans_ext(plink_prefixes)
    } else if (!is.null(input$bedFile)) {
      plink_prefixes <- tools::file_path_sans_ext(input$bedFile$datapath)
    }

    temp_snplist <- file.path(tempdir(), "temp_snplist.txt")
    input_type <- if (!is.null(input$markerFile)) {
      if (grepl("\\.bcf$", input$markerFile$name, ignore.case = TRUE)) {
        "bcf"
      } else {
        "vcf"
      }
    } else {
      "bfile"
    }

    snps <- character(0)

    if (input_type %in% c("vcf", "bcf")) {
      cmd <- paste(
        shQuote(plink2_path),
        paste0("--", input_type), shQuote(input$markerFile$datapath),
        "--write-snplist allow-dups",
        "--out", shQuote(temp_snplist)
      )
      system(cmd)
    } else {
      cmd <- paste(
        shQuote(plink2_path),
        paste0("--", input_type), shQuote(plink_prefixes),
        "--write-snplist allow-dups",
        "--out", shQuote(temp_snplist)
      )
      system(cmd)
    }

    snps <- tryCatch(readLines(paste0(temp_snplist, ".snplist")),
      error = function(e) character(0)
    )

    snps <- unique(snps)

    if (length(snps) == 0) snps <- "."

    variant_ids(unique(snps))
    rsID_available(length(snps) > 1 || snps[1] != ".")
  })

  output$markerOptionsUI <- renderUI({
    req(variant_ids())

    if (rsID_available()) {
      tagList(
        radioButtons("markerType", "Choose Marker Type",
          choices = c("rsID", "pos"), inline = TRUE
        ),
        conditionalPanel(
          condition = "input.markerType == 'rsID'",
          radioButtons("rsIDInputType", "rsID Input", choices = c("manual", "upload")),
          conditionalPanel(
            condition = "input.rsIDInputType == 'manual'",
            textAreaInput("typedrsIDs", "Enter rsIDs (one per line)", rows = 5)
          ),
          conditionalPanel(
            condition = "input.rsIDInputType == 'upload'",
            fileInput("markerList1", "Upload rsID List File")
          )
        ),
        conditionalPanel(
          condition = "input.markerType == 'pos'",
          fileInput("markerList2", "Upload POS List (.csv, .xlsx)"),
          checkboxInput("addrsID", "Add marker information/rsID to output?", value = FALSE)
        ),
        shinyjs::disabled(actionButton("extractBtn", "Run Marker Extraction", icon = icon("play")))
      )
    } else {
      tagList(
        h4("No rsIDs detected. Extraction will require a POS list."),
        fileInput("markerList2", "Upload POS List (.csv, .xlsx)"),
        checkboxInput("addrsID", "Add marker information/rsID to output?", value = FALSE),
        shinyjs::disabled(actionButton("extractBtn", "Run Marker Extraction", icon = icon("play")))
      )
    }
  })

  output$variantTable <- DT::renderDT({
    snps <- variant_ids()
    req(!is.null(snps))
    display_snps <- if (all(snps == ".")) character(0) else snps
    DT::datatable(data.frame(Variant_ID = display_snps),
      options = list(scrollX = TRUE, pageLength = 10)
    )
  })

  can_extract <- reactive({
    if (!rsID_available()) {
      !is.null(input$markerList2)
    } else {
      if (is.null(input$markerType)) {
        return(FALSE)
      }
      if (input$markerType == "rsID") {
        if (is.null(input$rsIDInputType)) {
          return(FALSE)
        }
        if (input$rsIDInputType == "manual") {
          nzchar(trimws(input$typedrsIDs))
        } else if (input$rsIDInputType == "upload") {
          !is.null(input$markerList1)
        } else {
          FALSE
        }
      } else if (input$markerType == "pos") {
        !is.null(input$markerList2)
      } else {
        FALSE
      }
    }
  })

  observe({
    shinyjs::toggleState("extractBtn", condition = can_extract())
  })

  observeEvent(input$extractBtn, {
    disable("extractBtn")
    temp_dir <- tempdir()

    tryCatch(
      {
        pgen_prefix <- file.path(temp_dir, "input_pgen")
        print(input$markerFile$name)
        print(input$markerFile$datapath)

        if (!is.null(input$markerFile)) {
          input_file <- input$markerFile$datapath
          convert_to_plink2(input_file,
            original_name = input$markerFile$name,
            isplink = FALSE,
            name = pgen_prefix
          )
        } else {
          bed_prefix <- tools::file_path_sans_ext(input$bedFile$datapath)
          convert_to_plink2(bed_prefix,
            original_name = NULL,
            isplink = TRUE, name = pgen_prefix
          )
        }

        merged_name <- "extracted_markers"

        if (rsID_available() && input$markerType == "rsID") {
          snps_list <- tempfile(fileext = ".txt")
          if (input$rsIDInputType == "manual") {
            rsIDs <- trimws(unlist(strsplit(input$typedrsIDs, "\n")))
            rsIDs <- rsIDs[nzchar(rsIDs)]
            writeLines(rsIDs, snps_list)
          } else {
            df <- load_csv_xlsx_files(input$markerList1$datapath)
            writeLines(as.character(df[[1]]), snps_list)
          }

          extracted <- extract_by_ID_pgen(
            pgen_prefix = pgen_prefix,
            snps_list = snps_list,
            output_dir = temp_dir,
            merged_file = merged_name
          )
        } else {
          req(input$markerList2)
          pos_list <- as.data.frame(load_csv_xlsx_files(input$markerList2$datapath))

          if (ncol(pos_list) < 2) {
            stop("Position file must contain: chr and pos columns")
          }

          if (isTRUE(input$addrsID)) {
            if (ncol(pos_list) < 3) {
              stop("Position file must contain: rsID, chr, and pos columns")
            }

            extracted <- extract_POStoID_pgen(
              pos_list = pos_list,
              pgen_prefix = pgen_prefix,
              output_dir = temp_dir
            )
          } else {
            extracted <- extract_by_pos_pgen(
              pos_list = pos_list,
              pgen_prefix = pgen_prefix,
              output_dir = temp_dir,
              merged_file = merged_name
            )
          }
        }
        extracted_file(extracted)
        showNotification(
          "Marker extraction completed successfully.",
          type = "message"
        )
      },
      error = function(e) {
        showNotification(
          paste("Extraction error:", e$message),
          type = "error",
          duration = 10
        )
      },
      finally = {
        enable("extractBtn")
      }
    )
  })

  output$downloadVCF <- downloadHandler(
    filename = function() {
      paste0("extracted_markers_", Sys.Date(), ".vcf")
    },
    content = function(file) {
      req(extracted_file())
      file.copy(extracted_file(), file, overwrite = TRUE)
    }
  )

  output$downloadVCF_UI_Extracted <- renderUI({
    req(extracted_file())
    downloadButton("downloadVCF", "Download Extracted File (VCF)")
  })


  # =============== CONCORDANCE ANALYSIS ==================#

  observe({
    toggleState("compareBtn", !is.null(input$concordanceFile1) && !is.null(input$concordanceFile2))
  })

  concordanceResult <- reactiveVal(NULL)
  concordancePlotPath <- reactiveVal(NULL)

  observeEvent(input$compareBtn, {
    disable("compareBtn")

    withProgress(message = "Analyzing files...", value = 0, {
      tryCatch(
        {
          req(input$concordanceFile1$datapath, input$concordanceFile2$datapath)

          phased_flag <- input$isPhased
          
          f1_ext <- tools::file_ext(input$concordanceFile1$name)
          if (f1_ext == "vcf") {
            file1 <- vcf_to_csv(input$concordanceFile1$datapath)
            file1 <- file1$with_meta
          } else if (f1_ext %in% c("csv", "xlsx")) {
            file1 <- load_csv_xlsx_files(input$concordanceFile1$datapath)
          }
          
          f2_ext <- tools::file_ext(input$concordanceFile2$name)
          if (f2_ext == "vcf") {
            file2 <- vcf_to_csv(input$concordanceFile2$datapath)
            file2 <- file2$with_meta
          } else if (f2_ext %in% c("csv", "xlsx")) {
            file2 <- load_csv_xlsx_files(input$concordanceFile2$datapath)
          }
          
          result <- calc_concordance(file1, file2, phased = phased_flag)
          plot <- plot_concordance(result)
          concordanceResult(plot$results)
          concordancePlotPath(plot$plot)
          showNotification("Concordance analysis complete, rendering outputs.", type = "message", duration = 30)
          print(Sys.time())
        },
        error = function(e) {
          showNotification(paste("Error during analysis:", e$message), type = "error", duration = 10)
        }
      )
    })
    
    enable("compareBtn")
  })

  output$concordanceResults <- DT::renderDataTable(
    {
      req(concordanceResult())
      concordanceResult()
    },
    options = list(
      scrollX = TRUE,
      pageLength = 10
    )
  )

  output$concordancePlot <- renderPlot({
    req(concordancePlotPath())
    concordancePlotPath()
  })

  output$downloadConcordance <- downloadHandler(
    filename = function() {
      "concordance.csv"
    },
    content = function(file) {
      readr::write_csv(concordanceResult(), file)
    }
  )

  output$downloadConcordancePlot <- downloadHandler(
    filename = function() {
      paste0("concordance_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      req(concordancePlotPath())
      ggsave(file, plot = concordancePlotPath(), width = 8, height = 6, dpi = 600)
    }, contentType = "image/png"
  )

  output$downloadConcordance_UI <- renderUI({
    req(concordanceResult())
    downloadButton("downloadConcordance", "Download Concordance Results")
  })

  output$downloadConcordancePlot_UI <- renderUI({
    req(concordancePlotPath())
    downloadButton("downloadConcordancePlot", "Download Concordance Plot")
  })
}
