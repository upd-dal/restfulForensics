file_conversion_server <- function(input, output, session, rv) {
  # ================= FILE CONVERSION =====================#
  convertedVCF <- reactiveVal(NULL)
  convertedFASTA <- reactiveVal(NULL)
  convertedtoCSV <- reactiveVal(NULL)
  convertedPLINK <- reactiveVal(NULL)
  convertedBreakdown <- reactiveVal(NULL)

  observe({
    hasfile <- !is.null(input$VCFFile) || !is.null(input$BCFFile) || !is.null(input$CSVFile) || (!is.null(input$bedFile) && !is.null(input$bimFile) && !is.null(input$famFile))

    breakdown_selected <- !is.null(input$breakdown_vcf) || !is.null(input$breakdown_bcf) || !is.null(input$breakdown_plink)
    breakdown_ready <- !breakdown_selected ||
      (!is.null(input$breakdown_column_vcf) || !is.null(input$breakdown_column_bcf) || !is.null(input$breakdown_column_plink))

    toggleState("ConvertFILES", condition = hasfile && breakdown_ready)
  })

  marker_info_format <- data.frame(
    SNP = c("rs01", "rs02", "rs03", "rs04", "..."),
    chromosome = c("chr1", "chr4", "chr5", "chr5", "..."),
    position = c("1004", "90986", "5768", "9384982", "..."),
    genetic_distance = c("0", "0", "0", "0", "..."),
    ref_allele = c("A", "T", "G", "G", "C"),
    alt_allele = c("T", "A", "C", "C", "G")
  )

  output$markerInfoFormat <- DT::renderDataTable(
    {
      req(marker_info_format)
      marker_info_format
    },
    options = list(
      scrollX = TRUE,
      pageLength = 5
    )
  )

  output.dir <- tempdir()
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  outputName <- paste0("converted_", timestamp, ".csv")

  output_type <- function(input) {
    switch(input$inputType1,
      "vcf1" = input$inputType2_vcf,
      "bcf1" = input$inputType2_bcf,
      "plink1" = input$inputType2_plink,
      "csv1" = "vcf2"
    )
  }

  observeEvent(input$ConvertFILES, {
    req(input$ConvertFILES)
    disable("ConvertFILES")
    outputType <- output_type(input)
    req(outputType)

    input_file <- switch(input$inputType1,
      "vcf1" = input$VCFFile$datapath,
      "bcf1" = input$BCFFile$datapath,
      "csv1" = input$CSVFile$datapath,
      "plink1" = input$bedFile$datapath
    )

    if (input$inputType1 == "csv1") {
      csv_ext <- tools::file_ext(input$CSVFile$name)
      if (csv_ext %in% c("zip", "tar")) {
        unpacked_files <- unpack_input_file(input$CSVFile$datapath, output.dir)
        file_names <- unpacked_files$data_files
        all.list <- list()

        for (x in file_names) {
          all.list[[x]] <- load_csv_xlsx_files(x)
        }
        csv_merged <- data.table::rbindlist(all.list, fill = TRUE)
      }

      csv_file <- if (csv_ext %in% c("zip", "tar")) {
        csv_merged
      } else {
        input_file
      }

      csv_to_gen_obj <- csv_to_gentibble(csv_file, loci.meta = input$lociMetaFile$datapath)
      vcf_file <- file.path(tempdir(), paste0("csv_to_vcf_", timestamp, ".vcf"))
      converted_file <- tidypopgen::gt_as_vcf(csv_to_gen_obj, file = vcf_file, overwrite = TRUE)
      convertedVCF(converted_file)
    } else {
      prepared <- convert_files(
        input_file = input_file,
        output.dir = output.dir
      )

      result <- convert_from_plink2(
        prefix = prepared,
        output_type = outputType,
        output.dir = output.dir,
        ref = NULL
      )

      if (outputType == "vcf2") {
        convertedVCF(result)
      }

      if (outputType %in% c("plink2", "plink1")) {
        convertedPLINK(result)
      }
      
      if (outputType == "csv2") {
        convertedtoCSV(result$with_meta)
      }
      
    }
    enable("ConvertFILES")
  })

  output$toCSVtable <- DT::renderDataTable(
    {
      req(convertedtoCSV())
      convertedtoCSV()
    },
    options = list(
      scrollX = TRUE,
      pageLength = 10
    )
  )
  
  output$downloadConvertedVCF <- downloadHandler(
    filename = function() {
      paste0("csv_to_vcf_", timestamp, ".vcf")
    },
    content = function(file) {
      req(convertedVCF())
      file.copy(convertedVCF(), file)
    }
  )
  
  output$downloadConvertedtoCSV <- downloadHandler(
    filename = function() {
      paste0("file_to_csv_", timestamp, ".csv")
    },
    content = function(file) {
      req(convertedtoCSV())
      readr::write_csv(convertedtoCSV(), file)
    }
  )
  
  output$downloadConvertedPLINK <- downloadHandler(
    filename = function() {
      if (output_type(input) == "plink1") {
        paste0("plink1_files_", timestamp, ".zip")
      } else {
        paste0("plink2_files_", timestamp, ".zip")
      }
    },
    content = function(file) {
      req(convertedPLINK())
      file.copy(convertedPLINK(), file)
    },
    contentType = "application/zip"
  )

  output$downloadVCF_UI <- renderUI({
    req(convertedVCF())
    downloadButton("downloadConvertedVCF", "Download VCF File")
  })

  output$downloadPLINK_UI <- renderUI({
    req(convertedPLINK())
    downloadButton("downloadConvertedPLINK", "Download PLINK File")
  })
  
  output$downloadtoCSV_UI <- renderUI({
    req(convertedtoCSV())
    downloadButton("downloadConvertedtoCSV", "Download CSV File")
  })
  
  

  # =================== Add Metadata ======================#
  exampleRefCSV <- data.frame(
    Sample.Name = c("sample1", "sample2", "sample3", "sample4", "..."),
    Population = c("Malaysia", "Mexico", "Greece", "South Korea", "..."),
    Superpopulation = c("Southeast Asia", "North and South America", "Europe", "East Asia", "...")
  )

  output$ExampleRefFile <- DT::renderDataTable(
    {
      req(exampleRefCSV)
      exampleRefCSV
    },
    options = list(
      scrollX = TRUE,
      pageLength = 5
    )
  )

  convertedCSV <- reactiveVal(NULL)
  missingData <- reactiveVal(NULL)
  convertedBreakdown <- reactiveVal(NULL)
  output.dir <- tempdir()
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

  observe({
    fileready <- !is.null(input$genotypeFile)
    singlePop <- nzchar(input$typePop_meta)
    multiPop <- !is.null(input$refMetadata)
    metaReady <- singlePop || multiPop

    toggleState("addMetadata", condition = fileready && metaReady)
  })

  columns_target <- reactive({
    req(input$refMetadata)
    load_csv_xlsx_files(input$refMetadata$datapath)
  })

  output$metaHeader <- renderUI({
    req(columns_target())
    checkboxGroupInput(
      inputId = "col_targets",
      label = "Choose columns to be merged with sample and genotype data",
      choices = names(columns_target()),
      selected = NULL
    )
  })

  column_for_breakdown <- reactive({
    req(input$col_targets)
    input$col_targets
  })

  output$selectMetaHeader <- renderUI({
    req(input$col_targets)
    req(input$breakdownPop == "YesBreakdown")
    selectInput(
      "popForBreakdown",
      label = "Choose column as basis for tally",
      choices = input$col_targets
    )
  })

  output.dir <- tempdir()

  observeEvent(input$addMetadata, {
    disable("addMetadata")

    tryCatch(
      {
        input_file <- input$genotypeFile$datapath

        if (!is.null(input$refMetadata)) {
          ref_file <- load_csv_xlsx_files(input$refMetadata$datapath)

          if (length(input$col_targets) <= 0) {
            stop("Select at least one column to merge")
          }
          
          id_col <- colnames(ref_file)[1]
          selected_cols <- unique(c(id_col, input$col_targets))
          for_merging <- data.frame(ref_file[
            ,
            selected_cols,
            drop = FALSE
          ])
        } else {
          for_merging <- input$typePop_meta
        }
        
        # Check extension
        file_extension <- tools::file_ext(input_file)
        
        if (file_extension %in% c("zip", "tar")) {
          # Unpack to determine the data type
          unpacked <- unpack_input_file(input_file, output.dir)
          files <- unpacked$data_files
          ext <- tools::file_ext(files[[1]])
          
          if (ext == "csv") {
            all.list <- list()
            
            for (x in files) {
              all.list[[x]] <- read.csv(x, check.names = FALSE, row.names = 1)
            }
            
            merged <- dplyr::bind_rows(all.list)
            merged <- data.frame(rownames(merged), merged)
            result <- add_metadata(merged, for_merging)
          } else {
            result <- convert_files(input_file, output.dir = output.dir) # returns plink 2.0 files
            
            # convert to vcf
            vcf_file <- convert_from_plink2(result,
                                            output_type = "vcf2",
                                            output.dir = output.dir)
            # vcf to csv and merge
            csv_file <- convert_from_plink2(vcf_file, 
                                            output_type = "csv2",
                                            output.dir = output.dir,
                                            ref = for_merging)
            result <- csv_file
          }
          
        } else {
          # if single files
          if (file_extension == "csv") {
            csv_df <- load_csv_xlsx_files(input_file)
            result <- add_metadata(csv_df, for_merging)
          } else {
            result <- convert_files(input_file, output.dir = output.dir) # returns plink 2.0 files
            
            # convert to vcf
            vcf_file <- convert_from_plink2(result,
                                            output_type = "vcf2",
                                            output.dir = output.dir)
            # vcf to csv and merge
            csv_file <- convert_from_plink2(vcf_file, 
                                            output_type = "csv2",
                                            output.dir = output.dir,
                                            ref = for_merging)
            result <- csv_file 
          }
        }

        convertedCSV(result$with_meta)
        missingData(result$missing)

        if (input$breakdownPop == "YesBreakdown") {

          if (length(input$popForBreakdown) <= 0) {
            stop("No column selected.")
          }

          breakdown_results <- pop_breakdown(convertedCSV(), input$popForBreakdown)
          convertedBreakdown(breakdown_results)
        }
      }, # end for try catch
      error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
      }
    )
    enable("addMetadata")
  })

  output$downloadConvertedCSV <- downloadHandler(
    filename = function() {
      paste0("converted_", timestamp, ".csv")
    },
    content = function(file) {
      req(convertedCSV())
      readr::write_csv(convertedCSV(), file)
    }
  )

  output$downloadPopBreakdown <- downloadHandler(
    filename = function() {
      paste0("pop_breakdown_", timestamp, ".csv")
    },
    content = function(file) {
      req(convertedBreakdown())
      readr::write_csv(convertedBreakdown(), file)
    }
  )

  output$downloadMissingMeta <- downloadHandler(
    filename = function() {
      paste0("missing_meta_", timestamp, ".csv")
    },
    content = function(file) {
      req(missingData())
      readr::write_tsv(missingData(), file)
    }
  )

  output$previewTable <- DT::renderDataTable(
    {
      req(convertedCSV())
      convertedCSV()
    },
    options = list(
      scrollX = TRUE,
      pageLength = 10
    )
  )

  output$previewTableBreakdown <- DT::renderDataTable(
    {
      req(convertedBreakdown())
      convertedBreakdown()
    },
    options = list(
      scrollX = TRUE,
      pageLength = 10
    )
  )

  output$previewMissingData <- DT::renderDataTable(
    {
      req(missingData())
      missingData()
    },
    options = list(
      scrollX = TRUE,
      pageLength = 10
    )
  )

  output$downloadCSV_UI <- renderUI({
    req(convertedCSV())
    downloadButton("downloadConvertedCSV", "Download CSV File")
  })

  output$downloadBreakdown <- renderUI({
    req(convertedBreakdown())
    downloadButton("downloadPopBreakdown", "Download Population Breakdown Count")
  })

  output$downloadMissing <- renderUI({
    req(missingData())
    downloadButton("downloadMissingMeta", "Download Samples without Metadata")
  })

  # ================== Widen GT file ====================#

  convertedUAS <- reactiveVal(NULL)
  observe({
    shinyjs::toggleState("run_uas2csv", !is.null(input$uas_zip))
  })

  output$exampleXLSX <- renderTable({
    data.frame(
      Sample.Name = c("sample1", "sample1", "sample1", "sample1", "sample1", "sample1", "sample2", "sample3", "sample3", "sample3", "..."),
      Locus = c("rs01", "rs01", "rs02", "rs02", "rs02", "rs03", "rs01", "rs01", "rs02", "rs03", "..."),
      Allele = c("A", "T", "C", "A", "G", "T", "A", "T", "G", "A", "...")
    )
  })

  temp_dir <- tempdir()
  observeEvent(input$run_uas2csv, {
    req(input$uas_zip)
    disable("run_uas2csv")
    input_path <- file.path(temp_dir, input$uas_zip$name)
    file.copy(input$uas_zip$datapath, input_path, overwrite = TRUE)

    withProgress(message = "Converting file...", value = 0, {
      tryCatch(
        {
          widened.file <- widen_genotype_file(
            files = input_path,
            #population = ref_value,
            output.dir = temp_dir
          )
          convertedUAS(widened.file)
          showNotification("Conversion complete!", type = "message")
        },
        error = function(e) {
          showNotification(paste("Error:", e$message), type = "error")
        }
      )
    })
    enable("run_uas2csv")
  })

  outputName <- "merged_typed_data.csv"

  output$downloadUAScsv <- downloadHandler(
    filename = function() {
      outputName
    },
    content = function(file) {
      readr::write_csv(convertedUAS(), file)
    }
  )

  output$downloadUAScsv_UI <- renderUI({
    req(convertedUAS())
    downloadButton("downloadUAScsv", "Download CSV File")
  })

  output$previewTableUAS <- DT::renderDataTable(
    {
      req(convertedUAS())
      convertedUAS()
    },
    options = list(
      scrollX = TRUE,
      pageLength = 10
    )
  )

  # ==================== To SNIPPER =======================#
  exampleTableSnipper1 <- data.frame(
    Ind = c("sample1", "sample2", "sample3", "sample4", "..."),
    rs101 = c("A/A", "A/T", "T/T", "A/T", "..."),
    rs102 = c("G/C", "G/C", "G/G", "G/C", "..."),
    rs103 = c("C/C", "C/G", "G/G", "G/G", "..."),
    rs_n = c("...", "...", "...", "...", "...")
  )

  output$exampleTableSnipper1 <- DT::renderDataTable(
    {
      req(exampleTableSnipper1)
      exampleTableSnipper1
    },
    options = list(
      scrollX = TRUE,
      pageLength = 10
    )
  )

  exampleTableSnipper2 <- data.frame(
    Ind = c("sample1", "sample2", "sample3", "sample4", "..."),
    Pop = c("sub_pop1", "sub_pop2", "sub_pop3", "sub_pop4", "..."),
    Superpop = c("region1", "region1", "region2", "region3", "..."),
    rs101 = c("A/A", "A/T", "T/T", "A/T", "..."),
    rs102 = c("G/C", "G/C", "G/G", "G/C", "..."),
    rs103 = c("C/C", "C/G", "G/G", "G/G", "..."),
    rs_n = c("...", "...", "...", "...", "...")
  )

  output$exampleTableSnipper2 <- DT::renderDataTable(
    {
      req(exampleTableSnipper2)
      exampleTableSnipper2
    },
    options = list(
      scrollX = TRUE,
      pageLength = 10
    )
  )

  exampleRefSnipper <- data.frame(
    Sample.Name = c("sample1", "sample2", "sample3", "sample4", "..."),
    Population = c("Malaysia", "Mexico", "Greece", "South Korea", "..."),
    Superpopulation = c("Southeast Asia", "North and South America", "Europe", "East Asia", "...")
  )

  output$exampleRefSnipper <- DT::renderDataTable(
    {
      req(exampleRefSnipper)
      exampleRefSnipper
    },
    options = list(
      scrollX = TRUE,
      pageLength = 10
    )
  )

  convertedSNIPPER <- reactiveVal(NULL)
  outputName <- "snipper.xlsx"

  observe({
    hasFile <- !is.null(input$convertFile) && nrow(input$convertFile) > 0
    ready <- isTRUE(hasFile)
    shinyjs::toggleState("convertBtn", ready)
  })

  observeEvent(input$convertBtn, {
    disable("convertBtn")

    tosnipper_file <- load_csv_xlsx_files(input$convertFile$datapath)
    tosnipper_file <- dplyr::rename(tosnipper_file, Sample = 1)

    if (!is.null(input$refProvided)) {
      inputPath <- tosnipper_file[, -c(2, 3)]
      refPath <- tosnipper_file[, 2:3]
    } else {
      inputPath <- tosnipper_file
      refPath <- load_csv_xlsx_files(input$refFile$datapath)
    }

    refPath <- dplyr::rename(refPath, Sample = 1)
    targetSet <- input$targetPop
    targetName <- if (targetSet) input$targetPopName else NULL
    inputData <- colnames(inputPath)
    numMarkers <- length(inputData) - 1

    withProgress(message = "Converting to SNIPPER-analysis ready file...", value = 0, {
      snipper.file <- tryCatch(
        {
          to_snipper(
            input = inputPath,
            references = refPath,
            target.pop = targetSet,
            population.name = targetName,
            markers = numMarkers
          )
        },
        error = function(e) {
          showNotification(paste("Conversion failed:", e$message), type = "error")
          NULL
        }
      )

      if (!is.null(snipper.file)) {
        convertedSNIPPER(snipper.file)
        enable("convertBtn")
      }
      enable("convertBtn")
    })
  })

  output$downloadConverted <- downloadHandler(
    filename = function() {
      outputName
    },
    content = function(file) {
      openxlsx::write.xlsx(convertedSNIPPER(), file)
    }
  )

  output$downloadSNIPPER <- renderUI({
    req(convertedSNIPPER())
    downloadButton("downloadConverted", "Download SNIPPER-ready file")
  })

  output$previewTableSNIPPER <- DT::renderDataTable(
    {
      req(convertedSNIPPER())
      convertedSNIPPER()
    },
    options = list(
      scrollX = TRUE,
      pageLength = 10
    )
  )

  # ==================== CSV to STR =======================#
  examplePop_STR <- data.frame(
    Sample = c("Sample1", "Sample2", "Sample3", "Sample4", "..."),
    Population = c("POP1", "POP2", "POP3", "POP4", "..."),
    rs101 = c("A/A", "A/T", "A/A", "T/T", "..."),
    rs102 = c("G/G", "C/C", "G/C", "G/G", "..."),
    rs_n = c("...", "...", "...", "...", "...")
  )

  output$examplePop_STRUI <- DT::renderDataTable(
    {
      req(examplePop_STR)
      examplePop_STR
    },
    options = list(
      scrollX = TRUE,
      pageLength = 5
    )
  )

  observe({
    toggleState("csv2str", !is.null(input$tostrFile) && !is.null(input$systemFile))
  })

  csv_revised <- reactiveVal(NULL)
  strconvert <- reactiveVal(NULL)
  str_file <- reactiveVal(NULL)
  directory <- tempdir()

  observeEvent(input$csv2str, {
    disable("csv2str")

    withProgress(message = "Analyzing files...", value = 0, {
      tryCatch({
        req(input$tostrFile$datapath, input$systemFile)
        csv_file <- load_csv_xlsx_files(input$tostrFile$datapath)
        csv_file <- clean_input_data(csv_file)
        genind <- convert_to_genind(csv_file, to_str = TRUE, popinfo = TRUE)
        csv_revised(genind$csv_revised)
        strconvert(genind$fsnps_gen)
        str_path <- revise_structure_file(strconvert(), directory, system = input$systemFile)
        str_file(str_path)
      }, error = function(e) {
        showNotification(paste("Error during STRUCTURE conversion", e$message), type = "error", duration = 20)
      }, finally = {
        enable("csv2str")
      })
    })
    shinyjs::enable("csv2str")
  })

  output$revisedCSV <- DT::renderDataTable(
    {
      req(csv_revised())
      csv_revised()
    },
    options = list(pageLength = 10, scrollX = TRUE)
  )

  output$strFile <- DT::renderDataTable(
    {
      req(str_file())
      str_file <- read.table(str_file(), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      DT::datatable(str_file)
    },
    options = list(pageLength = 10, scrollX = TRUE)
  )

  output$downloadrevised <- downloadHandler(
    filename = function() {
      "revised_input.csv"
    },
    content = function(file) {
      req(csv_revised())
      readr::write_csv(csv_revised(), file)
    }
  )

  output$downloadSTRfile <- downloadHandler(
    filename = function() {
      "structure_file.str"
    },
    content = function(file) {
      req(str_file())
      file.copy(str_file(), file)
    }
  )

  output$downloadrevised_UI <- renderUI({
    req(csv_revised())
    downloadButton("downloadrevised", "Download Revised CSV File")
  })

  output$downloadSTRfile_UI <- renderUI({
    req(str_file())
    downloadButton("downloadSTRfile", "Download .str File")
  })

  # ==================== To Arlequin-compatible file =======================#
  exampleForArlecore <- data.frame(
    Sample = c("Sample1", "Sample2", "Sample3", "Sample4", "..."),
    Population = c("POP1", "POP2", "POP3", "POP4", "..."),
    rs101 = c("A/A", "A/T", "A/A", "T/T", "..."),
    rs102 = c("G/G", "C/C", "G/C", "G/G", "..."),
    rs_n = c("...", "...", "...", "...", "...")
  )

  output$exampleForArlecore_UI <- DT::renderDataTable(
    {
      req(exampleForArlecore)
      exampleForArlecore
    },
    options = list(
      scrollX = TRUE,
      pageLength = 10
    )
  )

  arleFile <- reactiveVal(NULL)

  observe({
    shinyjs::toggleState("convert2Arle", !is.null(input$toArleFile))
  })

  observeEvent(input$convert2Arle, {
    disable("convert2Arle")

    for_arp <- load_csv_xlsx_files(input$toArleFile$datapath)
    for_arp <- clean_input_data(for_arp)

    # All null values are "N", set to ""
    for_arp <- for_arp %>%
      mutate(across(everything(), ~ case_when(
        . == "N" ~ "",
        TRUE ~ .x
      )))
    for_arp <- as.data.frame(for_arp)

    # create the arp file
    arp_file <- build_arp_per_population(for_arp,
      genotypic_data = as.numeric(input$genotypicData),
      gametic_phase = as.numeric(input$gameticPhase),
      recessive_data = as.numeric(input$recessiveData),
      locus_sep = input$locusSep,
      output.prefix = "arp_file",
      data_type = "STANDARD"
    )
    
    arleFile(arp_file)

    enable("convert2Arle")
  })

  output$downloadArpFile <- downloadHandler(
    filename = function() {
      "arp_file.arp"
    },
    content = function(file) {
      req(arleFile())
      file.copy(arleFile(), file)
    }
  )

  output$downloadArpFile_UI <- renderUI({
    req(arleFile())
    downloadButton("downloadArpFile", "Download Arlequin-compatible input file")
  })
}
