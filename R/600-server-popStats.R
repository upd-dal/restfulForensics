pop_stats_server <- function(input, output, session, rv) {
  # =============== POPULATION STATISTICS =================#

  examplePop <- data.frame(
    Sample = c("Sample1", "Sample2", "Sample3", "Sample4", "..."),
    Population = c("POP1", "POP2", "POP3", "POP4", "..."),
    rs101 = c("A/A", "A/T", "A/A", "T/T", "..."),
    rs102 = c("G/G", "C/C", "G/C", "G/G", "..."),
    rs_n = c("...", "...", "...", "...", "...")
  )

  output$examplePop_UI <- DT::renderDataTable(
    {
      req(examplePop)
      examplePop
    },
    options = list(
      scrollX = TRUE,
      pageLength = 5
    )
  )

  observe({
    file_ready <- !is.null(input$popStatsFile)
    shinyjs::toggleState("runPopStats", condition = file_ready)
  })

  privAlleles <- reactiveVal(NULL)
  popStats <- reactiveVal(NULL)
  hardyWeinberg <- reactiveVal(NULL)
  fstStats <- reactiveVal(NULL)
  fstData <- reactiveVal(NULL)
  afData <- reactiveVal(NULL)
  statsMatrix <- reactiveVal(NULL)
  hwMatrix <- reactiveVal(NULL)
  fstMatrix <- reactiveVal(NULL)

  fsnps_gen <- reactive({
    req(input$popStatsFile)
    df <- load_csv_xlsx_files(input$popStatsFile$datapath)
    cleaned <- clean_input_data(df)
    convert_to_genind(cleaned, to_str = FALSE, popinfo = TRUE)
  })

  observeEvent(input$runPopStats, {
    disable("runPopStats")
    req(input$popStatsFile)

    withProgress(message = "Running population analysis...", value = 0, {
      tryCatch(
        {
          req(fsnps_gen())

          incProgress(0.4, detail = "Computing private alleles...")
          priv_alleles <- poppr::private_alleles(fsnps_gen())
          if (is.null(priv_alleles)) {
            priv_alleles <- list(message = "No private alleles detected")
          } else {
            priv_alleles <- data.frame(rownames(priv_alleles), priv_alleles)
            priv_alleles <- dplyr::rename(priv_alleles, Pop = 1)
            rownames(priv_alleles) <- NULL
          }
          privAlleles(priv_alleles)

          incProgress(0.6, detail = "Computing population statistics...")
          population_stats <- compute_pop_stats(fsnps_gen())
          popStats(population_stats)

          ## AF
          af_stats <- compute_af(fsnps_gen())
          afData(af_stats)

          ## HWE
          incProgress(0.8, detail = "Running HWE and FST calculations...")
          hardy_weinberg_stats <- compute_hwe(fsnps_gen(), correction = input$correctionModel, alpha = input$alphaValue, chains = input$markovHWE)
          hardyWeinberg(hardy_weinberg_stats)

          ## FST
          incProgress(1.0, detail = "Still running HWE and FST calculations...")
          fst_stats <- compute_fst(fsnps_gen())
          fstStats(fst_stats)
          fst_data <- fstStats()$fst_dataframe
          fstData(fst_data)

          showNotification("Rendering outputs, this might take some time...", type = "message", duration = 30)
          print(Sys.time())

          enable("runPopStats")
        },
        error = function(e) {
          showNotification(paste("Population stats error:", e$message), type = "error")
          enable("runPopStats")
        }
      )

      stats_matrix <- compute_pop_stats(fsnps_gen())
      hw_matrix <- compute_hwe(fsnps_gen())
      fst_matrix <- compute_fst(fsnps_gen())
      statsMatrix(stats_matrix)
      hwMatrix(hw_matrix)
      fstMatrix(fst_matrix)
    })
  })

  output$privateAlleleTable <- DT::renderDataTable(
    {
      as.data.frame(privAlleles())
    },
    options = list(pageLength = 10, scrollX = TRUE)
  )

  output$meanallelic <- DT::renderDataTable({
    popStats()$mar_list
  })

  output$heterozygosity_table <- DT::renderDataTable({
    popStats()$heterozygosity
  })

  output$heterozygosity_plot <- renderImage(
    {
      req(popStats()$heterozygosity)

      plot_path <- plot_heterozygosity(
        Het_fsnps_df = popStats()$heterozygosity,
        out_dir = tempdir()
      )
      list(
        src = plot_path,
        contentType = "image/png",
        alt = "Heterozygosity Plot",
        width = "100%"
      )
    },
    deleteFile = TRUE
  )

  output$inbreeding_table <- DT::renderDataTable({
    popStats()$inbreeding_coeff
  })

  output$ttest_table <- DT::renderDataTable({
    popStats()$ttest
  })

  output$allele_freq_table <- DT::renderDataTable(
    {
      afData()
    },
    options = list(scrollX = TRUE)
  )


  output$hwe_summary_text <- DT::renderDataTable(
    {
      hardyWeinberg()$hw_summary
    },
    options = list(scrollX = TRUE)
  )

  output$hwe_chisq_table <- DT::renderDataTable(
    {
      hardyWeinberg()$hw_dataframe
    },
    options = list(scrollX = TRUE)
  )

  output$hwe_loci <- DT::renderDataTable(
    {
      hardyWeinberg()$loci_HWE_failure
    },
    options = list(scrollX = TRUE)
  )

  output$hwe_pop <- DT::renderDataTable(
    {
      hardyWeinberg()$pops_out_of_HWE
    },
    options = list(scrollX = TRUE)
  )

  output$fstMatrixUI <- renderUI({
    fst <- fstStats()$fst_matrix
    if (is.list(fst) && "message" %in% names(fst)) {
      tags$p(style = "color:gray;", fst$message)
    } else {
      DT::dataTableOutput("fstMatrixTable")
    }
  })

  output$fstMatrixTable <- DT::renderDataTable(
    {
      matrix_data <- matrix(unlist(fstStats()$fst_matrix),
        nrow = sqrt(length(fstStats()$fst_matrix)),
        byrow = TRUE
      )
      rownames(matrix_data) <- colnames(matrix_data) <- attr(fsnps_gen(), "pop.names")
      as.data.frame(matrix_data)
    },
    options = list(scrollX = TRUE)
  )

  output$fstDfTable <- DT::renderDataTable({
    fstStats()$fst_dataframe
  })

  output$fst_heatmap_plot <- renderImage(
    {
      req(fstData())

      plot_path <- plot_fst(
        fst_df = fstData(),
        out_dir = tempdir()
      )
      list(
        src = plot_path,
        contentType = "image/png",
        alt = "FST Heatmap",
        width = "100%"
      )
    },
    deleteFile = TRUE
  )

  output$downloadHeterozygosityPlot <- downloadHandler(
    filename = function() {
      "heterozygosity_plot.png"
    },
    content = function(file) {
      plot_path <- plot_heterozygosity(
        Het_fsnps_df = popStats()$heterozygosity,
        out_dir = tempdir()
      )
      file.copy(plot_path, file)
    }
  )

  output$downloadHeterozygosityPlot_UI <- renderUI({
    downloadButton("downloadHeterozygosityPlot", "Download Heterozygosity Plot")
  })

  output$downloadFstHeatmap <- downloadHandler(
    filename = function() {
      "fst_heatmap.png"
    },
    content = function(file) {
      plot_path <- plot_fst(
        fst_df = fstData(),
        out_dir = tempdir()
      )
      file.copy(plot_path, file)
    }
  )

  output$downloadFstHeatmap_UI <- renderUI({
    downloadButton("downloadFstHeatmap", "Download Fst Plot")
  })

  ## download all results
  output$downloadStatsXLSX <- downloadHandler(
    filename = function() {
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
      paste0("population-statistics-results_", timestamp, ".xlsx")
    },
    content = function(file) {
      req(statsMatrix(), hwMatrix(), fstMatrix(), privAlleles(), afData())
      path <- export_pop_results(
        allele_freq = afData(),
        priv_alleles = privAlleles(),
        stats_matrix = statsMatrix(),
        hw_matrix = hwMatrix(),
        fst_matrix = fstMatrix(), dir = tempdir()
      )

      file.copy(path, file)
    }
  )

  output$downloadStatsXLSX_UI <- renderUI({
    req(statsMatrix(), hwMatrix(), fstMatrix(), privAlleles())
    downloadButton("downloadStatsXLSX", "Download Results (excel)")
  })

  # ==================== To Arlequin-compatible file =======================#
  exampleForArlecore <- data.frame(
    Sample = c("Sample1", "Sample2", "Sample3", "Sample4", "..."),
    Population = c("POP1", "POP2", "POP3", "POP4", "..."),
    rs101 = c("A/A", "A/T", "A/A", "T/T", "..."),
    rs102 = c("G/G", "C/C", "G/C", "G/G", "..."),
    rs_n = c("...", "...", "...", "...", "...")
  )

  output$exampleForArlecore <- DT::renderDataTable(
    {
      req(exampleForArlecore)
      exampleForArlecore
    },
    options = list(
      scrollX = TRUE,
      pageLength = 10
    )
  )

  arlequinFile <- reactiveVal(NULL)
  arlequinResult <- reactiveVal(NULL)
  arlequinPopLabels <- reactiveVal(NULL)
  arlequinHeterozygosity <- reactiveVal(NULL)
  arlequinFstMatrix <- reactiveVal(NULL)
  arlequinCoancestry <- reactiveVal(NULL)
  arlequinPairwise <- reactiveVal(NULL)
  arlequinPopDiversity <- reactiveVal(NULL)
  arlequinLD <- reactiveVal(NULL)
  locusData <- reactiveVal(NULL)

  observe({
    shinyjs::toggleState("runArlecore", !is.null(input$fileForArlecore))
  })

  observeEvent(input$runArlecore, {
    disable("runArlecore")

    withProgress(message = "Analysis ongoing...", {
      incProgress(0.2, detail = "Loading input file...")
      for_arp <- load_csv_xlsx_files(input$fileForArlecore$datapath)
      for_arp <- clean_input_data(for_arp)

      rsids <- as.data.frame(colnames(for_arp)[-c(1, 2)])
      rsids <- data.frame(rownames(rsids), rsids)
      locusData(rsids)

      # All null values are "N", set to ""
      for_arp <- for_arp %>%
        mutate(across(everything(), ~ case_when(
          . == "N" ~ "",
          TRUE ~ .x
        )))
      for_arp <- as.data.frame(for_arp)

      incProgress(0.4, detail = "Creating input file...")
      arp_file <- build_arp_per_population(for_arp,
        genotypic_data = as.numeric(input$genotypicData),
        gametic_phase = as.numeric(input$gameticPhase),
        recessive_data = as.numeric(input$recessiveData),
        locus_sep = input$locusSep,
        output.prefix = "arp_file",
        data_type = "STANDARD"
      )
      arlequinFile(arp_file)

      incProgress(0.6, detail = "Running arlecore...")
      # returns path of res folder
      ld_value <- input$calcLD
      hwe_value <- input$calcHWE
      print(ld_value)
      print(hwe_value)
      results <- run_arlequin(arp_file, ld = ld_value, hwe = hwe_value)

      incProgress(0.8, detail = "Loading report...")
      # get report file
      file_dir <- dirname(results)
      xml_file <- file.path(results, "arp_file.xml")

      if (!file.exists(xml_file)) {
        stop("Resulting 'arp_file.xml' report not found.")
      }

      arlequinResult(xml_file)

      doc <- xml2::read_xml(xml_file)

      # run it separately
      pop_labels <- parse_pop_labels(doc)
      heterozygosity <- parse_sections_arlequin(doc, "sumExpHeterozygosity", sumExpectedHeterozygosity) # bar plot
      hwe <- heterozygosity[["data"]]
      hwe$Locus <- rsids[[2]][
        match(as.numeric(hwe$Locus), as.numeric(rsids[[1]]))
      ]
      heterozygosity[["data"]] <- hwe

      fst_matrix <- parse_sections_arlequin(doc, "PairFstMat", pairFstMatrix) # matrix of pairwise fst
      coancestry_coeff <- parse_sections_arlequin(doc, "coancestryCoefficients", coancestryCoeff) # pairwise of fst and reynolds
      pairwise_matrix <- parse_sections_arlequin(doc, "pairwiseDifferenceMatrix", pairwiseDiffMatrix)

      
      if (isTRUE(input$calcHWE)) {
        population_diversity <- parse_pop_diversity(doc)
        population_diversity$Locus <- rsids[[2]][
          match(as.numeric(population_diversity$Locus), as.numeric(rsids[[1]]))
        ]
      } else {
        population_diversity <- NULL
      }

      if (isTRUE(input$calcLD)) {
        rsids_zero <- rsids
        rsids_zero[[1]] <- seq(0, nrow(rsids_zero) - 1)
        ld_vals <- parse_ld(doc)
        print("LD parsing")
        print(str(ld_vals))
        print(dim(ld_vals))
        print(head(ld_vals))
        ld_vals$Locus1 <- rsids_zero[[2]][
          match(ld_vals$Locus1, rsids_zero[[1]])
        ]
        ld_vals$Locus2 <- rsids_zero[[2]][
          match(ld_vals$Locus2, rsids_zero[[1]])
        ]
      } else {
        ld_vals <- NULL
      }

      arlequinPopLabels(pop_labels)
      arlequinHeterozygosity(heterozygosity)
      arlequinFstMatrix(fst_matrix)
      arlequinCoancestry(coancestry_coeff)
      arlequinPairwise(pairwise_matrix)
      arlequinPopDiversity(population_diversity)
      arlequinLD(ld_vals)
    })

    enable("runArlecore")
  }) # end of observe event

  # separate the tables
  output$population_tables <- renderUI({
    req(arlequinPopDiversity())
    populations <- unique(arlequinPopDiversity()$Population)

    tagList(
      lapply(
        seq_along(populations),
        function(i) {
          population <- populations[i]

          tagList(
            h3(population),
            DT::DTOutput(paste0("diversity_table_", i)),
            br()
          )
        }
      )
    )
  })

  observe({
    req(arlequinPopDiversity())
    data <- arlequinPopDiversity()
    populations <- unique(data$Population)

    lapply(seq_along(populations), function(i) {
      population <- populations[i]
      output_id <- paste0("diversity_table_", i)

      local({
        pop <- population
        output[[output_id]] <- DT::renderDT({
          pop_data <- arlequinPopDiversity() %>%
            dplyr::filter(Population == pop) %>%
            dplyr::select(-Population)

          DT::datatable(pop_data,
            rownames = FALSE,
            options = list(
              pageLength = 10,
              scrollX = TRUE
            )
          )
        })
      })
    })
  })

  output$hwe_arlecore <- DT::renderDataTable(
    {
      req(arlequinHeterozygosity())
      arlequinHeterozygosity()[["data"]]
    },
    options = list(scrollX = TRUE)
  )

  output$hwe_arlecore_plots <- renderUI({
    hwe <- arlequinHeterozygosity()[["data"]]
    req(hwe)
    populations <- setdiff(names(hwe), "Locus")

    tagList(
      lapply(seq_along(populations), function(i) {
        tagList(
          h3(populations[i]),
          plotly::plotlyOutput(
            paste0("het_plot_", i),
            height = "400px"
          ),
          br()
        )
      })
    )
  })

  observe({
    hwe <- arlequinHeterozygosity()[["data"]]
    req(hwe)
    populations <- setdiff(names(hwe), "Locus")

    lapply(seq_along(populations), function(i) {
      local({
        pop <- populations[i]
        output_id <- paste0("het_plot_", i)

        output[[output_id]] <- plotly::renderPlotly({
          pop_data <- hwe %>%
            dplyr::select(Locus, Heterozygosity = dplyr::all_of(pop))

          plotly::plot_ly(
            data = pop_data,
            x = ~Locus,
            y = ~Heterozygosity,
            type = "bar"
          ) %>%
            plotly::layout(
              title = list(text = pop),
              xaxis = list(title = "Locus"),
              yaxis = list(title = "Expected Heterozygosity"),
              showlegend = FALSE
            )
        })
      })
    })
  })

  output$hwe_arlecore_plot <- plotly::renderPlotly({
    data <- arlequinHeterozygosity()[["data"]]
    req(data)
    data_long <- data %>% tidyr::pivot_longer(
      cols = -1,
      names_to = "Population",
      values_to = "Heterozygosity"
    )
    fig <- plotly::plot_ly(
      data = data_long,
      x = ~Locus,
      y = ~Heterozygosity,
      type = "bar",
      color = ~Population,
      colors = "Set2",
      text = ~Population
    ) %>%
      plotly::layout(
        yaxis = list(title = "Heterozygosity"),
        xaxis = list(title = "Locus"),
        barmode = "group"
      )
  })

  output$fst_arlecore <- DT::renderDataTable(
    {
      req(arlequinFstMatrix())
      fst_data <- arlequinFstMatrix()[["data"]]
      fst_data[lower.tri(fst_data)] <- t(fst_data)[lower.tri(fst_data)]
      pop_names <- arlequinPopLabels()$Population
      rownames(fst_data) <- pop_names
      colnames(fst_data) <- pop_names
      fst_data
    },
    options = list(scrollX = TRUE)
  )

  output$fst_heatmap_plot_arlequin <- plotly::renderPlotly({
    req(arlequinFstMatrix())
    plot_heatmap_arlecore(
      arlequinFstMatrix()[["long"]],
      arlequinPopLabels(),
      legend_name = "FST"
    )
  })

  output$fst_pairwise_heatmap_plot <- plotly::renderPlotly({
    req(arlequinPairwise())
    data <- plot_pairwise_data_prep(
      arlequinPairwise(),
      arlequinPopLabels()
    )
    plot_pairwise_heatmap(
      data,
      arlequinPopLabels()
    )
  })

  output$fst_pairwise_heatmap_plot_overlap <- renderPlot({
    req(arlequinPairwise())
    data <- plot_pairwise_data_prep(
      arlequinPairwise(),
      arlequinPopLabels()
    )
    plot_pairwise_heatmap_overlap(
      data,
      arlequinPopLabels()
    )
  })

  output$coancestry_arlecore <- DT::renderDataTable(
    {
      req(arlequinCoancestry())
      coancestry_data <- arlequinCoancestry()[["data"]]
      coancestry_data[lower.tri(coancestry_data)] <- t(coancestry_data)[lower.tri(coancestry_data)]
      pop_names <- arlequinPopLabels()$Population
      rownames(coancestry_data) <- pop_names
      colnames(coancestry_data) <- pop_names
      coancestry_data
    },
    options = list(scrollX = TRUE)
  )

  output$coancestry_heatmap_plot <- plotly::renderPlotly({
    req(arlequinCoancestry())
    plot_heatmap_arlecore(
      arlequinCoancestry()[["long"]],
      arlequinPopLabels(),
      legend_name = "Coancestry Coefficient"
    )
  })


  output$ld_tables <- renderUI({
    req(arlequinLD())
    populations <- unique(arlequinLD()$Population)

    tagList(
      lapply(
        seq_along(populations),
        function(i) {
          population <- populations[i]

          tagList(
            h3(population),
            DT::DTOutput(paste0("ld_table_", i)),
            br()
          )
        }
      )
    )
  })

  observe({
    req(arlequinLD())
    data <- arlequinLD()
    populations <- unique(data$Population)

    lapply(seq_along(populations), function(i) {
      population <- populations[i]
      output_id <- paste0("ld_table_", i)

      local({
        pop <- population
        output[[output_id]] <- DT::renderDT({
          pop_data <- arlequinLD() %>%
            dplyr::filter(Population == pop) %>%
            dplyr::select(-Population)

          DT::datatable(pop_data,
            rownames = FALSE,
            options = list(
              pageLength = 10,
              scrollX = TRUE
            )
          )
        })
      })
    })
  })

  output$download_arlecore_results <- downloadHandler(
    filename = function() {
      paste0("Arlecore_results_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      wb <- openxlsx::createWorkbook()

      if (!is.null(arlequinHeterozygosity())) {
        hwe <- arlequinHeterozygosity()[["data"]]

        if (!is.null(hwe)) {
          openxlsx::addWorksheet(wb, "Heterozygosity")

          openxlsx::writeData(wb, "Heterozygosity", hwe)
        }
      }

      if (!is.null(arlequinFstMatrix())) {
        fst <- arlequinFstMatrix()[["data"]]

        if (!is.null(fst)) {
          # fst <- as.matrix(fst)
          fst[lower.tri(fst)] <- t(fst)[lower.tri(fst)]
          pop_names <- arlequinPopLabels()$Population

          rownames(fst) <- pop_names
          colnames(fst) <- pop_names

          fst_export <- cbind(
            Population = rownames(fst),
            as.data.frame(fst)
          )

          openxlsx::addWorksheet(wb, "FST Matrix")

          openxlsx::writeData(wb, "FST Matrix", fst_export)
        }
      }

      if (!is.null(arlequinCoancestry())) {
        coanc <- arlequinCoancestry()[["data"]]

        if (!is.null(coanc)) {
          # coanc <- as.matrix(coanc)
          coanc[lower.tri(coanc)] <- t(coanc)[lower.tri(coanc)]
          pop_names <- arlequinPopLabels()$Population

          rownames(coanc) <- pop_names
          colnames(coanc) <- pop_names

          coanc_export <- cbind(
            Population = rownames(coanc),
            as.data.frame(coanc)
          )

          openxlsx::addWorksheet(wb, "Coancestry Coefficient")

          openxlsx::writeData(wb, "Coancestry Coefficient", coanc_export)
        }
      }

      if (!is.null(arlequinPairwise())) {
        pairwise <- arlequinPairwise()
        used_sheet_names <- character(0)

        for (i in seq_along(pairwise)) {
          obj <- pairwise[[i]]

          if (!is.null(obj[["data"]])) {
            pair_data <- obj[["data"]]
            pair_data[lower.tri(pair_data)] <- t(pair_data)[lower.tri(pair_data)]
            pop_names <- arlequinPopLabels()$Population

            rownames(pair_data) <- pop_names
            colnames(pair_data) <- pop_names

            pair_data_export <- cbind(
              Population = rownames(pair_data),
              as.data.frame(pair_data)
            )

            sheet_name <- obj[["title"]]

            if (is.null(sheet_name) || sheet_name == "") {
              sheet_name <- paste0(
                "Pairwise ", i
              )
            }

            sheet_name <- gsub("[\\/:*?\\[\\]]", "-", sheet_name)
            sheet_name <- substr(sheet_name, 1, 31)
            original_name <- sheet_name
            count <- 2

            while (tolower(sheet_name) %in% tolower(used_sheet_names)) {
              suffix <- paste0(" (", count, ")")
              sheet_name <- paste0(
                substr(original_name, 1, 31 - nchar(suffix)),
                suffix
              )
              count <- count + 1
            }

            used_sheet_names <- c(used_sheet_names, sheet_name)
            openxlsx::addWorksheet(wb, sheet_name)

            openxlsx::writeData(wb, sheet_name, pair_data_export)
          }
        }
      }

      if (!is.null(arlequinPopDiversity())) {
        diversity_stats <- arlequinPopDiversity()

        populations <- unique(diversity_stats$Population)

        for (i in seq_along(populations)) {
          pop <- populations[i]

          pop_data <- diversity_stats %>%
            dplyr::filter(
              Population == pop
            ) %>%
            dplyr::select(-Population)

          sheet_name <- paste0("Diversity - ", pop)
          sheet_name <- substr(sheet_name, 1, 31)

          openxlsx::addWorksheet(wb, sheet_name)

          openxlsx::writeData(wb, sheet_name, pop_data)
        }
      }


      if (!is.null(arlequinLD())) {
        ld <- arlequinLD()

        populations <- unique(ld$Population)
        for (i in seq_along(populations)) {
          pop <- populations[i]

          pop_data <- ld %>%
            dplyr::filter(
              Population == pop
            ) %>%
            dplyr::select(-Population)

          sheet_name <- paste0("LD - ", pop)
          sheet_name <- substr(sheet_name, 1, 31)

          openxlsx::addWorksheet(wb, sheet_name)

          openxlsx::writeData(wb, sheet_name, pop_data)
        }
      }

      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    }
  )

  output$download_arp_file <- downloadHandler(
    filename = function() {
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
      paste0("arp_file", timestamp, ".arp")
    },
    content = function(file) {
      req(arlequinFile())
      file.copy(arlequinFile(), file)
    }
  )

  output$download_arp_result <- downloadHandler(
    filename = function() {
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
      paste0("arp_file", timestamp, ".arp")
    },
    content = function(file) {
      req(arlequinResult())
      file.copy(arlequinResult(), file)
    }
  )

  output$download_arlecore_results_UI <- renderUI({
    req(
      arlequinPopLabels(),
      arlequinHeterozygosity(),
      arlequinFstMatrix(),
      arlequinCoancestry(),
      arlequinPairwise(),
      arlequinPopDiversity(),
      arlequinResult()
    )
    tagList(
      downloadButton("download_arlecore_results", "Download Results (.xlsx)"),
      downloadButton("download_arp_file", "Download .arp file"),
      downloadButton("download_arp_result", "Download resulting .xml file")
    )
  })
}
