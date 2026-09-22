exploratory_analysis_server <- function(input, output, session, rv) {
  # ======================= PCA ===========================#
  GenindData <- reactiveVal(NULL)
  PCAResults <- reactiveVal(NULL)
  LabelColors <- reactiveVal(NULL)
  Populations <- reactiveVal(NULL)
  SelectedPops <- reactiveVal(NULL)

  output$examplePCA <- renderTable({
    data.frame(
      Sample = c("Sample1", "Sample2", "Sample3", "Sample4", "..."),
      Population = c("POP1", "POP2", "POP3", "POP4", "..."),
      rs101 = c("A/A", "A/T", "A/A", "T/T", "..."),
      rs102 = c("G/G", "C/C", "G/C", "G/G", "..."),
      rs_n = c("...", "...", "...", "...", "...")
    )
  })

  observe({
    hasDataFile <- !is.null(input$pcaFile)
    toggleState("runPCA", hasDataFile)
  })

  observe({
    toggleState("recalcPCA", !is.null(input$highlightPops))
  })

  output$selectedPopulation <- renderUI({
    req(Populations())

    shinyWidgets::prettyCheckboxGroup(
      inputId = "highlightPops",
      label = "Highlight Populations",
      choices = Populations(),
      inline = TRUE,
      selected = NULL,
      icon = icon("xmark")
    )
  })

  observeEvent(input$runPCA, {
    disable("runPCA")
    req(input$pcaFile)

    withProgress(message = "Running PCA...", {
      tryCatch(
        {
          df <- load_csv_xlsx_files(input$pcaFile$datapath)
          cleaned <- clean_input_data(df)

          val <- cleaned[1, 2]
          is_a_char <- stringr::str_count(val, "[A-Za-z]") == 2

          with_popinfo <- !isTRUE(is_a_char)

          fsnps_gen <- convert_to_genind(cleaned, to_str = FALSE, popinfo = with_popinfo)
          GenindData(fsnps_gen)

          label_file <- NULL

          if (!input$useDefaultColors) {
            req(input$pcaStyleFile)
            label_file <- input$pcaStyleFile$datapath
          }

          labels_colors <- get_labels(
            fsnps_gen = fsnps_gen,
            use_default = input$useDefaultColors,
            label_file = label_file,
            popinfo = with_popinfo
          )
          LabelColors(labels_colors)

          if (with_popinfo) {
            pops <- unique(as.character(adegenet::pop(fsnps_gen)))
            Populations(pops)
            SelectedPops(pops)
          } else {
            Populations(NULL)
            SelectedPops(NULL)
          }

          pca_results1 <- compute_pca(fsnps_gen, popinfo = with_popinfo)
          PCAResults(pca_results1)

          enable("runPCA")
        },
        error = function(e) {
          showNotification(paste("PCA Error:", e$message), type = "error")
          enable("runPCA")
        }
      )
    })
  })

  observeEvent(input$recalcPCA, {
    req(GenindData())
    req(input$highlightPops)
    disable("recalcPCA")

    withProgress(message = "Recalculating PCA...", {
      tryCatch(
        {
          filtered_pops <- subset_genind_pop(GenindData(), input$highlightPops)
          pca_results2 <- compute_pca(filtered_pops, popinfo = TRUE)
          PCAResults(pca_results2)

          Populations(unique(as.character(adegenet::pop(filtered_pops))))
        },
        error = function(e) {
          showNotification(paste("PCA recalculation error:", e$message),
            type = "error"
          )
        }
      )
    })
    enable("recalcPCA")
  })

  output$barPlot <- renderPlot({
    req(PCAResults())

    barplot(PCAResults()$percent,
      ylab = "Genetic variance explained by eigenvectors (%)", ylim = c(0, 25),
      names.arg = round(PCAResults()$percent, 1)
    )
  })

  output$pcaPlot <- plotly::renderPlotly({
    req(PCAResults(), LabelColors())

    p <- plot_pca(
      ind_coords = PCAResults()$ind_coords,
      centroid = PCAResults()$centroid,
      percent = PCAResults()$percent,
      labels_colors = LabelColors(),
      pc_x = input$pcX,
      pc_y = input$pcY,
      highlight_pop = input$highlightPops
    )
    plotly::ggplotly(p) %>% 
      plotly::layout(autosize = TRUE)
  })

  output$downloadbarPlot <- downloadHandler(
    filename = function() {
      paste0("bar_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      png(file, width = 800, height = 800, res = 300)
      graphics::barplot(
        PCAResults()$percent,
        ylab = "Genetic variance explained by eigenvectors (%)",
        ylim = c(0, 25),
        names.arg = round(PCAResults()$percent, 1)
      )
      dev.off()
    },
    contentType = "image/png"
  )

  output$downloadbarPlot_UI <- renderUI({
    req(PCAResults())
    downloadButton("downloadbarPlot", "Download Bar Plot")
  })

  output$downloadPCAPlot <- downloadHandler(
    filename = function() {
      paste0("pca_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      plot <- plot_pca(
        ind_coords = PCAResults()$ind_coords,
        centroid = PCAResults()$centroid,
        percent = PCAResults()$percent,
        labels_colors = LabelColors(),
        pc_x = input$pcX,
        pc_y = input$pcY,
        highlight_pop = input$highlightPops
      )

      ggsave(filename = file, plot = plot, width = 8, height = 8, dpi = 600)
    },
    contentType = "image/png"
  )

  output$downloadPCAPlot_UI <- renderUI({
    req(PCAResults())
    downloadButton("downloadPCAPlot", "Download PCA Plot")
  })
}
