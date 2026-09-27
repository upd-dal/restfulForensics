phylogeny_server <- function(input, output, session, rv) {
  # =========== PHYLOGENETIC TREE CONSTRUCTION ============#

  tree_plot <- reactiveVal(NULL)
  tree_path <- reactiveVal(NULL)
  tree_model <- reactiveVal(NULL)

  observeEvent(input$buildTree, {
    disable("buildTree")

    if (isTRUE(input$uploadMSA)) {
      alignment <- switch(input$treeAlignmentType,
        initial = rv$alignmentMSA,
        adjusted = rv$alignmentAdjusted,
        staggered = rv$alignmentStaggered
      )
      alignment_file <- alignment
    } else if (!is.null(input$msaFileforPhylogen)) {
      file_ext_ref <- tools::file_ext(input$msaFileforPhylogen$name)
      if (file_ext_ref == "msa") {
        alignment_file <- Biostrings::ReadDNAStringSet(input$msaFileforPhylogen$datapath)
      } else if (file_ext_ref %in% c("fasta", "msf", "aln")) {
        if (file_ext_ref == "aln") {
          file_ext_ref <- "clustal"
        }
        alignment_file <- seqinr::read.alignment(file = input$msaFileforPhylogen$datapath, format = file_ext_ref)
      }
    }

    withProgress(message = "Building phylogenetic tree...", value = 0, {
      tryCatch(
        {
          tree_type <- input$treeType
          outgroup <- input$outgroup
          bs <- input$bootstrapSamples
          model <- input$model
          aligned <- alignment_file
          directory <- tempdir()

          if (!is.null(input$seed)) {
            seed <- input$seed
          } else {
            seed <- "123"
          }

          if (tree_type == "NJ") {
            plot_obj <- build_nj_tree(aligned, outgroup = outgroup, model = model, seed = seed)
            tree_plot(plot_obj)
            tree_path(NULL)
            tree_model(paste("NJ C", model, ")", sep = ""))
          } else if (tree_type == "UPGMA") {
            plot_obj <- build_upgma_tree(aligned, outgroup = outgroup, model = model, seed = seed)
            tree_plot(plot_obj)
            tree_path(NULL)
            tree_model(paste("UPGMA (", model, ")", sep = ""))
          } else if (tree_type == "Parsimony") {
            path <- build_max_parsimony(aligned, outgroup = outgroup, directory = directory, seed = seed)
            tree_plot(NULL)
            tree_path(path)
            tree_model("Parsimony")
          } else if (tree_type == "Maximum Likelihood") {
            results <- build_ml_tree(aligned, outgroup = outgroup, directory = directory, seed = seed, bs = bs)
            tree_plot(NULL)
            tree_path(results$filename)
            tree_model(results$best_model)
          }
          showNotification("Tree construction complete.", type = "message")
        },
        error = function(e) {
          showNotification(paste("Error during tree construction:", e$message), type = "error", duration = 20)
        }
      )
      enable("buildTree")
    })
  })

  output$treeImage <- renderUI({
    if (!is.null(tree_plot())) {
      plotOutput("treePlot")
    } else if (!is.null(tree_path())) {
      imageOutput("treePNG")
    }
  })

  output$treePlot <- renderPlot({
    req(tree_plot())
    tree_plot()
  })

  output$treePNG <- renderImage(
    {
      req(tree_path())
      list(src = tree_path(), contentType = "image/png", width = 800, height = 600)
    },
    deleteFile = FALSE
  )

  output$downloadTree <- downloadHandler(
    filename = function() {
      "tree.png"
    },
    content = function(file) {
      if (!is.null(tree_plot())) {
        ggsave(file, plot = tree_plot(), width = 8, height = 6, dpi = 600)
      } else if (!is.null(tree_path())) {
        file.copy(tree_path(), file)
      }
    }
  )

  output$downloadTree_UI <- renderUI({
    if (!is.null(tree_plot()) || !is.null(req(tree_path()))) {
      downloadButton("downloadTree", "Download Phylogenetic Tree")
    }
  })
}
