phylogeny_server <- function(input, output, session, rv) {
  # =========== PHYLOGENETIC TREE CONSTRUCTION ============#

  tree_plot <- reactiveVal(NULL)
  tree_path <- reactiveVal(NULL)
  tree_model <- reactiveVal(NULL)
  distanceEmpty <- reactiveVal(NULL)
  alignmentFile <- reactiveVal(NULL)
  distanceNAcount <- reactiveVal(NULL)

  observe({
    if (isTRUE(input$uploadMSA)) {
      toggleState("validateAlignment", !is.null(input$treeAlignmentType))
    } else {
      toggleState("validateAlignment", !is.null(input$msaFileforPhylogen))
    }
  })
  
  # separate validation of the distance matrix
  observeEvent(input$validateAlignment, {
    tryCatch({
      if(isTRUE(input$uploadMSA)) {
        alignment_file <- switch(
          input$treeAlignmentType,
          initial = rv$alignmentMSA,
          adjusted = rv$alignmentAdjusted,
          staggered = rv$alignmentStaggered
        )
      } else if (!is.null(input$msaFileforPhylogen)) {
        file_ext_ref <- tools::file_ext(
          input$msaFileforPhylogen$datapath
        )
        
        file_ext_ref <- tools::file_ext(input$msaFileforPhylogen$name)
        if (file_ext_ref == "msa") {
          alignment_file <- Biostrings::ReadDNAStringSet(input$msaFileforPhylogen$datapath)
        } else if (file_ext_ref %in% c("fasta", "fa", "fas", "msf", "aln")) {
          if (file_ext_ref == "aln") {
            file_ext_ref <- "clustal"
          }
          alignment_file <- seqinr::read.alignment(file = input$msaFileforPhylogen$datapath, format = file_ext_ref)
        } else {
          stop("Unsupported file format. ")
        }
      } else {
        stop("Please select/upload an alignment.")
      }
      
      alignmentFile(alignment_file)
      alignment_bins <- ape::as.DNAbin(alignment_file)
      distance_bins <- ape::dist.dna(alignment_bins, model = "K80")
      missing_count <- sum(is.na(distance_bins))
      distanceNAcount(missing_count)
      
      showNotification(paste(
        "Alignment validated.", missing_count, "empty values found."
      ), type = "message")
      
    }, error = function(e) {
      distanceNAcount(NULL)
      showNotification(paste(
        "Validation error", e$message), type = "error")
    })
  })
  
  output$phylogeneticOptionsUI <- renderUI({
    req(!is.null(distanceNAcount()))
    
    missing_count <- distanceNAcount()
    if (missing_count >= 1) {
      tree_choices <- c("Neighbor Joining (NJ)" = "NJ")
      
      missing_message <- h4(paste(missing_count, "missing distances found.",
                                  "Only NJ is available. Trim the alignment to use other methods."))
    } else {
      tree_choices <- c("Neighbor Joining (NJ)" = "NJ", 
                        "UPGMA" = "UPGMA", 
                        "Parsimony" = "Parsimony", 
                        "Maximum Likelihood" = "Maximum Likelihood")
      missing_message <- h4(paste(missing_count, "missing distances found."))
    }
    
    tagList(
      missing_message,
      selectInput("treeType", "Choose Method for Tree Construction", choices = tree_choices),
      conditionalPanel(
        condition = "input.treeType == 'NJ' || input.treeType == 'UPGMA'",
        selectInput("model", "Choose Substitution Model",
                    choices = c("N", "TS", "TV", "JC69", "K80", "F81", "K81", "F84", "BH87", "T92", "TN93", "GG95"),
                    selected = "K80"
        ),
        checkboxInput("pairwiseDeletion", "Perform pairwise deletion when calculating distance", value = FALSE)
      ),
      conditionalPanel(
        condition = "input.treeType == 'Maximum Likelihood'",
        textInput("bootstrapSamples", "Set number of bootstrap samples", value = 100)
      ),
      textInput("outgroup", "Outgroup (optional)", placeholder = "e.g. Sample1"),
      textInput("seed", "Set Seed Value", value = 123),
      actionButton("buildTree", "Build Tree", icon = icon("tree"))
    )
  })
  
  
  observeEvent(input$buildTree, {
    disable("buildTree")
    tree_type <- input$treeType
    outgroup <- input$outgroup
    bs <- input$bootstrapSamples
    model <- input$model
    aligned <- alignmentFile()
    directory <- tempdir()
    
    withProgress(message = "Building phylogenetic tree...", value = 0, {
      tryCatch(
        {
          if (!is.null(input$seed)) {
            seed <- input$seed
          } else {
            seed <- 123
          }

          if (tree_type == "NJ") {
            plot_obj <- build_nj_tree(aligned, outgroup = outgroup, model = model, seed = seed, pairwise = input$pairwiseDeletion)
            tree_plot(plot_obj)
            tree_path(NULL)
            tree_model(paste("NJ C", model, ")", sep = ""))
          } else if (tree_type == "UPGMA") {
            plot_obj <- build_upgma_tree(aligned, outgroup = outgroup, model = model, seed = seed, pairwise = input$pairwiseDeletion)
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
