msa_server <- function(input, output, session, rv) {
  # ============ MULTIPLE SEQUENCE ALIGNMENT ==============#

  fasta_data <- reactiveVal(NULL)
  alignment_scores <- reactiveVal(NULL)
  directory <- tempdir()

  observeEvent(input$runMSA, {
    req(input$fastaFile)
    
    tryCatch({
      fasta <- read_fasta(input$fastaFile$datapath, directory)
    fasta_data(fasta)

    aligned <- calc_msa(fasta,
      algorithm = input$substitutionMatrix
    )
    showNotification("Alignment complete!", type = "message")
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    }
      )

    rv$alignmentMSA <- aligned$alignment
    rv$alignmentAdjusted <- aligned$adjusted
    rv$alignmentStaggered <- aligned$staggered
    alignment_scores(aligned$scores)
  })

  output$msaView <- msaR::renderMsaR({
    req(rv$alignmentMSA)
    msaR::msaR(rv$alignmentMSA)
  })

  output$adjustedAlignmentText <- renderPrint({
    req(rv$alignmentAdjusted)
    rv$alignmentAdjusted
  })

  output$staggeredAlignmentText <- renderPrint({
    req(rv$alignmentStaggered)
    rv$alignmentStaggered
  })

  output$alignmentScoresPreview <- renderPrint({
    req(alignment_scores())
    alignment_scores()
  })

  output$downloadAlignedFASTA <- downloadHandler(
    filename = function() {
      "aligned_sequences.fa"
    },
    content = function(file) {
      alignment <- switch(input$msaDownloadType,
        initial = rv$alignmentMSA,
        adjusted = rv$alignmentAdjusted,
        staggered = rv$alignmentStaggered
      )
      if (inherits(alignment, "MsaDNAMultipleAlignment")) {
        aln <- msa::msaConvert(alignment, type = "seqinr::alignment")
        seqs <- aln$seq
        names(seqs) <- aln$nam
        seqinr::write.fasta(
          sequences = as.list(seqs),
          names = names(seqs),
          file.out = file
        )
      } else if (inherits(alignment, "DNAStringSet")) {
        Biostrings::writeXStringSet(alignment, filepath = file, format = "fasta")
      } else {
        stop("Unsupported alignment type")
      }
    }
  )

  output$downloadAlignmentScores <- downloadHandler(
    filename = function() {
      "alignment_scores.txt"
    },
    content = function(file) {
      writeLines(utils::capture.output(print(alignment_scores())), file)
    }
  )

  output$downloadAlignmentPDF <- downloadHandler(
    filename = function() {
      "aligned_seqs.pdf"
    },
    contentType = "application/pdf",
    content = function(file) {
      alignment <- switch(input$msaDownloadType,
        initial = rv$alignmentMSA,
        adjusted = rv$alignmentAdjusted,
        staggered = rv$alignmentStaggered
      )

      tmp_dir <- tempfile("msa_")
      dir.create(tmp_dir)

      old_wd <- getwd()
      setwd(tmp_dir)

      on.exit(setwd(old_wd), add = TRUE)

      msa_file <- "alignment.pdf"

      msa::msaPrettyPrint(
        alignment,
        file = msa_file,
        output = "pdf",
        showNames = "left",
        showLogo = "none",
        askForOverwrite = FALSE
      )
      file.copy(msa_file, file, overwrite = TRUE)
      unlink(tmp_dir, recursive = TRUE)
    }
  )

  output$downloadAlignedFASTA_UI <- renderUI({
    downloadButton("downloadAlignedFASTA", "Download Aligned Sequences")
  })

  output$downloadAlignmentScores_UI <- renderUI({
    req(alignment_scores())
    downloadButton("downloadAlignmentScores", "Download Alignment Scores")
  })
}
