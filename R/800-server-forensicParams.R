forensic_params_server <- function(input, output, session, rv) {
  # =============== FORENSIC PARAMETERS ===================#

  referenceData <- data.frame(
    Sample = c("Sample1", "Sample2", "Sample3", "Sample4", "..."),
    Population = c("POP1", "POP2", "POP3", "POP4", "..."),
    rs101 = c("A/A", "A/T", "A/A", "T/T", "..."),
    rs102 = c("G/G", "C/C", "G/C", "G/G", "..."),
    rs_n = c("...", "...", "...", "...", "...")
  )
  
  output$referenceData_UI <- DT::renderDataTable(
    {
      req(referenceData)
      referenceData
    },
    options = list(
      scrollX = TRUE,
      pageLength = 5
    )
  )

  forenParams <- reactiveVal(NULL)
  genoFreq <- reactiveVal(NULL)
  afTable <- reactiveVal(NULL)
  
  observe({
    shinyjs::toggleState("calcIISNPs", !is.null(input$iisnpsFile))
  })

  observeEvent(input$calcIISNPs, {
    shinyjs::disable("calcIISNPs")
    req(input$iisnpsFile)

    withProgress(message = "Ongoing: ", {
      
    incProgress(0.2, detail = "Loading input file...")
    fileUploaded <- load_csv_xlsx_files(input$iisnpsFile$datapath)
    cleaned_data <- clean_input_data(fileUploaded)
    tryCatch(
      {
    genind_input <- convert_to_genind(cleaned_data, to_str = FALSE, popinfo = TRUE)
    
    incProgress(0.4, detail = "Calculating frequencies...")
    af_table <- compute_af(genind_input)
    af_expected <- calc_expected_genotype_freq(af_table)
    
    gt_freqs <- calc_observed_genotype_freq(cleaned_data) # returns list of per population gt
    
    incProgress(0.6, detail = "Calculating parameters...")
    params_res <- calc_iisnps_params(gt_freqs, af_expected)
    genoFreq(gt_freqs)
    forenParams(params_res)
    afTable(af_table)
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
    })
    shinyjs::enable("calcIISNPs")
  })

  observeEvent(forenParams(), {
    pops <- names(forenParams())
    req(length(pops) > 0)

    updateSelectInput(
      session,
      "selected_pop",
      choices = pops,
      selected = pops[1]
    )
  })

  observeEvent(genoFreq(), {
    pops <- names(genoFreq())
    req(length(pops) > 0)

    updateSelectInput(
      session,
      "selected_pop_gt",
      choices = pops,
      selected = pops[1]
    )
  })

  output$genotypeFreqs_UI <- DT::renderDataTable({
    req(genoFreq())
    df <- genoFreq()[[input$selected_pop_gt]]
    DT::datatable(df, rownames = FALSE, selection = "multiple",
                  options = list(
                    pageLength = 10,
                    scrollX = TRUE
                  ))
  })

  output$popTable <- DT::renderDataTable({
    req(forenParams())
    df <- forenParams()[[input$selected_pop]]
    DT::datatable(df, rownames = FALSE, selection = "multiple",
                  options = list(
                    pageLength = 10,
                    scrollX = TRUE
                  ))
  })

  output$downloadMetrics <- downloadHandler(
    filename = function() {
      paste0("forensic_metrics_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      req(forenParams(), genoFreq())
      sheets <- list()
      fp_sheets <- forenParams()
      
      names(fp_sheets) <- paste0("FP_", substr(gsub(
        "[^A-Za-z0-9]", "_",
        names(fp_sheets)
      ), 1, 25))
      sheets <- c(sheets, fp_sheets)

      gt_sheets <- genoFreq()
      names(gt_sheets) <- paste0("GT_", substr(gsub(
        "[^A-Za-z0-9]", "_",
        names(gt_sheets)
      ), 1, 25))
      sheets <- c(sheets, gt_sheets)
      writexl::write_xlsx(sheets, path = file)
    }
  )

  output$downloadMetrics_UI <- renderUI({
    req(forenParams(), genoFreq())
    downloadButton("downloadMetrics", "Download Forensic Parameters")
  })


  #================= Matching Profile
  forensicParamInput <- data.frame(
    markers = c("rs101.A", "rs101.T", "rs102.C", "rs102.G", "..."),
    pop1 = c("0.185185185	", "0.814814815", "0.777777778", "0.222222222", "..."),
    pop2 = c("0.89285714", "0.10714286", "0.89285714", "0.10714286", "...")
  )
  
  output$forensicParamInput_UI <- DT::renderDataTable(
    {
      req(forensicParamInput)
      forensicParamInput
    },
    options = list(
      scrollX = TRUE,
      pageLength = 5
    )
  )
  
  sampleProfile <- data.frame(
    markers = c("rs101", "rs102", "rs103", "rs104", "..."),
    genotype = c("A/A	", "G/T", "A/T", "C/G", "...")
  )
  
  output$sampleProfile_UI <- DT::renderDataTable(
    {
      req(sampleProfile)
      sampleProfile
    },
    options = list(
      scrollX = TRUE,
      pageLength = 5
    )
  )
  
  rmpResult <- reactiveVal(NULL)
  
  observeEvent(input$calcRMP, {
    disable("calcRMP")
    
    withProgress(message = "Ongoing: ", {
      
    incProgress(0.2, detail = "Loading input file...")
    profile <- load_csv_xlsx_files(input$fileProfile$datapath)
    profile <- clean_input_data(profile)
    
    tryCatch({
      
      incProgress(0.4, detail = "Loading reference database...")
      if (isFALSE(input$newPopDatabase)) {
      req(afTable())
      af_long <- af_to_long(afTable()) %>%
        dplyr::filter(population == input$rmp_population)
    } else {
      req(input$newPopDataFile)
      af_long <- load_csv_xlsx_files(input$newPopDataFile$datapath)
      af_long <- af_to_long(af_long)
    }
    
    af_long <- prep_af_table(af_long)
    incProgress(0.6, detail = "Calculating RMP..")
    result <- calc_rmp(
      profile = profile,
      af_table = af_long,
      n_ref = input$totalPop,
      theta = input$thetaValue
    )
    rmpResult(result)
    showNotification("Calculation complete!", type = "message")
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
    
    enable("calcRMP")
    })
    
    
    })
  
  output$rmp_summary <- renderTable({
    req(rmpResult())
    result <- rmpResult()
    data.frame(
      Metric = c(
        "Random Match Probability (RMP)",
        "Likelihood Ratio"
      ),
      Value = c(
        format(result$rmp, scientific = TRUE, digits = 6),
        format(result$lr, scientific = TRUE, digits = 6)
      )
    )
  })
  
  output$locus_information <- renderTable({
    req(rmpResult())
    result <- rmpResult()
    
    result$locus_results %>%
      dplyr::select(
        marker_std,
        genotype_std,
        allele1,
        allele2,
        p,
        q,
        gen_prob
      )
  }, digits = 6)
  
  output$profileData <- renderTable({
    req(rmpResult())
    result <- rmpResult()
    result$locus_results %>%
      dplyr::select(marker, genotype, marker_std, genotype_std)
  })
  
  output$rmp_population_UI <- renderUI({
    if (isFALSE(input$newPopDatabase)) {
      req(afTable())
      af_long <- af_to_long(afTable())
      selectInput("rmp_population",
                  "Select Reference Population",
                  choices = unique(af_long$population))
    } else { NULL }
  })

  
  }
