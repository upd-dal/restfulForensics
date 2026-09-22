forensic_params_server <- function(input, output, session, rv) {
  # =============== FORENSIC PARAMETERS ===================#

  referenceData <- data.frame(
    Sample = c("Sample1", "Sample2", "Sample3", "Sample4", "..."),
    Population = c("POP1", "POP2", "POP3", "POP4", "..."),
    rs101 = c("A/A", "A/T", "A/A", "T/T", "..."),
    rs102 = c("G/G", "C/C", "G/C", "G/G", "..."),
    rs_n = c("...", "...", "...", "...", "...")
  )

  afSample <- data.frame(
    markers = c("rs101.A", "rs101.T", "rs102.C", "rs102.G", "..."),
    POP1 = c("0.18518", "0.81481", ".77777", "0.22222", "..."),
    POP2 = c("0.89285", "0.10714", "0.89285", "0.10714", "..."),
    POP3 = c("0.15789", "0.84210", "0.87894", "0.12105", "..."),
    POPn = c("...", "...", "...", "...", "...")
  )

  profileSample <- data.frame(
    markers = c("rs101", "rs102", "rs103", "rs104", "..."),
    profile = c("A/T", "G/C", "G/A", "T/T", "...")
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

  output$afSample_UI <- DT::renderDataTable(
    {
      req(afSample)
      afSample
    },
    options = list(
      scrollX = TRUE,
      pageLength = 5
    )
  )

  output$profileSample_UI <- DT::renderDataTable(
    {
      req(profileSample)
      profileSample
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

    fileUploaded <- load_csv_xlsx_files(input$iisnpsFile$datapath)
    cleaned_data <- clean_input_data(fileUploaded)
    genind_input <- convert_to_genind(cleaned_data, to_str = FALSE, popinfo = TRUE)
    af_table <- compute_af(genind_input)
    af_expected <- calc_expected_genotype_freq(af_table)
    
    gt_freqs <- calc_observed_genotype_freq(cleaned_data) # returns list of per population gt
    params_res <- calc_iisnps_params(gt_freqs, af_expected)
    print(params_res)
    print(class(params_res))
    genoFreq(gt_freqs) # list of df per population containing the observed freq
    forenParams(params_res) # list of df per population containing the forensic param metrices
    afTable(af_table)
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
  
  rmpResult <- reactiveVal(NULL)
  
  observeEvent(input$calcRMP, {
    disable("calcRMP")
    profile <- load_csv_xlsx_files(input$fileProfile$datapath)
    if (isFALSE(input$newPopDatabase)) {
      req(afTable())
      af_long <- af_long(afTable()) %>%
        dplyr::filter(population == input$rmp_population)
    } else {
      req(input$newPopDataFile)
      af_long <- load_csv_xlsx_files(input$newPopDataFile$datapath)
    }
    result <- calc_rmp(
      profile = profile,
      af_table = af_long,
      n_ref = input$totalPop,
      theta = thetaValue
    )
    rmpResult(result)
    enable("calcRMP")
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
