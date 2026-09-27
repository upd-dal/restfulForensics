app_server <- function(input, output, session) {
  rv <- reactiveValues(
    structureRes = NULL,
    alignmentMSA = NULL,
    alignmentAdjusted = NULL,
    alignmentStaggered = NULL
  )

  navigation_server(input, output, session)

  dashboard_server(input, output, session, rv)

  file_conversion_server(input, output, session, rv)

  snp_extraction_server(input, output, session, rv)

  filtering_server(input, output, session, rv)

  msa_server(input, output, session, rv)

  phylogeny_server(input, output, session, rv)

  barcoding_server(input, output, session, rv)

  pop_stats_server(input, output, session, rv)

  forensic_params_server(input, output, session, rv)

  exploratory_analysis_server(input, output, session, rv)

  run_structure_analysis(input, output, session, rv)

  plot_structure_server(input, output, session, rv)

  classification_server(input, output, session, rv)
}
