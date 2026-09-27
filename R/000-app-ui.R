app_ui <- function(request = NULL) {
  dashboardPage(
    restful_header(),
    restful_sidebar(),
    dashboardBody(
      useShinyjs(),
      tags$head(
        includeCSS("www/custom.css")
      ),
      shinybusy::add_busy_spinner(spin = "fading-circle", position = "full-page"),
      tags$style(HTML("
                      .shinybusy-overlay {
                      background-color: rgba(0,0,0,0.25) !important;
                      }
                      ")),
      tabItems(
        dashboard_tab(),
        file_conversion_tab(),
        snp_extraction_tab(),
        filtering_tab(),
        exploratory_tab(),
        popstats_tab(),
        # structure_tab(), to separate into two
        structure_runs(),
        plot_structure_runs(),
        forensic_params_tab(),
        classification_tab(),
        msa_tab(),
        phylogeny_tab(),
        barcoding_tab(),
        references_tab(),
        about_tab()
      ) # end of tabitems
    ) # end of dashboard body
  ) # end of page
}
