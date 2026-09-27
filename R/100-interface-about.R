about_tab <- function() {
  tabItem(
    tabName = "About",
    fluidRow(
      tabBox(
        title = tagList(icon("book-open-reader"), "Background"),
        width = 12,
        # height = "250px",
        h4("The restful forensics toolkit is an output of the project", strong("Development and validation of an automated web-based tool for efficient genomic marker extraction to assist in genomic research.")),
        h4("This is based on the preliminary work on ancestry marker analysis at the DNA Analysis Laboratory, Natural Sciences Research Institute, University of the Philippines Diliman."),
        p(strong("Primary Developer:"), "Leda Celeste Samin (DAL-NSRI-UPD)"),
        p(strong("Project Leader:"), "Nelvie Fatima Jane Soliven (DAL-NSRI-UPD)"),
        p(strong("Contributors and Collaborators:")),
        p("Melvin Ambrocio Matias (Institute of Biology - UPD)"),
        p("Jazelyn Salvador (DAL-NSRI-UPD)"),
        p("Maria Corazon De Ungria (DAL-NSRI-UPD)")
      )
    ),
    fluidRow(
      tabBox(
        title = tagList(icon("wallet"), "Funding"),
        width = 6,
        height = "250px",
        h4("The project is funded by the Natural Sciences Research Institute at the University of the Philippines Diliman (Grant NSR-25-1-02)."),
        div(
          tags$img(
            src = "funding.png",
            width = "200px"
          ),
          style = "text-align: center;"
        )
      ),
      tabBox(
        title = tagList(icon("address-book"), "Contact"),
        width = 6,
        height = "250px",
        h4(strong("Project Leader: "), "nasoliven@up.edu.ph"),
        h4(strong("Institutional Contact: "), "dnalab.updiliman@up.edu.ph"),
        h4(strong("Location:"), "Miranda Hall, Natural Sciences Research Institute, University of the Philippines Diliman, Quezon City, Philippines")
      )
    )
  )
}
