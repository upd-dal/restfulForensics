classification_tab <- function() {
  # CLASSIFICATION MODULE ==========================================================================
  tabItem(
    tabName = "Classification",
    fluidRow(
      box(
        fileInput("forPredFile", "Upload file", accept = c(".csv", "xlsx", "xlsm", "xlsb", "xls")),
        actionButton("runNaiveBayes", "Classify", icon = icon("align-justify")),
        uiOutput("downloadClassification_UI")
      ),
      tabBox(
        tabPanel(
          title = "Instructions",
          h4("Perform NaÏve Bayes classification"),
          p(strong("Input file/s:"), "CSV/Excel file containing training data or merged training and test data."),
          p(strong("Expected output file:"), "Excel (.xlsx) file"),
          hr(),
          p(
            "This performs a NaÏve Bayes classification and leave-one-out cross-validation
                                  of individuals with SNP data (ancestry- or phenotype- informative) using the ",
            tags$a("e1071",
              href = "https://cran.r-project.org/web/packages/e1071/index.html",
              target = "_blank"
            ),
            " and ",
            tags$a("caret",
              href = "https://cran.r-project.org/web/packages/caret/index.html",
              target = "_blank"
            ),
            " R packages."
          ),
          p(
            "This tool may be used for rapid assessment of marker sets and training datasets. For
                                  more comprehensive forensic DNA inference of an individual, it is recommended to use this tool
                                  in conjunction with",
            tags$a(actionLink("topcatab", "PCA,")),
            tags$a(actionLink("tostructuretab", "STRUCTURE")), "and",
            tags$a("SNIPPER.",
              href = "https://mathgene.usc.es/snipper/index.php",
              target = "_blank"
            ),
            " Multiple resources are available detailing the limiation of each classification method
                                   (Cheung et al., 2016; Barash et al., 2024)."
          )
        ),
        tabPanel(
          title = "Sample Input File",
          DT::dataTableOutput("classificationRef_UI")
        ),
        tabPanel(
          title = "Download Sample File",
          tags$ul(
            tags$a("Sample CSV file", href = "sample.csv", download = "sample.csv")
          )
        )
      )
    ),
    fluidRow(
      tabBox(
        title = "Results",
        width = 12,
        tabPanel(
          title = "Prediction Table",
          DT::dataTableOutput("predictionTableResult")
        ),
        tabPanel(
          title = "Statistics by Population",
          DT::dataTableOutput("statbyClassResult")
        ),
        tabPanel(
          title = "Overall Statistics",
          DT::dataTableOutput("overallStatResult")
        ),
        tabPanel(
           title = "Prediction by Individual",
           DT::dataTableOutput("predictionList")
        )
      )
    )
  )
}
