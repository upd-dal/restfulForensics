classification_tab <- function() {
  # CLASSIFICATION MODULE ==========================================================================
  tabItem(
    tabName = "Classification",
    tabsetPanel(
      tabPanel(
        title = "Leave-one-out",
        profvis::profvis_ui("profiler_loo"),
        fluidRow(
          box(
            fileInput("forPredFile", "Upload file", accept = c(".csv", "xlsx", "xlsm", "xlsb", "xls")),
            actionButton("runNaiveBayes", "Classify", icon = icon("align-justify")),
            uiOutput("downloadClassification_UI")
          ),
          tabBox(
            tabPanel(
              title = "Instructions",
              h4("Perform NaÏve Bayes classification using the leave-one-out method"),
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
      ),
      tabPanel(
        title = "Cross-validation",
        profvis::profvis_ui("profiler_cval"),
        fluidRow(
          box(
            fileInput("crossValFile", "Upload File", accept = c("xlsx", "xlsm", "xlsb", "xls", ".csv")),
            checkboxInput("optimizeVals", "Optimize "),
            numericInput("kValue", "Set K value for K-fold cross validation", value = 10),
            actionButton("performCrossval", "Perform cross-validation", icon = icon("shuffle"))
          ),
          tabBox(
            tabPanel(
              title = "Instructions",
              h4("Perform NaÏve Bayes classification using k-fold cross validation."),
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
              DT::dataTableOutput("classificationRefCV_UI")
            ),
            tabPanel(
              title = "Download Sample File",
              tags$ul(
                tags$a("Sample CSV file", href = "sample.csv", download = "sample.csv")
              )
            )
          )
        ), # end of first fluid row
        fluidRow(
          tabBox(
            title = "Results",
            width = 12,
            tabPanel(
              title = "Model Information",
              DT::dataTableOutput("modelInfoCV")
            ),
            tabPanel(
              title = "Prediction Table",
              uiOutput("downloadPredictionTableCV_UI"),
              DT::dataTableOutput("predictionTableResultCV")
            ),
            tabPanel(
              title = "Statistics by Population",
              uiOutput("downloadstatbyClassResultCV_UI"),
              DT::dataTableOutput("statbyClassResultCV")
            ),
            tabPanel(
              title = "Overall Statistics",
              uiOutput("downloadoverallStatResultCV_UI"),
              DT::dataTableOutput("overallStatResultCV")
            ),
            tabPanel(
              title = "Prediction by Individual",
              uiOutput("downloadpredictionListCV_UI"),
              DT::dataTableOutput("predictionListCV")
            )
          )
        )
      )
    )
  )
}
