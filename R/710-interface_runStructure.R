structure_runs <- function() {
  # POPULATION STRUCTURE MODULE ===========================================================================
  tabItem(
    tabName = "StructureRun",
    fluidRow(
      box(
        fileInput("structureFile", "Upload File", accept = c(".csv", "xlsx", "xlsm", "xlsb", "xls")),
        helpText("Input file should be similar to the output of the 'Convert to CSV' tab under 'File Conversion'"),
        numericInput("kMin", "Min K", value = 2, min = 1),
        numericInput("kMax", "Max K", value = 5, min = 1),
        numericInput("numKRep", "Replicates per K", value = 10, min = 1),
        numericInput("burnin", "Burn-in Period", value = 100000),
        numericInput("numreps", "MCMC Reps After Burn-in", value = 100000),
        br(),
        checkboxInput("advancedStructure", "See Advanced Parameters", value = FALSE),
        conditionalPanel(
          condition = "input.advancedStructure == true",
          checkboxInput("freqScore", "Use 'Correlated Frequencies' Model", value = FALSE),
          #checkboxInput("popAlpha", "Infer separate alpha per population", value = TRUE),
          checkboxInput("noadmix", "Use 'No Admixture' Model", value = TRUE),
          conditionalPanel(
            condition = "input.noadmix == false",
            checkboxInput("inferAlpha", "Infer the value of model parameters", value = TRUE),
            conditionalPanel(
              condition = "input.inferAlpha == false",
              numericInput("alphaValStructure", "Alpha Value", value = 0.05, max = 1)
            )
          )
        ),
        actionButton("runStructure", "Run STRUCTURE", icon = icon("play")),
        uiOutput("downloadButtons")
      ),
      tabBox(
        tabPanel(
          title = "Instructions",
          h4("Run population structure analysis using STRUCTURE v2.3.4"),
          p(strong("Input file:"), "CSV or Excel file"),
          p(strong("Expected output file:"), "Zipped qmatrices and individual files"),
          hr(),
          p(
            "This runs the basic Windows implementation of ",
            tags$a("STRUCTURE v2.3.4",
              href = "https://web.stanford.edu/group/pritchardlab/structure_software/release_versions/v2.3.4/structure_doc.pdf",
              target = "_blank"
            ),
            "without a front-end using the",
            tags$a("strataG",
              href = "https://github.com/EricArcher/strataG/tree/master",
              target = "_blank"
            ), "R package. This generates STRUCTURE input files and qmatrices files compatible for
                                   other programs such as pong (Behr et al., 2016) or CLUMPP (Jakobsson & Rosenberg, 2007)."
          )
        ),
        tabPanel(
          title = "Sample Input File/s",
          DT::dataTableOutput("examplePop_STR2UI")
        ),
        tabPanel(
          title = "Download Sample File",
          tags$ul(
            tags$a("Sample CSV file", href = "sample.csv", download = "sample.csv")
          )
        )
      )
    )
  )
}
