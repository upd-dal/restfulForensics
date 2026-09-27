plot_structure_runs <- function() {
  tabItem(
    tabName = "PlotStructure",
    fluidRow(
      box(
        checkboxInput("plotSR", "Use STRUCTURE results", value = TRUE),
        conditionalPanel(
          condition = "input.plotSR == false",
          fileInput("structureFilesZipped", "Upload zipped file of STRUCTURE results", accept = c("zip", "tar"))
        ),
        numericInput("kChoice", "Target K Value", value = 2, min = 1),
        helpText("To see other K values, run CLUMPP again."),
        selectInput("greedyAlgo", "Choose Algorithm",
          choices = c("Full Search" = "full.search", "Greedy" = "greedy", "Large K" = "large.k"),
          selected = "greedy"
        ),
        radioButtons("simStat", "Pairwise matrix similarity statistics",
          choices = c("g" = "g", "g.prime" = "g.prime")
        ),
        conditionalPanel(
          condition = "input.greedyAlgo == 'greedy' || input.greedyAlgo == 'large.k'",
          radioButtons("greedyOption", "Input order of the runs to be tested",
            choices = c("All" = "all", "ran order" = "ran.order"),
            selected = "ran.order"
          ),
          numericInput("repeatsClumpp", "Number of input orders of runs for testing", value = 100, min = 0)
        ),
        numericInput("orderRun", "Permute clusters according to cluster order?", value = 0),
        actionButton("plotStructureResults", "Run CLUMPP", icon = icon("play"))
      ),
      tabBox(
        tabPanel(
          "Instructions",
          h4("Run CLUMPP and plot STRUCTURE results"),
          p(strong("Input file:"), "Zipped structure files or structure result from Run STRUCTURE v2.3.4 tab."),
          p(strong("Expected output file:"), "PNG plot"),
          hr(),
          p(
            "Runs ",
            tags$a("CLUMPP",
              href = "https://rosenberglab.stanford.edu/clumpp.html",
              target = "_blank"
            ),
            "(Jakobsson & Rosenberg, 2007) to visualize STRUCTURE results using the ",
            tags$a("strataG",
              href = "https://github.com/EricArcher/strataG/tree/master",
              target = "_blank"
            ), "R package."
          )
        ),
        tabPanel(
          "Download Sample File",
          tags$ul(
            tags$a("Sample STRUCTURE results", href = "structure_results.zip", download = "structure_results.zip")
          )
        )
      )
    ),
    fluidRow(
      tabBox(
        title = "Results",
        width = 12,
        tabPanel(
          "CLUMPP results",
          uiOutput("downloadStructure_UI"),
          plotOutput("structure_result_plots")
        )
      )
    )
  )
}
