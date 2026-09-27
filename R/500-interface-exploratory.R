exploratory_tab <- function() {
  # PRINCIPAL COMPONENT ANALYSIS MODULE =============================================================
  tabItem(
    tabName = "PCAtab",
    fluidRow(
      box(
        fileInput("pcaFile", "Upload Genotype File", accept = c(".csv", ".txt", "xlsx", "xlsm", "xlsb", "xls")),
        checkboxInput("useDefaultColors", "Use Default Colors and Labels", TRUE),
        conditionalPanel(
          condition = "!input.useDefaultColors",
          fileInput("pcaStyleFile", "Customize population colors and shapes.", accept = c(".csv", "xlsx", "xlsm", "xlsb", "xls", ".txt")),
          helpText("Columns should contain: [1] Population name, [2] Color (name or hex code), [3] Shapes"),
          p("The order of the colors would match the order of PCA labels")
        ),
        br(),
        numericInput("pcX", "PC Axis X", value = 1, min = 1),
        numericInput("pcY", "PC Axis Y", value = 2, min = 1),
        uiOutput("selectedPopulation"),
        actionButton("runPCA", "Run PCA Analysis", icon = icon("play")),
        actionButton("recalcPCA", "Recalculate PCA for Selected Populations", icon = icon("filter"))
      ),
      tabBox(
        tabPanel(
          title = "Instructions",
          h4("Run principal component analysis using the ade4 (Dray and Dufour, 2007) package in R"),
          p(strong("Input file:"), "CSV or Excel file and color labels (optional)"),
          p(strong("Optional additional input file/s:"),
            "If using custom visualizations, upload a file (CSV/Excel) with columns containing 
            [1] Unique population name/s that matches the input file,
            [2] Color for a given population (name or hex code), and
            [3] Desired point", tags$a("shapes",
                href = "https://ggplot2.tidyverse.org/reference/scale_shape.html",
                target = "_blank"
              ), "for a given population"
            ),
          p(strong("Expected output file:"), "PNG plots"),
          br(),
          p(strong("Note (if analyzing many populations):")),
          p("If the input contains many populations, not all legends will be immediately visible under the Plots tab.
          Scroll through the legends to see other population labels. The 'Download PCA Plot' will produce a PNG file
          with labels positioned at the center of the population cluster. If using the 'Download plot as png' via the camera button,
          note that only a fixed number of populations will be included in the PNG file.")
        ),
        tabPanel(
          title = "Sample Input Format/s",
          h4("Example Input Format"),
          tableOutput("examplePCA")
        ),
        tabPanel(
          title = "Download sample files",
          tags$ul(
            tags$a("Sample file", href = "sample.csv", download = "sample.csv")
          )
        )
      )
    ),
    fluidRow(
      tabBox(
        title = "PCA Results",
        width = 12,
        tabPanel(
          title = "Plots",
          uiOutput("downloadPCAPlot_UI"),
          div(style = "height: 70vh; overflow-y: auto;", 
              plotly::plotlyOutput("pcaPlot", height = "100%") 
          )
        ),
        tabPanel(
          title = "Bar Plot",
          plotOutput("barPlot"),
          uiOutput("downloadbarPlot_UI")
        )
      )
    )
  )
}
