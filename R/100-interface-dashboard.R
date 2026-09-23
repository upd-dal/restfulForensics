dashboard_tab <- function() {
  tabItem(
    tabName = "dashboard",
    fluidRow(
      tabBox(
        title = "Introduction",
        side = "right",
        width = 12,
        h4("restful forensics (Reproducible and Efficient Sequence Toolkit) is a toolkit dedicated for
            the forensic analysis of single nucleotide polymorphisms (SNPs) and DNA barcodes. It compiles
                          reference population datasets extracted from publicly available databases
                          for direct evaluation of forensic marker panels. It allows the user to
                          analyze their dataset with compiled reference datasets, perform exploratory
                          data analysis (i.e. principal component analysis), and calculate population
                          summary statistics and forensic parameters. This version also contains modules
                          for population structure analysis (i.e. STRUCTURE), forensic DNA
                          inference/classification using ancestry or phenotype-informative SNPs, and DNA
                          barcoding (i.e. multiple sequencing alignment).")
      ),
      tabBox(
        title = "Citation",
        side = "right",
        width = 12,
        h4("Cite the application:"),
        h4("DNA Analysis Laboratory. (2025). restfulForensics (Version 1.0) [Computer software]. GitHub. github.com"),
        br(),
        h4("Manuscript pending.")
      ),
    ),
    fluidRow(
      tabBox(
        title = "Overview of Features",
        width = 12,
        tabPanel(
          "Workflow",
          # Zoom features adapted from: https://forum.posit.co/t/zoom-in-zoom-out-in-r-shiny-while-working-with-images/183567
          div(style = "display:flex; justify-content: space-evenly; margin-bottom:10px;"),
          actionButton("smaller_workplan", "-"),
          actionButton("bigger_workplan", "+"),
          div(
            style = "overflow:auto; text-align:center;",
            uiOutput("workplanImg")
          )
        ),
        tabPanel(
          title = "File Conversion Options",
          div(style = "display:flex; justify-content: space-evenly; margin-bottom:10px;"),
          actionButton("smaller_fc", "-"),
          actionButton("bigger_fc", "+"),
          div(
            style = "overflow:auto; text-align:center;",
            uiOutput("fileConvTable")
          )
        )
      )
    )
  )
}
