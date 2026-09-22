popstats_tab <- function() {
  # POPULATION STATISTICS MODULE=====================================================================
  tabItem(
    tabName = "PopStatistics",
    tabPanel(
      title = "Population Statistics",
      tabsetPanel(
        tabPanel(
          title = "R-based Calculations",
          fluidRow(
            box(
              fileInput("popStatsFile", "Upload CSV or XLSX Dataset", accept = c(".xlsx", ".csv")),
              numericInput("markovHWE", "Set Markov Chains for HWE calculations", value = 1000, min = 1000),
              selectInput("correctionModel", "Select Correction Model", choices = c("Bonferroni" = "Bonferroni", "FDR" = "FDR")),
              numericInput("alphaValue", "Set Alpha Value", value = 0.05, min = 0.00, max = 1),
              actionButton("runPopStats", "Analyze", icon = icon("magnifying-glass-chart")),
              uiOutput("downloadStatsXLSX_UI")
            ),
            tabBox(
              tabPanel(
                title = "Instructions",
                div(
                  style = "height: 40vh; overflow-y:scroll;",
                  uiOutput("popstatRef")
                )
              ),
              tabPanel(
                title = "Sample Input Format/s",
                h4("Example: Population File Format"),
                DT::dataTableOutput("examplePop_UI")
              ),
              tabPanel(
                title = "Download Sample Files",
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
                title = "1. Private Alleles",
                h4("Private Alleles Summary"),
                DT::dataTableOutput("privateAlleleTable")
              ),
              tabPanel(
                title = "2. Mean Allelic Richness",
                h4("Mean Allelic Richness per site"),
                DT::dataTableOutput("meanallelic")
              ),
              tabPanel(
                title = "3. Heterozygosity",
                h4("Observed vs Expected Heterozygosity"),
                DT::dataTableOutput("heterozygosity_table"),
                br(),
                h4("Heterozygosity Plot"),
                imageOutput("heterozygosity_plot"),
                uiOutput("downloadHeterozygosityPlot_UI")
              ),
              tabPanel(
                title = "4. Inbreeding Coefficients",
                h4("Inbreeding Coefficient by Population"),
                DT::dataTableOutput("inbreeding_table")
              ),
              tabPanel(
                title = "5. Allele Frequencies",
                h4("Allele Frequency Table"),
                DT::dataTableOutput("allele_freq_table")
              ),
              tabPanel(
                title = "6. Hardy-Weinberg Equilibrium",
                h4("HWE P-value Summary"),
                DT::dataTableOutput("hwe_summary_text"),
                h4("Population-wise HWE Chi-Square Table"),
                DT::dataTableOutput("hwe_chisq_table"),
                h4("Loci out of HWE"),
                DT::dataTableOutput("hwe_loci"),
                h4("Populations out of HWE"),
                DT::dataTableOutput("hwe_pop")
              ),
              tabPanel(
                title = "7. Fst Values",
                h4("Pairwise Fst Matrix"),
                uiOutput("fstMatrixUI"),
                h4("Tidy Pairwise Fst Data"),
                DT::dataTableOutput("fstDfTable"),
                br(),
                h4("Fst Heatmap"),
                uiOutput("downloadFstHeatmap_UI"),
                imageOutput("fst_heatmap_plot", width = "100%")
              )
            )
          )
        ),
        tabPanel(
          title = "Arlecore",
          fluidRow(
            box(
              fileInput("fileForArlecore", "Input file (CSV/XLSX)", accept = c("xlsx", "csv")),
              checkboxInput("calcLD", "Perform linkage disequilibrium test?", value = FALSE),
              conditionalPanel(
                condition = "input.calcLD",
                helpText("The Markov Chains are automatically set to 10,000 steps")
              ),
              
              checkboxInput("calcHWE", "Perform exact test of HWE?", value = FALSE),
              conditionalPanel(
                condition = "input.calcHWE",
                helpText("The Markov Chains are automatically set to 1,000,000 steps")
              ),
              
              conditionalPanel(
                condition = "input.calcHWE || input.calcLD",
                helpText("Running will take some time and will take longer with more populations/samples.")
              ),

              actionButton("runArlecore", "Run Arlecore", icon = icon("arrow-up-right-from-square")),
              uiOutput("download_arlecore_results_UI")
            ),
            tabBox(
              tabPanel(
                title = "Instructions",
                h4("Calculate common population statistics using Arlecore (terminal-based version of Arlequin)"),
                p(strong("Input file/s:"), "CSV file containing marker and population data.
                                              Each row should represent multi-locus data for an individual sample."),
                p(strong("Expected output file/s:"), ".xlsx and .ars file"),
                hr(),
                p(tags$a("Arlequin",
                  href = "https://cmpg.unibe.ch/software/arlequin35/",
                  target = "_blank"
                ), " is a free software package for population genetic analysis, calculating intra-population and inter-population metrices (Excoffier & Lischer, 2010).")
              ),
              tabPanel(
                title = "Sample Input Format/s",
                DT::dataTableOutput("exampleForArlecore")
              ),
              tabPanel(
                title = "Download Sample Files",
                h4("Sample File"),
                tags$ul(
                  tags$a("Sample CSV file", href = "sample.csv", download = "sample.csv")
                )
              )
            ),
            fluidRow(
              tabBox(
                title = "Results",
                width = 12,
                tabPanel(
                  title = "Expected Heterozygosity",
                  DT::DTOutput("hwe_arlecore"),
                  plotly::plotlyOutput("hwe_arlecore_plot"),
                  uiOutput("hwe_arlecore_plots")
                ),
                tabPanel(
                  title = "FST Matrix",
                  DT::DTOutput("fst_arlecore"),
                  br(),
                  plotly::plotlyOutput("fst_heatmap_plot_arlequin"),
                  br(),
                  plotOutput("fst_pairwise_heatmap_plot_overlap"),
                  br(),
                  plotly::plotlyOutput("fst_pairwise_heatmap_plot", height = "600px")
                ),
                tabPanel(
                  title = "Coancestry Coefficient",
                  DT::DTOutput("coancestry_arlecore"),
                  plotly::plotlyOutput("coancestry_heatmap_plot")
                ),
                tabPanel(
                  title = "Diversity and HWE calculations",
                  uiOutput("population_tables")
                ),
                tabPanel(
                  title = "Loci in LD",
                  DT::DTOutput("ld_tables")
                )
              )
            )
          )
        )
      )
    )
  )
}
