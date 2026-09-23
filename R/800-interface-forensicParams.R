forensic_params_tab <- function() {
  # CALCULATION OF FORENSIC PARAMETERS MODULE =======================================================
  tabItem(
    tabName = "ForensicParams",
    tabsetPanel(
      tabPanel(
        title = "Population Database",
        fluidRow(
          box(
            fileInput("iisnpsFile", "Upload File", accept = c(".csv", "xlsx", "xlsm", "xlsb", "xls")),
            helpText("See 'Sample Input File' for accepted formats. Frequency table and genotype files are accepted."),

            actionButton("calcIISNPs", "Calculate", icon = icon("calculator")),
            uiOutput("downloadMetrics_UI")
          ),
          tabBox(
            tabPanel(
              title = "Instructions",
              h4("Calculate forensic parameters specific for identity-informative SNPs"),
              tags$ul(
                tags$li("Random match probability (PM)"),
                tags$li("Power of discrimination (PD)"),
                tags$li("Polymorphism Information Content (PIC)"),
                tags$li("Power of Exclusion (PE)"),
                tags$li("Typical Paternity Index (TPI)")
              ),
              p(strong("Input file:"), "CSV or Excel file containing population and genotype information."),
              p(strong("Expected output file:"), "CSV file"),
              hr(),
              p(
                "Guidelines on statistical calculations for casework: ",
                tags$a("Guidelines (for STR):",
                  href = "https://dfs.dc.gov/sites/default/files/dc/sites/dfs/page_content/attachments/FBS22%20-%20STR%20Statistical%20Calculations%20Guidelines.pdf",
                  target = "_blank"
                )
              ),
              p(
                "Guidelines and interpretations: ",
                tags$a("Based on the STRAF book",
                  href = "https://agouy.github.io/straf_book/forensic-parameters.html",
                  target = "_blank"
                )
              ),
            ),
            tabPanel(
              title = "Sample Input File",
              h4("Acceptable file inputs: genotype files or an allele frequency table:"),
              p("Genotype file"),
              DT::dataTableOutput("referenceData_UI"),
              p("Allele frequency table"),
              DT::dataTableOutput("afSample_UI"),
              h4("Sample profile to match"),
              DT::dataTableOutput("profileSample_UI")
            ),
            tabPanel(
              title = "Download Sample Files",
              tags$ul(
                tags$a("Sample CSV file", href = "sample.csv", download = "sample.csv")
              ),
              tags$ul(
                tags$a("Sample Allele Frequency Table", href = "pop_stat.xlsx", download = "pop_stat.xlsx")
              )
            )
          )
        ),
        fluidRow(
          tabBox(
            title = "Results",
            width = 12,
            tabPanel(
              title = "Overall Forensic Params",
              selectInput("selected_pop", "Select Population", choices = NULL),
              div(
                style = "overflow-x: auto;",
                DT::dataTableOutput("popTable")
              )
            ),
            tabPanel(
              title = "Genotype Frequencies",
              selectInput("selected_pop_gt", "Select Population", choices = NULL),
              div(
                style = "overflow-x: auto;",
                DT::dataTableOutput("genotypeFreqs_UI")
              )
            )
          )
        )
      ), 
      tabPanel(
        title = "Calculate RMP for a Profile",
        box(
          checkboxInput("newPopDatabase", "Upload allele frequency for a population?", value = FALSE),
          helpText("If FALSE, the data used for the 'Population Database' tab will be used"),
          conditionalPanel(
            condition = "!input.newPopDatabase",
            uiOutput("rmp_population_UI")
            ),
          conditionalPanel(
            condition = "input.newPopDatabase",
            fileInput("newPopDataFile", "Upload population reference database (XLSX/CSV)", accept = c("xlsx", "xlsm", "xlsb", "xls", ".csv"))
          ),
          fileInput("fileProfile", "Upload Profile", accept = c(".csv", "xlsx", "xlsm", "xlsb", "xls", ".txt")),
          checkboxInput("floorCeiling", "Use 5/2n rule in calculating genotype frequency?", FALSE),
          helpText("Defining the total individuals/samples assumes that the population data uploaded is from a single population
                   with subpopulations; hence the use of theta."),
          conditionalPanel(
            condition = "input.floorCeiling",
            numericInput("totalPop", "Total individuals/samples", value = 100, min = 10, )
          ),
          numericInput("thetaValue", "Theta Value", min = 0, value = 0.01, max = 1),
          actionButton("calcRMP", "Calculate", icon = icon("calculator"))
          
        ),
        tabBox(
          tabPanel(
            title = "Instructions",
            p("Calculate the RMP value for a profile."),
            p(strong("Input file:")),
            tags$ul(
              tags$li("(Reference database) CSV or Excel file containing allele frequency for a set of loci in a population"),
              tags$li("(Reference database) Allele frequencies calculated from the 'Population Database' submodule"),
              tags$li("Sample profile containing genotype information for the individual with one marker per row")
            ),
            p(strong("Expected output file:"), "Calculation results"),
          ),
          tabPanel(
            title = "Sample Input Format/s"
          ),
          tabPanel(
            title = "Download Sample Files",
          )
        ), # end of instructions box
        fluidRow(
          tabBox(
            title = "Results",
            width = 12,
            tabPanel(
              title = "Calculations",
              tableOutput("rmp_summary")
            ),
            tabPanel(
              title = "Locus",
              tableOutput("locus_information")
            ),
            tabPanel(
              title = "Profile Data",
              tableOutput("profileData")
            )
          )
        )
      ) # end of profile
    )
  )
}
