snp_extraction_tab <- function() {
  # SNP EXTRACTION MODULE ===========================================================================
  tabItem(
    tabName = "markerExtract",
    tabsetPanel(
      # SNP extraction submodule =============================================================
      tabPanel(
        title = "SNP Data Extraction",
        fluidRow(
          box(
            width = 5,
            radioButtons("inputFileType", "A. Input file format",
              choices = c("VCF/VCF.GZ/BCF", "PLINK"), inline = TRUE
            ),
            conditionalPanel(
              condition = "input.inputFileType == 'VCF/VCF.GZ/BCF'",
              fileInput("markerFile", "Upload Genotype File", accept = c(".vcf", ".bcf", ".vcf.gz"))
            ),
            conditionalPanel(
              condition = "input.inputFileType == 'PLINK'",
              fileInput("bedFile", "PLINK BED file", accept = c(".bed")),
              fileInput("bimFile", "PLINK BIM file", accept = c(".bim")),
              fileInput("famFile", "PLINK FAM file", accept = c(".fam")),
              helpText("If multiple PLINK files will be used, upload as a zipped file. Note that files will be merged."),
              fileInput("zippedPLINK", "Zipped PLINK Files", accept = c(".zip", ".tar"))
            ),
            actionButton("validateBtn", "Validate Input File Format", icon = icon("check")),
            uiOutput("markerOptionsUI")
          ),
          tabBox(
            tabPanel(
              title = "Instructions",
              h4("Extract SNPs based on reference SNP cluster ID (rsID) or GRCh37/GRCh38 position"),
              p(strong("Input file/s:")),
              p("(1) VCF, BCF, or PLINK (.bed, .bim, .fam) files."),
              p("(2) Markers/position list — type rsIDs manually, upload a list, or use a POS CSV/Excel/Text file."),
              p("Position list format:"),
              tags$ul(
                tags$li("[1] (optional) rsID/marker name"),
                tags$li("[2] Chromosome number"),
                tags$li("[3] Position (bp)")
              ),
              p(strong("Expected output file:"), "VCF File"),
              br(),
              p("Maximum accepted file size: 5GB. It is recommended to split files with sizes larger than 5GB into multiple smaller files.")
            ),
            tabPanel(
              title = "Sample Input Format/s",
              h4("rsID Format"),
              tableOutput("examplersID"),
              h4("Position Format"),
              tableOutput("examplePOS"),
              helpText("rsID column is required if rsID will be added to output.")
            ),
            tabPanel(
              title = "Download Sample Files",
              tags$ul(
                tags$a("A. Sample VCF file", href = "sample_hgdp.vcf", download = "sample_hgdp.vcf"),
                br(),
                tags$a("B. Sample marker metadata file (with rsID)", href = "marker_info2.csv", download = "marker_info2.csv"),
                br(),
                tags$a("C. Sample marker metadata file (wihout rsID)", href = "marker_noid.csv", download = "marker_noid.csv")
              )
            )
          )
        ),
        fluidRow(
          tabBox(
            title = "Results",
            width = 12,
            tabPanel(
              title = "Detected rsIDs",
              DT::DTOutput("variantTable")
            ),
            tabPanel(
              title = "Extraction Results",
              uiOutput("downloadVCF_UI_Extracted"),
              helpText("Note: File can appear as vCard files, it is still a VCF file.")
            )
          )
        )
      ),

      # Concordance analysis submodule =======================================================
      tabPanel(
        title = "Concordance Analysis",
        fluidRow(
          box(
            fileInput("concordanceFile1", "Upload File A", accept = c("xlsx", "xlsm", "xlsb", "xls", ".csv", ".vcf")),
            fileInput("concordanceFile2", "Upload File B", accept = c("xlsx", "xlsm", "xlsb", "xls", ".csv", ".vcf")),
            checkboxInput("isPhased", "Phased genotypes", value = FALSE),
            actionButton("compareBtn", "Run Concordance Analysis", icon = icon("play"))
          ),
          tabBox(
            tabPanel(
              title = "Instructions",
              h4("Perform concordance analysis between files/datasets with overlapping samples"),
              p(strong("Input file/s:"), "CSV/Excel/VCF files with the same samples. File A and B can either be a CSV, Excel, or VCF file"),
              p(strong("Parameter/s:"), "Indicate if using phased genotypes"),
              helpText("Markers will be directly compared and the order of the alleles is considered when matching for concordance."),
              p(strong("Expected output/s:")),
              tags$ul(
                tags$li("Concordance table (.xlsx)"),
                tags$li("Concordance plot (.png)")
              )
            ),
            tabPanel(
              title = "Sample Input Format/s",
              h4("File Format (for concordance)"),
              tableOutput("exampleTable")
            ),
            tabPanel(
              title = "Download Sample Files",
              tags$ul(
                tags$a("Sample CSV file (1)", href = "sample1_for_concordance.csv", download = "sample1_for_concordance.csv"),
                br(),
                tags$a("Sample CSV file (2)", href = "sample2_for_concordance.csv", download = "sample2_for_concordance.csv")
              )
            )
          )
        ),
        fluidRow(
          tabBox(
            title = "Results",
            width = 12,
            tabPanel(
              title = "Summary Table",
              div(
                style = "overflow-x: auto;",
                DT::dataTableOutput("concordanceResults")
              ),
              br(),
              uiOutput("downloadConcordance_UI")
            ),
            tabPanel(
              title = "Concordance Plot",
              div(
                style = "overflow-x: auto;",
                plotOutput("concordancePlot", height = "600px"),
                br(),
                uiOutput("downloadConcordancePlot_UI")
              )
            )
          )
        )
      )
    )
  )
}
