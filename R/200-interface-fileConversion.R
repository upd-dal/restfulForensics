file_conversion_tab <- function() {
  # FILE CONVERSION MODULE===========================================================================
  tabItem(
    tabName = "FileConv",
    tabsetPanel(
      # Format conversion submodule ==========================================================
      tabPanel(
        title = "Convert Files",
        fluidRow(
          box(
            title = "File Conversion Options",
            width = 5,
            radioButtons("inputType1", "Choose starting file type",
              choices = c(
                "VCF file" = "vcf1",
                "BCF file" = "bcf1",
                "PLINK files (.bed/.bim/.fam)" = "plink1",
                "CSV/Excel file" = "csv1"
              )
            ),
            conditionalPanel(
              condition = "input.inputType1 == 'vcf1'",
              fileInput("VCFFile", "Upload VCF File", accept = c(".vcf", ".zip", ".tar")),
              radioButtons("inputType2_vcf", "Choose final file type",
                choices = c(
                  "PLINK2 files (.psam/.pvar/.pgen)" = "plink2",
                  "PLINK1.9 files (.bed/.bim/.fam)" = "plink1",
                  "CSV file" = "csv2"
                )
              )
            ),
            conditionalPanel(
              condition = "input.inputType1 == 'bcf1'",
              fileInput("BCFFile", "Upload BCF File", accept = c(".bcf", ".zip", ".tar")),
              radioButtons("inputType2_bcf", "Choose final file type",
                choices = c(
                  "VCF file" = "vcf2",
                  "PLINK2 files (.psam/.pvar/.pgen)" = "plink2",
                  "PLINK1.9 files (.bed/.bim/.fam)" = "plink1",
                  "CSV file" = "csv2"
                )
              )
            ),
            conditionalPanel(
              condition = "input.inputType1 == 'plink1'",
              fileInput("bedFile", "Upload BED File", accept = c(".bed")),
              fileInput("bimFile", "Upload BIM File", accept = c(".bim")),
              fileInput("famFile", "Upload FAM File", accept = c(".fam")),
              radioButtons("inputType2_plink", "Choose final file type",
                choices = c(
                  "VCF file" = "vcf2",
                  "PLINK2 files (.psam/.pvar/.pgen)" = "plink2",
                  "CSV file" = "csv2"
                )
              )
            ),
            conditionalPanel(
              condition = "input.inputType1 == 'csv1'",
              p("This feature automatically converts a CSV file to VCF"),
              fileInput("CSVFile", "Upload File (CSV/Excel file)", accept = c(".csv", "xlsx", "xlsm", "xlsb", "xls")),
              fileInput("lociMetaFile", "Upload loci/marker information (CSV/Excel file)", accept = c("xlsx", "xlsm", "xlsb", "xls", ".csv"))
            ),
            actionButton("ConvertFILES", "Convert files", icon = icon("file-csv"))
          ),
          tabBox(
            tabPanel(
              title = "Instructions",
              h4("Interconvert common genetic files"),
              p(strong("Input file/s:")),
              p("Required: VCF, BCF, or PLINK (.bed, .bim, .fam) files. It also accepts zipped files as long as it contains the same file type, except if using PLINK files."),
              p("Optional input file:"),
              tags$ul(
                tags$li("(CSV/Excel to VCF) Marker information with the following columns: [1] SNP, [2] CHR, [3] POS, [4] Genetic distance, [5] REF Allele [6] ALT Allele.")
              ),
              p(strong("Expected output file/s:"), "VCF, PLINK-associated files, or CSV."),
              br(),
              p("Maximum accepted file size: 5GB. It is recommended to split files with sizes larger than 5GB into multiple smaller files.")
            ),
            tabPanel(
              title = "Sample Input Format/s",
              h4("For CSV/Excel File to VCF conversion, a separate file on marker information is needed."),
              h4("Required marker info format:"),
              DT::dataTableOutput("markerInfoFormat")
            ),
            tabPanel(
              title = "Download Sample Files",
              tags$a("A. Sample VCF", href = "sample_hgdp.vcf", download = "sample_hgdp.vcf"),
              br(),
              tags$a("B. Sample zipped file (VCF files)", href = "vcf_sample_files.zip", download = "vcf_sample_files.zip"),
              br(),
              tags$a("C. Sample CSV file (for VCF conversion)", href = "sample.csv", download = "sample.csv"),
              br(),
              tags$a("D. Sample marker metadata file (for CSV-VCF conversion)", href = "marker_info.csv", download = "marker_info.csv")
            )
          )
        ),
        fluidRow(
          tabBox(
            title = "Results",
            width = 12,
            tabPanel(
              title = "Download Results",
              uiOutput("downloadVCF_UI"),
              uiOutput("downloadPLINK_UI"),
              uiOutput("downloadtoCSV_UI")
            ),
            tabPanel(
              title = "Genotype table (to CSV)",
              DT::dataTableOutput("toCSVtable")
            )
          )
        )
      ),
      tabPanel(
        title = "Add Metadata",
        fluidRow(
          box(
            fileInput("genotypeFile", "Genotype File or zipped files", accept = c(".vcf", ".bcf", ".gz", ".zip", ".tar")),
            radioButtons("populationType", "Do samples come from a single population?",
              choices = c("Yes" = "single", "No" = "multiplepop")
            ),
            conditionalPanel(
              condition = "input.populationType == 'multiplepop'",
              fileInput("refMetadata", "Upload file containing the metadata", accept = c(".csv", "xlsx", "xlsm", "xlsb", "xls")),
              uiOutput("metaHeader"),
              helpText("*Accepts CSV or Excel files"),

              # --- Breakdowns
              radioButtons("breakdownPop", "Calculate population breakdown?",
                choices = c("Yes" = "YesBreakdown", "No" = "NoBreakdown"), selected = "No"
              ),
              conditionalPanel(
                condition = "input.breakdownPop == 'YesBreakdown'",
                helpText("Specify column name to serve as a basis for the summary count."),
                uiOutput("selectMetaHeader")
              )
            ),
            conditionalPanel(
              condition = "input.populationType == 'single'",
              textAreaInput("typePop_meta", "Enter population", rows = 1)
            ),
            actionButton("addMetadata", "Add Metadata", icon = icon("file-circle-plus"))
          ),
          tabBox(
            tabPanel(
              title = "Instructions",
              h4("Merge select metadata to genotype/sample data"),
              p(strong("Input file/s:")),
              tags$ul(
                tags$li("Genotype data (.vcf, .vcf.gz, .bcf, or PLINK 1.9 files) or CSV/Excel genotype data. Zipped files are also accepted as long as it consists a single file type."),
                tags$li("Sample metadata (CSV/Excel file). Select column names to be merged with the genotype data. Ensure sample IDs are the same.")
              ),
              p(strong("Expected output file/s:"), "CSV file of sample with metadata.")
            ),
            tabPanel(
              title = "Sample Input Format/s",
              h4("To convert to a CSV file with population metadata:"),
              DT::dataTableOutput("ExampleRefFile")
            ),
            tabPanel(
              title = "Download Sample Files",
              tags$a("A. Sample VCF", href = "sample_hgdp.vcf", download = "sample_hgdp.vcf"),
              br(),
              tags$a("B. Sample zipped file (VCF files)", href = "vcf_sample_files.zip", download = "vcf_sample_files.zip")
            )
          )
        ),
        fluidRow(
          tabBox(
            title = "Results",
            width = 12,
            tabPanel(
              title = "Preview File and Download",
              div(
                style = "overflow-x: auto;",
                DT::dataTableOutput("previewTable")
              ),
              br(),
              uiOutput("downloadCSV_UI")
            ),
            tabPanel(
              title = "Population Breakdown",
              div(
                style = "overflow-x: auto;",
                DT::dataTableOutput("previewTableBreakdown")
              ),
              br(),
              uiOutput("downloadBreakdown")
            ),
            tabPanel(
              title = "Samples without Metadata",
              div(
                style = "overflow-x: auto;",
                DT::dataTableOutput("previewMissingData")
              ),
              uiOutput("downloadMissing")
            )
          )
        )
      ),

      # Widen long genotype file submodule ===================================================
      tabPanel(
        title = "Widen SNP calls",
        fluidRow(
          box(
            fileInput("uas_zip", "Upload ZIP or TAR file",
              accept = c(".zip", ".tar")
            ),
            helpText("*Accepts compressed files containing CSV/XLSX files."),
#            fileInput("ref_file", "Optional Reference File (CSV or XLSX)",
#              accept = c(".csv", ".xlsx", ".zip", ".tar")
#            ),
            actionButton("run_uas2csv", "Run Conversion")
          ),
          tabBox(
            tabPanel(
              title = "Instructions",
              h4("Convert zipped files (CSV/Excel file) containing SNP calls in long format into a single .xslx file in a wide format"),
              p(strong("Input file/s:"), "Compressed folder (.zip or .tar) containing CSV/Excel files."),
              p(strong("Expected output file/s:"), "Single CSV file (merged CSV/Excel files).")
            ),
            tabPanel(
              title = "Sample Input Format/s",
              h4("Sample input file. All alleles of available SNPs per sample are listed in a long format."),
              tableOutput("exampleXLSX")
            )
          )
        ),
        fluidRow(
          tabBox(
            title = "Results",
            width = 12,
            tabPanel(
              title = "Preview and Download",
              div(
                style = "overflow-x: auto;",
                DT::dataTableOutput("previewTableUAS")
              ),
              br(),
              uiOutput("downloadUAScsv_UI")
            )
          )
        )
      ),

      # Convert file to SNIPPER compatible file submodule ====================================
      tabPanel(
        title = "To SNIPPER file",
        fluidRow(
          box(
            width = 6,
            fileInput("convertFile", "Upload Genotype File", accept = c("xlsx", "xlsm", "xlsb", "xls", ".csv")),
            checkboxInput("refProvided", "Reference populations included in the file?", value = TRUE),
            helpText("See 'sample input formats' for guidance. If 'TRUE', assuming population metadata is included."),
            conditionalPanel(
              condition = "input.refProvided == false",
              fileInput("refFile", "Upload Reference File", accept = c("xlsx", "xlsm", "xlsb", "xls", ".csv"))
            ),
            checkboxInput("targetPop", "Classify a Target Population?", value = FALSE),
            conditionalPanel(
              condition = "input.targetPop == true",
              textInput("targetPopName", "Target Population Name"),
              helpText("Indicate population name to classify. Population should exist in the input file.")
            ),
            actionButton("convertBtn", "Convert Format", icon = icon("arrow-up-right-from-square"))
          ),
          tabBox(
            tabPanel(
              title = "Instructions",
              h4("Convert file into a SNIPPER-compatible input file for individual classification using ancestry-informative markers"),
              p(strong("Input file/s:"), "CSV or Excel file"),
              p(strong("Parameter/s:"), "(optional) Target population name for classification."),
              p(strong("Expected output file/s:"), "Excel (.xslx) file"),
              hr(),
              p(
                "The ",
                tags$a("SNIPPER app suite",
                  href = "https://mathgene.usc.es/snipper/index.php",
                  target = "_blank"
                ), " is a web portal for classification of individuals into phenotypes and biogeographical ancestry (BGA), developed and hosted by Universidade de Santiago de Compostela. SNIPPER is a companion site to several ",
                tags$a("publications.",
                  href = "https://mathgene.usc.es/snipper4/papers.php",
                  target = "_blank"
                )
              )
            ),
            tabPanel(
              title = "Sample Input Format/s",
              h4("Sample Input File"),
              p("Format if input file does", strong("not"), "contain population metadata:"),
              DT::dataTableOutput("exampleTableSnipper1"),
              br(),
              p("Format if input file", strong("contains"), "population metadata:"),
              DT::dataTableOutput("exampleTableSnipper2"),
              br(),
              h4("Sample reference file"),
              DT::dataTableOutput("exampleRefSnipper")
            ),
            tabPanel(
              title = "Download Sample File",
              tags$ul(
                tags$a("A. Sample CSV file", href = "sample.csv", download = "sample.csv"),
                br(),
                tags$a("B. Sample zipped file", href = "sample_snipper.csv", download = "sample_snipper.csv")
              )
            )
          )
        ),
        fluidRow(
          tabBox(
            title = "Results",
            width = 12,
            tabPanel(
              title = "Preview and Download",
              div(
                style = "overflow-x: auto;",
                DT::dataTableOutput("previewTableSNIPPER")
              ),
              br(),
              uiOutput("downloadSNIPPER")
            )
          )
        )
      ),

      # Convert file to standard STRUCTURE compatible file ===================================
      tabPanel(
        title = "To STRUCTURE file",
        fluidRow(
          box(
            width = 6,
            fileInput("tostrFile", "Upload Genotype file", accept = c("xlsx", "xlsm", "xlsb", "xls", ".csv")),
            helpText("Use the 'Convert files' submodule if using VCF, BCF, or PLINK files. Population data is necessary."),
            radioButtons("systemFile", "Choose the operating system where STRUCTURE v2.3.4 is installed",
              choices = c("Linux" = "Linux", "Windows" = "Windows")
            ),
            actionButton("csv2str", "Generate File", icon = icon("arrow-up-right-from-square"))
          ),
          tabBox(
            tabPanel(
              title = "Instructions",
              h4("Convert file into a standard STRUCTURE v2.3.4-compatible file"),
              p(strong("Input file/s:"), "CSV/Excel file containing marker and population data.
                                              Each row should represent multi-locus data for an individual sample."),
              p(strong("Parameter/s:"), "User's operating system (for STRUCTURE input compatibility)"),
              p(strong("Expected output file/s:")),
              tags$ul(
                tags$li("structure (.str) file"),
                tags$li("revised input file")
              ),
              hr(),
              p(tags$a("STRUCTURE",
                href = "https://web.stanford.edu/group/pritchardlab/structure.html",
                target = "_blank"
              ), " is a free software package for investigating population structure using
                                              multi-locus genotype data (Stephen and Donnelly, 2000; Falush et al., 2003;
                                              Falush et al., 2007; Hubisz et al., 2009)."),
              p("STRUCTURE generally can't handle sample labels with alphabets, the function converts sample labels to their associated row number."),
              p(
                "For users who opt to use STRUCTURE via the terminal or GUI, instructions can be found here: ",
                tags$a("STRUCTURE v2.3.4 documentation",
                  href = "https://web.stanford.edu/group/pritchardlab/structure_software/release_versions/v2.3.4/html/structure.html",
                  target = "_blank"
                )
              )
            ),
            tabPanel(
              title = "Sample Input Format/s",
              DT::dataTableOutput("examplePop_STRUI")
            ),
            tabPanel(
              title = "Download Sample File",
              tags$a("Sample file", href = "sample.csv", download = "sample.csv")
            )
          )
        ),
        fluidRow(
          tabBox(
            title = "Results",
            width = 12,
            tabPanel(
              title = "Preview .str file",
              DT::DTOutput("strFile"),
              uiOutput("downloadSTRfile_UI")
            ),
            tabPanel(
              title = "Revised Input file",
              DT::DTOutput("revisedCSV"),
              uiOutput("downloadrevised_UI")
            )
          )
        )
      ),

      # Convert file to Arlequin compatible file ===================================
      tabPanel(
        title = "To Arlequin file",
        fluidRow(
          box(
            width = 6,
            fileInput("toArleFile", "Upload Genotype file", accept = c("xlsx", "xlsm", "xlsb", "xls", ".csv")),
            helpText("Use the 'Convert files' submodule if using VCF, BCF, or PLINK files. Population data is necessary."),
            checkboxInput("genotypicData", "Using Genotypic Data? (Untick if haplotypic data)", value = TRUE),
            checkboxInput("gameticPhase", "Gametic Phase known?", value = FALSE),
            checkboxInput("recessiveData", "Using co-dominant data?", value = FALSE),
            textInput("locusSep", "Symbol that separates alleles", value = "/"),
            # selectInput("dataTypeArlecore", "Type of Data to be Analyzed",
            #            choices = c(
            #               "DNA" = "DNA",
            #               "STANDARD" = "STANDARD",
            #               "FREQUENCY" = "FREQUENCY",
            #               "MICROSAT" = "MICROSAT",
            #               "RFLP" = "RFLP"
            #            ),
            #            selected = "STANDARD"),
            actionButton("convert2Arle", "Generate File", icon = icon("arrow-up-right-from-square"))
          ),
          tabBox(
            tabPanel(
              title = "Instructions",
              h4("Convert genotype and population data (CSV/Excel) to Arlequin-compatible file"),
              p(strong("Input file/s:"), "File (CSV/Excel) containing marker and population data.
                                            Each row should represent multi-locus data for an individual sample."),
              p(strong("Expected output file/s:"), ".arp file"),
              hr(),
              p(tags$a("Arlequin",
                href = "https://cmpg.unibe.ch/software/arlequin35/",
                target = "_blank"
              ), " is a free software package for population genetic analysis, calculating intra-population and inter-population metrices (Excoffier & Lischer, 2010).")
            ),
            tabPanel(
              title = "Sample Input Format/s",
              DT::dataTableOutput("exampleForArlecore_UI")
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
              title = "Download Results",
              uiOutput("downloadArpFile_UI")
            )
          )
        )
      )
    )
  )
}
