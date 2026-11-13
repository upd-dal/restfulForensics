msa_tab <- function() {
  # BARCODING MODULE ================================================================================
  tabItem(
    tabName = "MSAtab",
    # Alignment submodule =====================================================================
    tabsetPanel(
      tabPanel(
        title = "Multiple Sequence Alignment",
        profvis::profvis_ui("profiler_msa"),
        fluidRow(
          box(
            fileInput("fastaFile", "Upload zipped FASTA files", accept = c(".zip", ".tar")),
            radioButtons("substitutionMatrix", "Choose Substitution Matrix for MSA",
                         choices = c("ClustalW" = "ClustalW", "ClustalOmega" = "ClustalOmega", "MUSCLE" = "Muscle")
            ),
            actionButton("runMSA", "Align", icon = icon("align-justify")),
            br(),
            selectInput("msaDownloadType", "Choose alignment (FASTA and PDF) version to download:",
                        choices = c(
                          "Initial" = "initial",
                          "Adjusted" = "adjusted",
                          "Staggered" = "staggered"
                        ),
                        selected = "initial"
            )
          ),
          tabBox(
            tabPanel(
              title = "Instructions",
              h4("Perform multiple sequence alignment"),
              p(strong("Input file/s:"), "Zipped folder (.zip) containing FASTA files."),
              p(strong("Parameter/s:"), "Substitution matrix for the alignment (ClustalW, ClustalOmega, MUSCLE)"),
              p(strong("Expected output/s:")),
              tags$ul(
                tags$li("Aligned sequences"),
                tags$li("Alignment scores"),
                tags$li("Alignment PDF")
              ),
              hr(),
              p("This performs multiple sequence alignment using the msa R package (Bodenhofer et al., 2015)
                                               and post-processing using DECIPHER (Wright, 2015)."),
              p(
                "The aligned sequences can be used in the tab ",
                tags$a(actionLink("tophylogentab", "'Phylogenetic Tree Analysis'."))
              )
            ),
            tabPanel(
              title = "Download Sample File",
              p("The zipped data contains pangolin case samples (Luzcon et al., 2016)."),
              tags$ul(
                tags$a("Sample zipped FASTA file", href = "pangolins.zip", download = "pangolins.zip")
              ),
            )
          ),
          tabBox(
            title = "Results",
            width = 12,
            tabPanel(
              title = "Preview of Alignments",
              verbatimTextOutput("initialAlignmentText"),
              br(),
              verbatimTextOutput("adjustedAlignmentText"),
              br(),
              verbatimTextOutput("staggeredAlignmentText"),
              br(),
              verbatimTextOutput("alignmentScoresPreview")
            ),
            tabPanel(
              title = "Download Results",
              p("Note: Alignment best visualized in the PDF file."),
              uiOutput("downloadAlignedFASTA_UI"),
              uiOutput("downloadAlignmentScores_UI"),
              downloadButton("downloadAlignmentPDF", "Download Alignment in PDF")
            )
          )
        )
      ),
      tabPanel(
        title = "Trim Sequences",
        profvis::profvis_ui("profiler_trim"),
        box(
          checkboxInput("uploadMSATrim", "Use results from the MSA tab?", value = FALSE),
          conditionalPanel(
            condition = "input.uploadMSATrim == true",
            selectInput("sequenceToTrim", "Alignment to trim",
                        choices = c(
                          "Initial" = "initial",
                          "Adjusted" = "adjusted",
                          "Staggered" = "staggered"
                        ),
                        selected = "initialTrim"
            )
          ),
          conditionalPanel(
            condition = "input.uploadMSATrim == false",
            fileInput("fileForTrim", "Upload MSA file for trimming")
          ),
          numericInput("startSeq", "Starting position", value = 50),
          numericInput("endSeq", "End position", value = 150),
          helpText("Indicate the start and end position of the sequence to keep"),
          actionButton("trimSequences", "Trim sequences", icon=icon("scissors")),
          uiOutput("downloadTrimmedFASTAUI")
          ),
        tabBox(
          tabPanel(
            title = "Instructions",
            h4("Trim sequence alignments"),
            p(strong("Input file/s:"), "Alignment file (msa, fasta, fas, fa, aln, msf)"),
            p(strong("Parameter/s:"), "Start and end position of the subsequence to keep."),
            p(strong("Expected output/s: Fasta file of trimmed sequences.")),
            hr(),
            p("The ape and phangorn packages for building phylogenetic trees require
              a distance matrix built from MSA without missing values, hence the
              aligned sequences should be trimmed as needed.")
            )
        )
      )
    ) # end of tabset panel
  )
}

phylogeny_tab <- function() {
  # Phylogenetic analysis submodule =================================================================
  tabItem(
    tabName = "PhylogenAnalysis",
    tabPanel(
      title = "Phylogenetic Tree Analysis",
      profvis::profvis_ui("profiler_tree"),
      fluidRow(
        box(
          checkboxInput("uploadMSA", "Use results from the MSA tab?", value = FALSE),
          conditionalPanel(
            condition = "input.uploadMSA == true",
            selectInput("treeAlignmentType", "Use alignment for tree construction:",
              choices = c(
                "Initial" = "initial",
                "Adjusted" = "adjusted",
                "Staggered" = "staggered",
                "Trimmed" = "trimmed"
              ),
              selected = "initial"
            )
          ),
          conditionalPanel(
            condition = "input.uploadMSA == false",
            fileInput("msaFileforPhylogen", "Upload MSA file")
          ),
          
          # 
          actionButton("validateAlignment", "Validate Input File (distance)", icon = icon("check")),
          uiOutput("phylogeneticOptionsUI")
        ),
        tabBox(
          title = "Instructions",
          h4("Perform phylogenetic tree reconstruction"),
          p(strong("Input file"), "is a multiple sequence alignment. Results from the 'MSA' tab are also accepted.
                                     If using outputs from the 'MSA' tab, there is an option to use the raw, adjusted, or staggered alignment for tree construction.
                                  Accepted MSA formats: .msa, .fasta, .msf, .aln, .faa, .fas"),
          p(strong("Parameters"), "vary based on the method."),
          p(strong("Expected output"), "is the phylogenetic tree in PNG format."),
          hr(),
          p(
            "This performs phylogenetic tree reconstruction using ape (Paradis and Strimmer, 2004) and phangorn (Schliep, 2011) R packages
                                      on multiple sequence alignments. Outputs generated from the ",
            tags$a(actionLink("tomsatab", "'Multiple Sequence Alignment'")), " tab are also accepted."
          ),
          p(
            "Approaches to tree construction are Neighbor Joining (NJ), Unweight Pair Group Method using Arithmetic averages (UPGMA), Maximum Parsimony, and Maximum Likelihood. Check the assumptions and constraints of each approach ",
            tags$a("(Zou et al., 2024)",
              href = "https://pmc.ncbi.nlm.nih.gov/articles/PMC11117635/",
              target = "_blank"
            ), "."
          )
        ),
        tabBox(
          title = "Results",
          width = 12,
          tabPanel(
            title = "View Results",
            h4("Phylogenetic Tree"),
            uiOutput("downloadTree_UI"),
            uiOutput("treeImage")
          )
        )
      )
    )
  )
}

barcoding_tab <- function() {
  tabItem(
    tabName = "BarcodingTab",
    tabPanel(
      title = "Barcoding",
      profvis::profvis_ui("profiler_barcoding"),
      tabsetPanel(
        # Species Identification submodule =========================================
        tabPanel(
          title = "Species Identification",
          fluidRow(
            box(
              fileInput("refBarcoding", "Upload Aligned Reference Sequences"),
              fileInput("queBarcoding", "Upload Aligned Query Sequences"),
              helpText("The reference and query sequences should have the same length."),
              checkboxInput("kmerSelect", "Use k-mer method?", value = FALSE),
              conditionalPanel(
                "input.kmerSelect == false",
                selectInput("barcodingMethod", "Select method to train model and infer membership:",
                  choices = c("fuzzyId", "bpNewTraining", "bpNewTrainingOnly", "bpUseTrained", "Bayesian"),
                  selected = "bpNewTraining"
                )
              ),
              conditionalPanel(
                "input.kmerSelect == true",
                radioButtons("kmerType", "Choose Method",
                  choices = c("Fuzzy-set Method and kmer", "BP-based Method and kmer")
                ),
                conditionalPanel(
                  "input.kmerType == 'Fuzzy-set Method and kmer'",
                  numericInput("kmerValueFuzzy", "K-mer value", value = 1, min = 0),
                  checkboxInput("optimizationKMER", "Use different kmer length?", value = FALSE)
                ),
                conditionalPanel(
                  "input.kmerType == 'BP-based Method and kmer'",
                  numericInput("kmerValueBP", "K-mer value", value = 1, min = 0),
                  checkboxInput("builtModel", "Use built model", value = FALSE),
                  numericInput("lrValue", "Parameter for weight decay", value = 0.00005),
                  numericInput("maxitValue", "Maximum number of iterations", value = 1000000)
                )
              ),
              actionButton("identifySpecies", "Identify Species", icon = icon("magnifying-glass"))
            ),
            tabBox(
              tabPanel(
                title = "Instructions",
                h4("Perform species identification using the R package 'BarcodingR' (Zhang et al., 2016)"),
                p(strong("Input file/s:")),
                tags$ul(
                  tags$li("Aligned reference sequences (.msa, .fasta, .msf, .aln, .faa, .fas)"),
                  tags$li("Aligned query sequences (.msa, .fasta, .msf, .aln, .faa, .fas)")
                ),
                p(strong("Parameter/s:")),
                tags$ul(
                  tags$li("(without kmer method) Training model: bpNewTraining, fuzzyId, bpNewTrainingOnly, bpUsedTrained, or Bayesian"),
                  tags$li("(with kmer method) Fuzzy-set Method or BP-based method")
                )
              ),
              tabPanel(
                title = "Sample Input Format/s",
                p("The sample alignment files are generated from the 'MSA' tab
                  using the pangolins (Luczon et al., 2016) case dataset aligned using ClustalW."),
                tags$ul(
                  tags$a("A. Sample files", href = "pangolins_samples.fa", download = "pangolins_samples.fa"),
                  br(),
                  tags$a("B. Sample reference files", href = "pangolins_ref.fa", download = "pangolins_ref.fa")
                )
              )
            ),
            box(
              title = "Results",
              width = 12,
              p("Species Identification Success Rate"),
              verbatimTextOutput("identificationResult"),
              br(),
              p("Species inferred and corresponding confidence levels:"),
              DT::DTOutput("identified_samples")
            )
          )
        ),
        # Optimization of kmer values submodule ====================================
        tabPanel(
          title = "Optimize kmer values",
          fluidRow(
            box(
              title = "Optimization Options",
              width = 6,
              fileInput("optimizeKmerRef", "Upload reference dataset"),
              numericInput("maxKmer", "Length of maximum kmer value", value = 5, min = 2),
              actionButton("calOptimumKmer", "Identify Optimum kmer value", icon = icon("upload"))
            ),
            tabBox(
              title = "Instructions",
              width = 6,
              tabPanel(
                title = "Overview",
                h4("Calculate the optimal kmer values using BarcodingR (Zhang et al., 2016)"),
                p(strong("Input file/s:"), "Aligned sequences (.msa, .fasta, .msf, .aln, .faa, .fas)"),
                p(strong("Parameter/s:"), "Length of maximum kmer value"),
                p(strong("Expected output file:"), "Kmer plot (.png)")
              )
            )
          ),
          fluidRow(
            tabBox(
              title = "Results",
              width = 12,
              tabPanel(
                title = "Outputs",
                verbatimTextOutput("kmerResult"),
                imageOutput("kmerPlot"),
                uiOutput("downloadKmerPlot_UI")
              )
            )
          )
        ),

        # Barcoding gap submodule ==================================================
        tabPanel(
          title = "Barcoding Gap",
          fluidRow(
            box(
              title = "Gap Calculation Options",
              width = 6,
              fileInput("barcodeRef", "Upload reference dataset"),
              selectInput("gapModel", "Choose Distance",
                choices = c("raw", "K80", "euclidean"),
                selected = "raw"
              ),
              actionButton("gapBarcodes", "Calculate gap", icon = icon("arrows-left-right-to-line"))
            ),
            tabBox(
              title = "Instructions",
              width = 6,
              tabPanel(
                title = "Overview",
                h4("Calculate the barcoding gap using BarcodingR (Zhang et al., 2016)"),
                p(strong("Input file:"), "Aligned reference sequences (.msa, .fasta, .msf, .aln, .faa, .fas)."),
                p(strong("Parameter/s:"), "Distance (raw, K80, euclidean)"),
                p(strong("Expected output file:"), "Barcoding gap plot (.png)")
              )
            )
          ),
          fluidRow(
            tabBox(
              title = "Results",
              width = 12,
              tabPanel(
                title = "Outputs",
                verbatimTextOutput("barcodingResult"),
                imageOutput("BarcodingGapPlot"),
                uiOutput("downloadGapPlot_UI")
              )
            )
          )
        ),

        # Evaluation of barcodes submodule =========================================
        tabPanel(
          title = "Evaluate Barcodes",
          fluidRow(
            box(
              fileInput("barcode1", "Upload Barcode 1"),
              fileInput("barcode2", "Upload Barcode 2"),
              numericInput("kmer1", "Length of kmer for barcode 1", value = 5, min = 1),
              numericInput("kmer2", "Length of kmer for barcode 2", value = 5, min = 1),
              actionButton("evalBarcodes", "Evaluate Barcodes", icon = icon("code-compare"))
            ),
            tabBox(
              title = "Instructions",
              h4("Evaluate barcodes using species identification success rate criteria (Zhang et al., 2016)"),
              p(strong("Input file/s:"), "Aligned sequences (.msa, .fasta, .msf, .aln, .faa, .fas)."),
              p(strong("Parameter/s:"), "Length of kmer for barcode 1 and barcode 2 (separate)")
            ),
            tabBox(
              title = "Results",
              tableOutput("evalBarcodesResult")
            )
          )
        ),

        # Calculation of species membership value using TDR submodule ==============
        tabPanel(
          title = "Species Membership Value (TDR)",
          fluidRow(
            box(
              p("Calculate the TDR2 value"),
              fileInput("oneSpe", "Upload DNA seq from a single query species"),
              fileInput("queSpe", "Upload DNA seq from different samples"),
              numericInput("bootValue1", "Bootstrap value for query species", value = 10, min = 1),
              numericInput("bootValue2", "Bootstrap value for reference samples", value = 10, min = 1),
              actionButton("calculateTDR2", "Calculate", icon = icon("calculator"))
            ),
            tabBox(
              title = "Instructions",
              h4("Calculate the Species Membership Value in terms of
                                                  Two-Dimensional non-parametric resampling (TDR) using BarcodingR (Zhang et al., 2016)"),
              p(strong("Input file/s:"), "Aligned sequences (.msa, .fasta, .msf, .aln, .faa, .fas)."),
              p(strong("Parameter/s:"), "Bootstrap value for query and reference samples.")
            ),
            tabBox(
              title = "Results",
              tableOutput("tdrValues")
            )
          )
        )
      )
    )
  )
}
