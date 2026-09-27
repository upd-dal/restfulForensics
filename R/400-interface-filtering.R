filtering_tab <- function() {
  # FILTERING MODULE ================================================================================
  tabItem(
    tabName = "FilterTab",
    fluidRow(
      box(
        title = "Upload Files",
        radioButtons("inputFileTypeFilter", "A. Input file format",
          choices = c("VCF/VCF.GZ/BCF", "PLINK"), inline = TRUE
        ),
        conditionalPanel(
          condition = "input.inputFileTypeFilter == 'VCF/VCF.GZ/BCF'",
          fileInput("markerFileFilter", "Upload Genotype File", accept = c(".vcf", ".bcf", ".vcf.gz"))
        ),
        conditionalPanel(
          condition = "input.inputFileTypeFilter == 'PLINK'",
          fileInput("bedFileFilter", "PLINK BED file", accept = c(".bed")),
          fileInput("bimFileFilter", "PLINK BIM file", accept = c(".bim")),
          fileInput("famFileFilter", "PLINK FAM file", accept = c(".fam"))
        ),
        fileInput("highlightRef", "Optional Reference file for highlighting (CSV/Excel file)", accept = c("xlsx", "xlsm", "xlsb", "xls", ".csv")),
        checkboxInput("enableDP", "Plot Depth of Coverage", value = TRUE),
        helpText("Depth of Coverage Plot only available if using a VCF file."),
        selectInput("colorPalette", "Color Palette", choices = rownames(RColorBrewer::brewer.pal.info), selected = "Set2"),
        hr(),
        h4("PLINK 2.0 Filtering Options"),
        checkboxInput("filterIndiv", "Filter Individuals (--mind)", value = FALSE),
        helpText("Exclude individuals with a missing genotype rate greater than the threshold"),
        conditionalPanel(
          "input.filterIndiv == true",
          numericInput("mindThresh", "Missingness Threshold (--mind)", value = 0.1, min = 0, max = 1, step = 0.01)
        ),
        checkboxInput("filterVariant", "Filter Variants (--geno)", value = FALSE),
        helpText("Exclude SNPs with a missing genotype rate greater than the threshold."),
        conditionalPanel(
          "input.filterVariant == true",
          numericInput("genoThresh", "Missingness Threshold (--geno)", value = 0.1, min = 0, max = 1, step = 0.01)
        ),
        checkboxInput("filterAllele", "Filter Variants (--maf)", value = FALSE),
        helpText("Exclude SNPs with a minor allele frequency less than the threshold."),
        conditionalPanel(
          "input.filterAllele == true",
          numericInput("mafThresh", "Minor Allele Frequency Threshold (--maf)", value = 0.1)
        ),
        checkboxInput("filterQuality", "Filter by Quality (--qual-threshold)", value = FALSE),
        helpText("Exclude variants with quality scores below the threshold."),
        conditionalPanel(
          "input.filterQuality == true",
          numericInput("qualThresh", "Quality Score Threshold (--qual-threshold)", value = 5)
        ),
        checkboxInput("filterHWE", "Filter Variants (--hwe)", value = FALSE),
        helpText("Exclude SNPs deviating from the Hardy-Weinberg Equilibrium."),
        conditionalPanel(
          "input.filterHWE == true",
          numericInput("qualHWE", "Hardy-Weinberg equilibrium exact test p-value Threshold (--hwe)", value = 0.000001, min = 0.0000000001),
          numericInput("kval", "K parameter (Greer et al. 2024) to adjust p-value threshold", value = 0.001)
        ),
        checkboxInput("filterLD", "Filter Variants (--indep-pairwise)", value = FALSE),
        helpText("Prune markers in approximate linkage equilibrium with each other."),
        conditionalPanel(
          "input.filterLD == true",
          numericInput("ldWindow", "Window Size (kb)", value = 500, min = 1, step = 1),
          numericInput("ldStep", "Step Size (variants)", value = 50, min = 1, step = 1),
          numericInput("ldR2", "r2 Threshold", value = 0.2, min = 0, max = 1, step = 0.01)
        ),
        checkboxInput("cutoffKing", "Filter based on Relationships (--king-cutoff)", value = FALSE),
        helpText("Exclude a member of a pair with a kinship coefficient greater than the threshold. 
                 Use '0.354' to screen for monozygotic twins and duplicate samples, '0.177' for 1st-degree, 
                 '0.0884' for 2nd-degree, and '0.0442' for 3rd-degree relationships."),
        conditionalPanel(
          "input.cutoffKing == true",
          selectInput("kingThresh", "Kinship Coefficient", choices = c("0.354", "0.177", "0.0884", "0.0442"), selected = "0.354")
        ),
        textInput("customFilter", "Additional PLINK flags", placeholder = "--keep filestokeep.txt"),
        fileInput("extraFile1", "Optional file for first flag", accept = c(".txt", ".ped", ".psam", ".pheno", ".xlsx", ".csv")),
        fileInput("extraFile2", "Optional file for second flag", accept = c(".txt", ".ped", ".psam", ".pheno", ".xlsx", ".csv")),
        helpText("Upload extra files only if required by additional PLINK flags."),
        actionButton("calcDP", "Run Filtering & Plotting", icon = icon("filter")),
        textOutput("filterWarning")
      ),
      tabBox(
        tabPanel(
          title = "Instructions",
          h4("Filter individuals and variants using standard options in PLINK 2.0 (Chang et al., 2015)"),
          p(strong("Input file/s:"), ".vcf, .vcf.gz, .bcf, or PLINK (bed/bim/fam) files"),
          p(strong("Parameter/s:")),
          tags$ul(
            tags$li("--mind [value]"),
            tags$li("--geno [value]"),
            tags$li("--maf [value]"),
            tags$li("--qual-threshold [value]"),
            tags$li("--hwe [value]"),
            tags$li("--indep-pairwise [value]"),
            tags$li("--king-cutoff [value]"),
            tags$li("Other additional PLINK flags such as '--within [phenotype_pop.txt] --loop-cats'")
          ),
          p(strong("Expected output/s:")),
          tags$ul(
            tags$li("VCF file"),
            tags$li("Depth of Coverage Plots")
          ),
          br(),
          p(
            "Standard filtering flags are indicated. For other PLINK flags, see the following for options to be specified in the 'Additional PLINK flags' text box: ",
            tags$a("PLINK 2.0 Documentation",
              href = "https://www.cog-genomics.org/plink/2.0/",
              target = "_blank"
            )
          ),
          br(),
          p("PLINK is used strictly for filtering variants/individuals. Analysis that produces files other than a VCF file will not be shown here.")
        ),
        tabPanel(
          title = "Download sample files",
          tags$a("Sample VCF", href = "sample_hgdp.vcf", download = "sample_hgdp.vcf")
        )
      ),
      tabBox(
        title = "Results",
        tabPanel(
          title = "PLINK Commands Preview",
          verbatimTextOutput("plinkCommandPreview"),
        ),
        tabPanel(
          title = "Depth Plots",
          imageOutput("depthMarkerPlot"),
          imageOutput("depthSamplePlot")
        ),
        tabPanel(
          title = "Download Files",
          uiOutput("depthMarkerPlot_UI"),
          uiOutput("depthSamplePlot_UI"),
          uiOutput("downloadFilteredFile_UI")
        )
      )
    )
  )
}
