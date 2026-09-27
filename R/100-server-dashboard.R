dashboard_server <- function(input, output, session, rv) {
  observeEvent(input$refreshApp, {
    session$reload()
  })


  zoom <- reactiveVal(1)

  observeEvent(input$smaller_workplan, {
    zoom(max(0.5, zoom() - 0.1))
  })

  observeEvent(input$bigger_workplan, {
    zoom(min(3, zoom() + 0.1))
  })

  output$workplanImg <- renderUI({
    tags$img(
      src = "readme/chart.png",
      style = paste0(
        "transform: scale(", zoom(), ");",
        "transform-origin: top center;",
        "transition: transform 0.2s;"
      )
    )
  })


  observeEvent(input$smaller_fc, {
    zoom(max(0.5, zoom() - 0.1))
  })

  observeEvent(input$bigger_fc, {
    zoom(min(3, zoom() + 0.1))
  })

  output$fileConvTable <- renderUI({
    tags$img(
      src = "fileconv.png",
      style = paste0(
        "transform: scale(", zoom(), ");",
        "transform-origin: top center;",
        "transition: transform 0.2s;"
      )
    )
  })

  output$referenceTexts <- renderUI({
    p_lists <- list(
      p(
        "Pritchard, J.K., Stephens, M., & Donnelly, P. (2000). Inference of Population Structure Using Multilocus Genotype Data. Genetics Society of America, 155, 945-959.",
        tags$a("https://doi.org/10.1093/genetics/155.2.945",
          href = "https://doi.org/10.1093/genetics/155.2.945",
          target = "_blank"
        )
      ),
      p(
        "Falush, D., Stephens, M., & Pritchard, J.K. (2003). Inference of Population Structure Using Multilocus Genotype Data: Linked Loci and Correlated Allele Frequencies. Genetics Society of America, 164, 1567-1587.",
        tags$a("https://doi.org/10.1093/genetics/164.4.1567",
          href = "https://doi.org/10.1093/genetics/164.4.1567",
          target = "_blank"
        )
      ),
      p(
        "Falush, D., Stephens, M., & Pritchard, J.K. (2007). Inference of population structure using multilocus genotype data: dominant markers and null alleles. Molecular Ecology Notes, 7(4), 574-578.",
        tags$a("https://doi.org/10.1111/j.1471-8286.2007.01758.x",
          href = "https://doi.org/10.1111/j.1471-8286.2007.01758.x",
          target = "_blank"
        )
      ),
      p(
        "Hubisz, M. J., Falush, D., Stephens, M., & Pritchard, J. K. (2009). Inferring weak population structure with the assistance of sample group information. Molecular ecology resources, 9(5), 1322–1332.",
        tags$a("https://doi.org/10.1111/j.1755-0998.2009.02591.x",
          href = "https://doi.org/10.1111/j.1755-0998.2009.02591.x",
          target = "_blank"
        )
      ),
      p(
        "Purcell, S., Neale, B., Todd-Brown, K., Thomas, L., Ferreira, M. A., Bender, D., Maller, J., Sklar, P., de Bakker, P. I., Daly, M. J., & Sham, P. C. (2007). PLINK: a tool set for whole-genome association and population-based linkage analyses. American journal of human genetics, 81(3), 559–575.",
        tags$a("https://doi.org/10.1086/519795",
          href = "https://doi.org/10.1086/519795",
          target = "_blank"
        )
      ),
      p(
        "Chang, C.C., Chow, C.C., Tellier, L.C., Vattikuti, S., Purcell, S.M., & Lee, J.J. (2015). Second-generation PLINK: rising to the challenge of larger and richer datasets. GigaScience, 4(7).",
        tags$a("https://doi.org/10.1186/s13742-015-0047-8",
          href = "https://doi.org/10.1186/s13742-015-0047-8",
          target = "_blank"
        )
      ),
      p(
        "Dray, S., & Dufour, A.-B. (2007). The ade4 Package: Implementing the Duality Diagram for Ecologists. Journal of Statistical Software, 22(4), 1–20.",
        tags$a("https://doi.org/10.18637/jss.v022.i04",
          href = "https://doi.org/10.18637/jss.v022.i04",
          target = "_blank"
        )
      ),
      p(
        "Kamvar, Z.N., Tabima, J.F., & Grunwald, N.J. (2014). Poppr: an R package for genetic analysis of populations with clonal, partially clonal, and/or sexual reproduction. PeerJ (2), e281,",
        tags$a("https://doi.org/10.7717/peerj.281",
          href = "https://doi.org/10.7717/peerj.281",
          target = "_blank"
        )
      ),
      p(
        "Goudet, J. (2004). hierfstat, a package for r to compute and test hierarchical F-statistics. Molecular Ecology Notes, 5(1), 184-186.",
        tags$a("https://doi.org/10.1111/j.1471-8286.2004.00828.x",
          href = "https://doi.org/10.1111/j.1471-8286.2004.00828.x",
          target = "_blank"
        )
      ),
      p(
        "Jombart, T. (2008). adegenet: a R package for the multivariate analysis of genetic markers. Bioinformatics, 24(11), 1403-1405.",
        tags$a("https://doi.org/10.1093/bioinformatics/btn129",
          href = "https://doi.org/10.1093/bioinformatics/btn129",
          target = "_blank"
        )
      ),
      p(
        "Jombart, T., & Ahmed, I. (2011). adegenet 1.3-1: new tools for the analysis of genome-wide SNP data. Bioinformatics (Oxford, England), 27(21), 3070–3071.",
        tags$a("https://doi.org/10.1093/bioinformatics/btr521",
          href = "https://doi.org/10.1093/bioinformatics/btr521",
          target = "_blank"
        )
      ),
      p(
        "Paradis E. (2010). pegas: an R package for population genetics with an integrated-modular approach. Bioinformatics (Oxford, England), 26(3), 419–420.",
        tags$a("https://doi.org/10.1093/bioinformatics/btp696",
          href = "https://doi.org/10.1093/bioinformatics/btp696",
          target = "_blank"
        )
      ),
      p("Petit, R.J., El Mousadik, A., & Pons, O. (1998). Identifying populations for conservation on the basis of genetic markers. Conservation Biology, 12:844-855"),
      p(
        "Foulley, J.F., & Ollivier, L. (2005). Estimating allelic richness and its diversity. Livestock Science, 101:150-158.",
        tags$a("https://doi.org/10.1016/j.livprodsci.2005.10.021",
          href = "https://doi.org/10.1016/j.livprodsci.2005.10.021",
          target = "_blank"
        )
      ),
      p(
        "Nei, M. (1978). Estimation of average heterozygosity and genetic distance from a small number of individuals. Genetics:89:583-590.",
        tags$a("https://doi.org/10.1093/genetics/89.3.583",
          href = "https://doi.org/10.1093/genetics/89.3.583",
          target = "_blank"
        )
      ),
      p("Rousset, F. (2002). Inbreeding and relatedness coefficients: what do they measure? Heredity, 88:371-380."),
      p(
        "Rezaei, N., & Hedayat, M. (2013). Allele Frequency. Brenner's Encyclopedia of Genetics (Second Edition).",
        tags$a("https://doi.org/10.1016/B978-0-12-374984-0.00032-2",
          href = "https://doi.org/10.1016/B978-0-12-374984-0.00032-2",
          target = "_blank"
        )
      ),
      p("Tiret, L., & Cambien, F. (1995). Departure from Hardy-Weinberg equilibrium should be systematically tested in studies of association between genetic markers and disease. Circulation, 92(11):3364-3365."),
      p(
        "Weir, B.S., & Cockerham, C.C. (1984). Estimating F-statistics for the analysis of population structure. Evolution; International Journal of Organic Evolution, 38(6): 1358-1370.",
        tags$a("https://doi.org/10.1111/j.1558-5646.1984.tb05657.x",
          href = "https://doi.org/10.1111/j.1558-5646.1984.tb05657.x",
          target = "_blank"
        )
      ),
      p(
        "Archer, F. (2025). strataG: An R package for manipulating, summarizing and analysing population genetic data. R package version 1.0.6. Zenodo.",
        tags$a("http://doi.org/10.5281/zenodo.60416",
          href = "http://doi.org/10.5281/zenodo.60416",
          target = "_blank"
        )
      ),
      p(
        "Jakobsson, M., Rosenberg, N.A (2007). CLUMPP: a cluster matching and permutation program for dealing with label switching and multimodality in analysis of population structure. Bioinformatics, 23(14), 1801–1806",
        tags$a("https://doi.org/10.1093/bioinformatics/btm233",
          href = "https://doi.org/10.1093/bioinformatics/btm233",
          target = "_blank"
        )
      ),
      p(
        "Gouy, A., & Zieger, M. (n.d.). The STRAF Book: Statistical Forensics made easy.",
        tags$a("https://agouy.github.io/straf_book/index.html",
          href = "https://agouy.github.io/straf_book/index.html",
          target = "_blank"
        )
      ),
      p(
        "Dimitriadou, E., Hornik, K., Leisch, F., Meyer, D., & Weingessel, A. (2009). E1071: Misc Functions of the Department of Statistics (E1071), TU Wien.",
        tags$a("https://www.researchgate.net/publication/221678005_E1071_Misc_Functions_of_the_Department_of_Statistics_E1071_TU_Wien",
          href = "https://www.researchgate.net/publication/221678005_E1071_Misc_Functions_of_the_Department_of_Statistics_E1071_TU_Wien",
          target = "_blank"
        )
      ),
      p(
        "Kuhn, M. (2008). Building Predictive Models in R Using the caret Package. Journal of Statistical Software, 28(5), 1–26.",
        tags$a("https://doi.org/10.18637/jss.v028.i05",
          href = "https://doi.org/10.18637/jss.v028.i05",
          target = "_blank"
        )
      ),
      p(
        "Cheung, E. Y. Y., Gahan, M. E., & McNevin, D. (2017). Prediction of biogeographical ancestry from genotype: a comparison of classifiers. International journal of legal medicine, 131(4), 901–912.",
        tags$a("https://doi.org/10.1007/s00414-016-1504-3",
          href = "https://doi.org/10.1007/s00414-016-1504-3",
          target = "_blank"
        )
      ),
      p(
        "Bodenhofer, U., Bonatesta, E., Horejš-Kainrath, C., & Hochreiter, S. (2015). msa: an R package for multiple sequence alignment. Bioinformatics (Oxford, England), 31(24), 3997–3999.",
        tags$a("https://doi.org/10.1093/bioinformatics/btv494",
          href = "https://doi.org/10.1093/bioinformatics/btv494",
          target = "_blank"
        )
      ),
      p(
        "Zou, Y., Zhang, Z., Zeng, Y., Hu, H., Hao, Y., Huang, S., & Li, B. (2024). Common Methods for Phylogenetic Tree Construction and Their Implementation in R. Bioengineering (Basel, Switzerland), 11(5), 480.",
        tags$a("https://doi.org/10.3390/bioengineering11050480",
          href = "https://doi.org/10.3390/bioengineering11050480",
          target = "_blank"
        )
      ),
      p(
        "Wright E. S. (2015). DECIPHER: harnessing local sequence context to improve protein multiple sequence alignment. BMC bioinformatics, 16, 322.",
        tags$a("https://doi.org/10.1186/s12859-015-0749-z",
          href = "https://doi.org/10.1186/s12859-015-0749-z",
          target = "_blank"
        )
      ),
      p(
        "Paradis, E., & Schliep, K. (2019). ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. Bioinformatics (Oxford, England), 35(3), 526–528.",
        tags$a("https://doi.org/10.1093/bioinformatics/bty633",
          href = "https://doi.org/10.1093/bioinformatics/bty633",
          target = "_blank"
        )
      ),
      p(
        "Schliep, K.P. (2011). phangorn: phylogenetic analysis in R. Bioinformatics, 27(4), 592-593.",
        tags$a("https://doi.org/10.1093/bioinformatics/btq706",
          href = "https://doi.org/10.1093/bioinformatics/btq706",
          target = "_blank"
        )
      ),
      p(
        "Zhang, A., Hao, M., Yang, C., & Shi, Z. (2016). BarcodingR: an integrated r package for species identification using DNA barcodes. Methods in Ecology and Evolution, 8, 627-634.",
        tags$a("https://doi.org/10.1111/2041-210X.12682",
          href = "https://doi.org/10.1111/2041-210X.12682",
          target = "_blank"
        )
      ),
      p(
        "Behr, A.A., Liu, K.Z., Liu-Fang, G., Nakka, P., & Ramachandran, S. (2016). pong: fast analysis and visualization of latent clusters in population genetic data. Bioinformatics, 32(18), 2817-2823.",
        tags$a("https://doi.org/10.1093/bioinformatics/btw327",
          href = "https://doi.org/10.1093/bioinformatics/btw327",
          target = "_blank"
        )
      ),
      p(
        "Gruber, B., Unmack, P. J., Berry, O. F., & Georges, A. (2018). dartr: An r package to facilitate analysis of SNP data generated from reduced representation genome sequencing. Molecular ecology resources, 18(3), 691–699.",
        tags$a("https://doi.org/10.1111/1755-0998.12745",
          href = "https://doi.org/10.1111/1755-0998.12745",
          target = "_blank"
        )
      ),
      p(
        "Archer, F. I., Adams, P. E., & Schneiders, B. B. (2017). stratag: An r package for manipulating, summarizing and analysing population genetic data. Molecular ecology resources, 17(1), 5–11.",
        tags$a("https://doi.org/10.1111/1755-0998.12559",
          href = "https://doi.org/10.1111/1755-0998.12559",
          target = "_blank"
        )
      ),
      p(
        "Barash, M., McNevin, D., Fedorenko, V., & Giverts, P. (2024). Machine learning applications in forensic DNA profiling: A critical review. Forensic Science International: Genetics (69).",
        tags$a("https://doi.org/10.1016/j.fsigen.2023.102994",
          href = "https://doi.org/10.1016/j.fsigen.2023.102994",
          target = "_blank"
        )
      )
    )

    do.call(tagList, p_lists)
  })

  output$popstatRef <- renderUI({
    p_lists <- list(
      h4("Calculate common population statistics using R-based packages"),
      p(strong("Input file:"), "CSV or Excel file"),
      p(strong("Expected output files:")),
      tags$ul(
        tags$li("Excel (.xlsx) file with all results"),
        tags$li("Heterozygosity Plot (.png)"),
        tags$li("Fst Plots (.png)")
      ),
      p("Population statistics included:"),
      tags$ul(
        tags$li(strong("Private alleles"), "[1] calculated using the poppr R package (Kamvar et al., 2014)"),
        tags$li(strong("Mean Allelic Richness"), "[2] using the hierfstat R package (Goudet, 2004)"),
        tags$li(strong("Heterozygosity"), "[3] using the hierfstat R package (Goudet, 2004). Exact test based on Monte Carlo permutations (Guo & Thompson, 1992) is computed using 1000 replicates."),
        tags$li(strong("Inbreeding Coefficient"), "[4] using the hierfstat R package (Goudet, 2004)"),
        tags$li(strong("Allele frequency"), "[5] using the adegenet R package (Jombart, 2008)"),
        tags$li(strong("Hardy-Weinberg equilibrium"), "[6] using the pegas R package (Paradis, 2010)"),
        tags$li(strong("FST values"), "[7] using the hierfstat R package (Goudet, 2004)")
      ),
      br(),
      h5("To learn more about the statistics:"),
      h5("References"),
      p("[1] Petit, R.J., El Mousadik, A., and Pons, O. (1998). Identifying populations for conservation on the basis of genetic markers. Conservation Biology, 12:844-855"),
      p("[2] Foulley, J.F., and Ollivier, L. (2005). Estimating allelic richness and its diversity. Livestock Science, 101:150-158. https://doi.org/10.1016/j.livprodsci.2005.10.021"),
      p("[3] Nei, M. (1978). Estimation of average heterozygosity and genetic distance from a small number of individuals. Genetics:89:583-590. https://doi.org/10.1093/genetics/89.3.583"),
      p("[4] Rousset, F. (2002). Inbreeding and relatedness coefficients: what do they measure? Heredity, 88:371-380."),
      p("[5] Rezaei, N., and Hedayat, M. (2013). Allele Frequency. Brenner's Encyclopedia of Genetics (Second Edition). https://doi.org/10.1016/B978-0-12-374984-0.00032-2"),
      p("[6] Tiret, L., and Cambien, F. (1995). Departure from Hardy-Weinberg equilibrium should be systematically tested in studies of association between genetic markers and disease. Circulation, 92(11):3364-3365."),
      p("[7] Weir, B.S., and Cockerham, C.C. (1984). Estimating F-statistics for the analysis of population structure. Evolution; International Journal of Organic Evolution, 38(6): 1358-1370. https://doi.org/10.1111/j.1558-5646.1984.tb05657.x")
    )

    do.call(tagList, p_lists)
  })
}
