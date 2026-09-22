# ======================================
# R function to export a genind object in STRUCTURE format
#
# Tom Jenkins t.l.jenkins@exeter.ac.uk
#
# July 2018
#
# ======================================

# data: genind object
# file: file name to write
# pops: whether to include population info
# markers: whether to include marker info
# unix: export as a unix text file (windows text file default)
# Function is flexible with regards to ploidy, although genotypes are
# considered to be unambiguous.
# Missing data must be recorded as NA.
# SNP data must be biallelic with all alleles present.

# Example use:
# library(adegenet)
# ind = as.character(paste("ind_", seq(1:100), sep=""))
# pop = as.character(c(rep("pop1",25), rep("pop2",25), rep("pop3",25), rep("pop4",25)))
# loci = list(c("AA","AC","CC"), c("GG","GC","CC"), c("TT","TA","AA"), c("CC","CT","TT"))
# loci = sample(loci, 100, replace=T)
# loci = lapply(loci, sample, size=100, replace=TRUE)
# geno = as.data.frame(loci, col.names= .genlab("loc",100))
# data = df2genind(geno, ploidy=2, ind.names=ind, pop=pop, sep="")
# genind2structure(data, file="example_structure.str")

# Convert Windows text file to Unix text file in Linux
# awk '{ sub("\r$", ""); print }' winfile.txt > unixfile.txt

genind2structure2 <- function(data, file = "", pops = TRUE, markers = TRUE, unix = FALSE) {
  ## Check input file a genind object
  if (!"genind" %in% class(data)) {
    warning("Function was designed for genind objects.")
  }
  
  ## Check adegenet, miscTools and stringr are installed
  if (!require(adegenet)) {
    install.packages("adegenet")
  }
  if (!require(miscTools)) {
    install.packages("miscTools")
  }
  if (!require(stringr)) {
    install.packages("stringr")
  }
  
  
  # ---------------- #
  #
  # Preamble
  #
  # ---------------- #
  
  ## Ploidy
  ploid <- max(data$ploidy)
  
  ## Number of individuals
  ind <- nInd(data)
  
  ## Create dataframe containing individual labels
  ## Number of duplicated labels depends on the ploidy of the dataset
  df <- data.frame(ind = rep(indNames(data), each = ploid))
  
  ## Locus names
  loci <- locNames(data)
  
  
  # ---------------- #
  #
  # Population IDs
  #
  # ---------------- #
  
  if (pops) {
    ## Create dataframe containing individual labels
    ## Number of duplicated labels depends on the ploidy of the dataset
    df.pop <- data.frame(pop = rep(as.numeric(data$pop), each = ploid))
    
    ## Add population IDs to dataframe
    df$Pop <- df.pop$pop
  }
  
  # ---------------- #
  #
  # Process genotypes
  #
  # ---------------- #
  
  ## Add columns for genotypes
  df <- cbind(df, matrix(-9, # -9 codes for missing data
                         nrow = dim(df)[1],
                         ncol = nLoc(data),
                         dimnames = list(NULL, loci)
  ))
  
  ## Loop through dataset to extract genotypes
  for (L in loci) {
    #=========================== revised by forcing 'thesedata' as a matrix in restful forensics 
    thesedata <- as.matrix(data$tab[, grep(paste("^", L, "\\.", sep = ""), dimnames(data$tab)[[2]])])# dataotypes by locus
    al <- 1:dim(thesedata)[2] # numbered alleles
    for (s in 1:ind) {
      if (all(!is.na(thesedata[s, ]))) {
        tabrows <- (1:dim(df)[1])[df[[1]] == indNames(data)[s]] # index of rows in output to write to
        tabrows <- tabrows[1:sum(thesedata[s, ])] # subset if this is lower ploidy than max ploidy
        df[tabrows, L] <- rep(al, times = thesedata[s, ])
      }
    }
  }
  
  
  # ---------------- #
  #
  # Marker IDs
  #
  # ---------------- #
  
  if (markers) {
    ## Add a row at the top containing loci names
    df <- as.data.frame(insertRow(as.matrix(df), 1, c(loci, "", "")))
    df <- df[-1,]
    df <- rbind(colnames(df), df)
  }
  
  # ---------------- #
  #
  # Export file
  #
  # ---------------- #
  
  ## Export dataframe
  write.table(df,
              file = file, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE
  )
  
  ## Export dataframe as unix text file
  if (unix) {
    output_file <- file(file, open = "wb")
    write.table(df,
                file = output_file, sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = FALSE
    )
    close(output_file)
  }
}


#' Generate structure R object
#' 
#' @param files The directory containing structure results.
#' @param pops The populations in the input file.
#'
#' @returns The structure object generated from actual STRUCTURE resulting files.
structureImport <- function(files, pops = NULL) {
  if (length(files) == 0) {
    stop("No STRUCTURE files found.")
  }
  
  if (is.null(pops)) {
    pops <- structure_get_populations(files[1])
  }
  
  sr <- lapply(files, function(f) {
    result <- strataG::structureRead(f, pops = pops)
    
    result$files <- f
    result$label <- basename(f)
    
    result
  })
  
  names(sr) <- vapply(sr, function(x) x$label, character(1))
  sr <- sr[order(
    as.numeric(sub(".*k([0-9]+).*", "\\1", names(sr))),
    as.numeric(sub(".*r([0-9]+).*", "\\1", names(sr)))
  )]
  class(sr) <- c("structure.result", "list")
  
  return(sr)
}

#' Extract populations from STRUCTURE files
#' 
#' @param file The directory containing structure results.
#'
#' @returns A list of populations.
structure_get_populations <- function(file){
  txt <- readLines(file)
  
  first <- grep("(%Miss)", txt, fixed = TRUE) + 1
  last <- grep("Estimated Allele", txt, fixed = TRUE) - 1
  
  if (length(first) == 0 || length(last) == 0) {
    stop("Cannot find STRUCTURE Q matrix in _f files.")
  }
  
  q.txt <- txt[first:last]
  
  # same as structureRead()
  q.txt <- sub("[\\*]+", "", q.txt)
  q.txt <- sub("[(]", "", q.txt)
  q.txt <- sub("[)]", "", q.txt)
  q.txt <- sub("[|][ ]+$", "", q.txt)
  
  pops <- vapply(strsplit(q.txt, "\\s+"), function(x){
    x <- x[x != ""]
    x[4]
  }, character(1))
  
  return(unique(pops))
}