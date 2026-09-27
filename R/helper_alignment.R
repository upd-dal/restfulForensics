#' Perform multiple sequence alignment
#'
#' @param files The DNA bin object.
#' @param algorithm The method for alignment (ClustalW, ClustalOmega, MUSCLE).
#'
#' @returns A list of aligned sequences.
calc_msa <- function(files, algorithm) {
  # Creating Substitution Matrix
  personal_matrix <- pwalign::nucleotideSubstitutionMatrix(
    match = 1, mismatch = 0, baseOnly = FALSE, type = "DNA"
  )

  gap_penalty <- -2
  personal_matrix <- rbind(personal_matrix, "-" = gap_penalty)
  personal_matrix <- cbind(personal_matrix, "-" = gap_penalty)
  colnames(personal_matrix) <- c("A", "C", "G", "T", "M", "R", "W", "S", "Y", "K", "V", "H", "D", "B", "N", "-")
  rownames(personal_matrix) <- c("A", "C", "G", "T", "M", "R", "W", "S", "Y", "K", "V", "H", "D", "B", "N", "-")
  personal_matrix <- as.matrix(personal_matrix)

  # perform msa
  aligned_sequences <- msa::msa(
    files,
    substitutionMatrix = personal_matrix,
    method = algorithm
  )

  # calculate alignment score
  alignment_scores <- msa::msaConservationScore(
    aligned_sequences,
    substitutionMatrix = personal_matrix
  )

  # Post-processing
  aligned_dnastrings <- msa::msaConvert(
    aligned_sequences,
    type = "seqinr::alignment"
  )
  aligned_dnastrings <- Biostrings::DNAStringSet(
    setNames(aligned_dnastrings$seq, aligned_dnastrings$nam)
  )

  adjusted <- DECIPHER::AdjustAlignment(aligned_dnastrings)
  staggered <- DECIPHER::StaggerAlignment(adjusted)

  return(list(
    alignment = aligned_sequences,
    scores = alignment_scores,
    adjusted = adjusted,
    staggered = staggered
  ))
}


#' Build phylogenetic tree using NJ method
#'
#' @param alignment The alignment object.
#' @param outgroup The sample name to serve as the outgroup.
#' @param seed The value for bootstrapping.
#' @param model The substitution model.
#'
#' @returns A phylogenetic tree plot.
build_nj_tree <- function(alignment, outgroup = NULL, seed = 123, model = model) {
  bins <- ape::as.DNAbin(alignment)
  distance <- ape::dist.dna(bins, model = model)
  nj_tree <- ape::nj(distance)

  if (!is.null(outgroup) && outgroup %in% nj_tree$tip.label) {
    nj_tree <- ape::root(nj_tree, outgroup = outgroup)
  }

  nj_tree <- ape::ladderize(nj_tree)
  num_sites <- ncol(bins)
  if (num_sites < 10) {
    warning("Alignment has fewer than 10 sites. Skipping bootstrap.")
    tree_plot <- ggtree(nj_tree, branch.length = "none") +
      theme_tree2() +
      geom_tiplab() +
      ggtitle("NJ Tree")
    return(tree_plot)
  }

  # Bootstrapping
  set.seed(seed)
  boots <- ape::boot.phylo(nj_tree, bins,
    FUN = function(x) {
      tree <- ape::nj(ape::dist.dna(x, model = model))

      if (!is.null(outgroup) && outgroup %in% tree$tip.label) {
        tree <- ape::root(tree, outgroup = outgroup)
      }

      ape::ladderize(tree)
    }, rooted = TRUE
  )

  boots[is.na(boots)] <- 0
  nj_tree$node.label <- as.character(boots)
  tree_plot <- ggtree(nj_tree, branch.length = "none") +
    theme_tree2() +
    geom_tiplab() +
    geom_text2(aes(subset = !isTip, label = label), hjust = -0.3) +
    ggtitle("NJ Tree") +
    xlim(0, 20)

  return(tree_plot)
}


#' Build phylogenetic tree using UPGMA method
#'
#' @inheritParams build_nj_tree
#'
#' @returns A phylogenetic tree plot.
build_upgma_tree <- function(alignment, outgroup = NULL, seed = 123, model = model) {
  bins <- ape::as.DNAbin(alignment)
  distance <- ape::dist.dna(bins, model = model)
  upgma_tree <- phangorn::upgma(distance)

  if (!is.null(outgroup) && outgroup %in% upgma_tree$tip.label) {
    upgma_tree <- ape::root(upgma_tree, outgroup = outgroup)
  }

  upgma_tree <- ape::ladderize(upgma_tree)
  num_sites <- ncol(bins)
  if (num_sites < 10) {
    warning("Alignment has fewer than 10 sites. Skipping bootstrap.")

    tree_plot <- ggtree(upgma_tree, branch.length = "none") +
      theme_tree2() +
      geom_tiplab() +
      ggtitle("UPGMA Tree")

    return(tree_plot)
  }

  set.seed(seed)
  boots <- ape::boot.phylo(upgma_tree, bins,
    FUN = function(x) {
      tree <- phangorn::upgma(ape::dist.dna(x, model = model))

      if (!is.null(outgroup) && outgroup %in% tree$tip.label) {
        tree <- ape::root(tree, outgroup = outgroup)
      }

      ape::ladderize(tree)
    }, rooted = TRUE,
  )

  boots[is.na(boots)] <- 0
  upgma_tree$node.label <- as.character(boots)
  tree_plot <- ggtree(upgma_tree, branch.length = "none") +
    theme_tree2() +
    geom_tiplab() +
    geom_text2(aes(subset = !isTip, label = label), hjust = -0.3) +
    ggtitle("UPGMA Tree") +
    xlim(0, 20)

  return(tree_plot)
}

#' Build phylogenetic tree using Maximum Parismony method
#'
#' @param alignment The alignment object.
#' @param outgroup The sample name to serve as the outgroup.
#' @param seed The value for bootstrapping.
#' @param directory The directory to save the png plot.
#'
#' @returns A phylogenetic tree plot.
build_max_parsimony <- function(alignment, outgroup = NULL, seed = 123, directory = ".") {
  bins <- ape::as.DNAbin(alignment)
  phy <- phangorn::phyDat(bins, type = "DNA")
  dm <- phangorn::dist.ml(phy)
  start_tree <- NJ(dm)
  parsimony_tree <- phangorn::optim.parsimony(start_tree, phy)

  # Rooting
  if (!is.null(outgroup) && outgroup %in% parsimony_tree$tip.label) {
    parsimony_tree <- root(parsimony_tree, outgroup = outgroup, resolve.root = TRUE)
  }

  # boostrapping
  set.seed(seed)
  bs_pars <- phangorn::bootstrap.phyDat(phy, \(x) phangorn::optim.parsimony(NJ(phangorn::dist.ml(x)), x))

  # plot
  filename <- paste(directory, "parsimony_tree.png")
  png(filename, width = 800, height = 600)
  plotBS(parsimony_tree, bs_pars, main = "Parsimony Tree")
  dev.off()

  return(filename)
}


#' Build phylogenetic tree using Maximum Likelihood method
#'
#' @param alignment The alignment object.
#' @param outgroup The sample name to serve as the outgroup.
#' @param seed The value for bootstrapping.
#' @param bs_reps The number of boostrap replicates to run.
#' @param directory The directory to save the png plot.
#'
#' @returns A list specifying the (1) model and (2) phylogenetic tree plot.
build_ml_tree <- function(alignment,
                          outgroup = NULL,
                          seed = 123,
                          bs_reps = 100,
                          directory = ".") {
  bins <- ape::as.DNAbin(alignment)
  phy <- phyDat(bins, type = "DNA")

  dm <- dist.ml(phy)
  start_tree <- NJ(dm)

  fit <- phangorn::pml(start_tree, data = phy)

  # find best-fit model
  model_test <- phangorn::modelTest(phy, tree = start_tree)
  best_model <- model_test$Model[which.min(model_test$BIC)]
  best_model <- sub("\\+.*", "", best_model)
  fit_opt <- phangorn::optim.pml(fit, model = best_model, optGamma = TRUE, optInv = TRUE, rearrangement = "stochastic")

  tree <- fit_opt$tree
  if (!is.null(outgroup) && outgroup %in% tree$tip.label) {
    tree <- root(tree, outgroup = outgroup, resolve.root = TRUE)
  }

  # bootstrapping
  set.seed(seed)
  bs <- phangorn::bootstrap.pml(fit_opt, bs = bs_reps, optNni = TRUE)

  filename <- paste(directory, "ml_tree.png")
  png(filename, width = 800, height = 600)
  plotBS(tree, bs, main = paste("ML Tree (", best_model, ")"))
  dev.off()

  return(list(
    best_model = best_model,
    filename = filename
  ))
}
