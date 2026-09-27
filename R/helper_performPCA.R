#' Find principal components
#'
#' @param fsnps_gen The genind data containing genotype and population information.
#' @param popinfo To indicate if population metadata is present. Default is TRUE.
#'
#' @returns A list of dataframes containing eigenvalues.
compute_pca <- function(fsnps_gen, popinfo = TRUE) {
  x <- adegenet::tab(fsnps_gen, NA.method = "mean")

  # remove SNPs with missing vals
  x <- x[, colSums(is.na(x)) == 0]

  if (ncol(x) < 2) {
    stop("Not enough variable SNPs for PCA after filtering.")
  }

  set.seed(9999)
  pca1 <- ade4::dudi.pca(x, scannf = FALSE, scale = FALSE, nf = 7)
  percent <- pca1$eig / sum(pca1$eig) * 100

  ind_coords <- as.data.frame(pca1$li)
  colnames(ind_coords) <- paste0("PC", seq_len(ncol(ind_coords)))
  ind_coords$Ind <- adegenet::indNames(fsnps_gen)

  if (isTRUE(popinfo)) {
    ind_coords$Site <- fsnps_gen@pop

    centroid <- stats::aggregate(
      ind_coords[, grep("^PC", names(ind_coords))],
      by = list(ind_coords$Site),
      FUN = mean
    )
    colnames(centroid)[1] <- "Site"
    centroid <- as.data.frame(centroid)
  } else {
    centroid <- NULL
  }

  return(list(
    pca1 = pca1,
    percent = percent,
    ind_coords = ind_coords,
    centroid = centroid
  ))
}

#' Set up labels for PCA plotting
#'
#' @param fsnps_gen The genind data containing genotype and population information.
#' @param use_default Indicates if default shapes and colors will be used. Default is TRUE.
#' @param label_file The file path to custom shapes and colors.
#'
#' @returns A list of custom labels, shapes, and colors.
get_labels <- function(fsnps_gen, use_default = TRUE, label_file = NULL, popinfo = TRUE) {
  if (use_default) {
    labels <- levels(as.factor(fsnps_gen@pop))
    n <- as.integer(length(labels))
    colors <- rep(
      RColorBrewer::brewer.pal(9, "Set1"),
      length.out = n
    )
    shapes <- rep(21:25, length.out = n)
  } else {
    if (is.null(label_file)) {
      stop("Please provide a label file.")
    }

    df <- load_csv_xlsx_files(label_file)
    if (ncol(df) < 3) {
      stop("Label file must contain at least three columns.")
    }

    labels <- trimws(as.character(df[[1]]))
    colors <- trimws(as.character(df[[2]]))
    colors <- setNames(colors, labels)
    shapes <- as.numeric(df[[3]])
    shapes <- setNames(shapes, labels)

    if (length(unique(labels)) != length(labels)) {
      stop("Duplicate population names found.")
    }

    expected_labels <- sort(unique(as.character(fsnps_gen@pop)))
    given_labels <- sort(unique(labels))
    missing <- BiocGenerics::setdiff(expected_labels, given_labels)
    extra <- BiocGenerics::setdiff(given_labels, expected_labels)

    if (length(missing) > 0 || length(extra) > 0) {
      stop("Populations in the dataset and given labels do not match.")
    }

    if (length(labels) != length(colors) ||
      length(labels) != length(shapes)
    ) {
      stop("Population, color, and shape columns must have the same length")
    }
  }
  return(list(labels = labels, colors = colors, shapes = shapes))
}


#' Plots PCA
#'
#' @param ind_coords The eigenvalues of the principal component.
#' @param centroid The the average eigenvalue per site/population.
#' @param percent The percentage of variance explained by the set x and y values.
#' @param labels_colors A list of custom shapes and colors.
#' @param width Output png's width.
#' @param height Output png's height.
#' @param pc_x Eigenvector to be plotted on the x-axis.
#' @param pc_y Eigenvector to be plotted on the y-axis.
#'
#' @returns A PCA plot.
plot_pca <- function(ind_coords,
                     centroid,
                     percent,
                     labels_colors,
                     width = 8,
                     height = 8,
                     pc_x = 1,
                     pc_y = 2,
                     highlight_pop = NULL) {
  # Ensure data frames
  if (!is.data.frame(ind_coords)) ind_coords <- as.data.frame(ind_coords)
  if (!is.data.frame(centroid)) centroid <- as.data.frame(centroid)

  has_pop <- !is.null(labels_colors) && !is.null(centroid)

  if (!is.null(highlight_pop)) {
    ind_coords$highlight <- ind_coords$Site %in% highlight_pop
  } else {
    ind_coords$highlight <- TRUE # or null?
  }

  # Axis labels
  xlab <- paste("PC", pc_x, " (", format(round(percent[pc_x], 1), nsmall = 1), "%)", sep = "")
  ylab <- paste("PC", pc_y, " (", format(round(percent[pc_y], 1), nsmall = 1), "%)", sep = "")

  # Theme
  ggtheme <- theme(
    axis.text.y = element_text(colour = "black", size = 12),
    axis.text.x = element_text(colour = "black", size = 12),
    axis.title = element_text(colour = "black", size = 12),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 15)
  )

  if (has_pop) {
    plot <- ggplot(ind_coords, aes(
      x = .data[[paste0("PC", pc_x)]],
      y = .data[[paste0("PC", pc_y)]],
      fill = Site,
      shape = Site
    )) +
      geom_hline(yintercept = 0) +
      geom_vline(xintercept = 0) +
      geom_point(aes(alpha = highlight),
        # color = "black",
        size = 3,
        show.legend = FALSE
      ) +
      ggrepel::geom_label_repel(
        data = centroid,
        aes(
          x = .data[[paste0("PC", pc_x)]],
          y = .data[[paste0("PC", pc_y)]],
          label = Site,
          fill = Site
        ),
        color = "black",
        size = 4, show.legend = FALSE
      ) +
      scale_fill_manual(values = labels_colors$colors) +
      scale_shape_manual(values = labels_colors$shapes) +
      scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.20)) +
      labs(x = xlab, y = ylab) +
      ggtheme
  } else {
    plot <- ggplot(ind_coords, aes(
      x = .data[[paste0("PC", pc_x)]],
      y = .data[[paste0("PC", pc_y)]]
    )) +
      geom_hline(yintercept = 0) +
      geom_vline(xintercept = 0) +
      geom_point(color = "black", size = 3, fill = "grey70", stroke = 0.4, shape = 21) +
      labs(x = xlab, y = ylab) +
      ggtheme
  }

  return(plot)
}
