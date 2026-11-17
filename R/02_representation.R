# Transcriptomic representation -----------------------------------------------

select_variable_genes <- function(input, n_genes = 10000) {
  validate_tumor2line_input(input)
  if (nrow(input$gene_stats) == 0) {
    stop("No genes with finite, non-zero expression variability are available for HVG selection.", call. = FALSE)
  }

  n_genes <- as.integer(n_genes)
  if (length(n_genes) != 1 || is.na(n_genes) || n_genes < 2) {
    stop("n_variable_genes must be an integer >= 2.", call. = FALSE)
  }

  n_genes <- min(n_genes, nrow(input$gene_stats))
  utils::head(input$gene_stats$gene, n_genes)
}

fit_pca_representation <- function(
    expression_log,
    selected_genes,
    scale_genes = TRUE) {

  selected_genes <- intersect(selected_genes, colnames(expression_log))
  if (length(selected_genes) < 2) {
    stop("At least two selected genes are required for PCA.", call. = FALSE)
  }
  if (nrow(expression_log) < 2) {
    stop("At least two samples are required for PCA.", call. = FALSE)
  }

  x <- expression_log[, selected_genes, drop = FALSE]
  if (scale_genes) {
    gene_center <- colMeans(x)
    gene_scale <- apply(x, 2, stats::sd)
    gene_scale[!is.finite(gene_scale) | gene_scale <= 0] <- 1
    x_scaled <- sweep(x, 2, gene_center, "-")
    x_scaled <- sweep(x_scaled, 2, gene_scale, "/")
  } else {
    gene_center <- stats::setNames(rep(0, ncol(x)), colnames(x))
    gene_scale <- stats::setNames(rep(1, ncol(x)), colnames(x))
    x_scaled <- x
  }

  x_scaled[!is.finite(x_scaled)] <- 0
  max_available_pcs <- min(if (scale_genes) nrow(x_scaled) - 1L else nrow(x_scaled), ncol(x_scaled))
  if (max_available_pcs < 1) {
    stop("Too few samples or features are available for PCA.", call. = FALSE)
  }

  # Fit all mathematically available components. The retained dimension is
  # selected afterwards from the requested cumulative-variance threshold.
  fit <- stats::prcomp(
    x_scaled,
    center = FALSE,
    scale. = FALSE,
    rank. = max_available_pcs
  )

  scores_raw <- fit$x
  pc_center <- colMeans(scores_raw)
  pc_scale <- apply(scores_raw, 2, stats::sd)
  pc_scale[!is.finite(pc_scale) | pc_scale <= 0] <- 1
  scores_scaled <- sweep(scores_raw, 2, pc_center, "-")
  scores_scaled <- sweep(scores_scaled, 2, pc_scale, "/")
  scores_scaled[!is.finite(scores_scaled)] <- 0

  explained_variance <- fit$sdev^2 / sum(fit$sdev^2)
  variance <- data.frame(
    PC_index = seq_along(explained_variance),
    explained_variance = explained_variance,
    cumulative_variance = cumsum(explained_variance),
    stringsAsFactors = FALSE
  )

  list(
    fit = fit,
    variance = variance,
    scores_raw = scores_raw,
    scores_scaled = scores_scaled,
    selected_genes = selected_genes,
    gene_center = gene_center,
    gene_scale = gene_scale,
    pc_center = pc_center,
    pc_scale = pc_scale,
    max_available_pcs = max_available_pcs,
    scale_genes = scale_genes
  )
}

select_pc_by_variance <- function(pca_res, threshold = 0.90) {
  threshold <- as.numeric(threshold)
  if (length(threshold) != 1 || !is.finite(threshold) || threshold <= 0 || threshold > 1) {
    stop("pca_variance_threshold must be in (0, 1].", call. = FALSE)
  }

  hit <- which(pca_res$variance$cumulative_variance >= threshold)
  if (length(hit) == 0) {
    selected <- nrow(pca_res$variance)
    warning(
      sprintf(
        "The requested cumulative variance threshold (%.1f%%) was not reached; using all %d available PCs (%.1f%%).",
        threshold * 100,
        selected,
        utils::tail(pca_res$variance$cumulative_variance, 1) * 100
      ),
      call. = FALSE
    )
    return(selected)
  }

  pca_res$variance$PC_index[min(hit)]
}

project_expression_to_reference <- function(expression_log, reference_representation) {
  required <- c(
    "selected_genes", "gene_center", "gene_scale", "pca_rotation",
    "pc_center", "pc_scale", "n_pcs"
  )
  missing_fields <- setdiff(required, names(reference_representation))
  if (length(missing_fields) > 0) {
    stop(
      "The reference representation is missing projection parameters: ",
      paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }

  selected_genes <- reference_representation$selected_genes
  missing_genes <- setdiff(selected_genes, colnames(expression_log))
  if (length(missing_genes) > 0) {
    stop(
      "The new expression data are missing reference genes: ",
      paste(utils::head(missing_genes, 10), collapse = ", "),
      if (length(missing_genes) > 10) " ..." else "",
      call. = FALSE
    )
  }

  x <- expression_log[, selected_genes, drop = FALSE]
  x_scaled <- sweep(x, 2, reference_representation$gene_center[selected_genes], "-")
  x_scaled <- sweep(x_scaled, 2, reference_representation$gene_scale[selected_genes], "/")
  x_scaled[!is.finite(x_scaled)] <- 0

  rotation <- reference_representation$pca_rotation[
    selected_genes,
    seq_len(reference_representation$n_pcs),
    drop = FALSE
  ]
  scores_raw <- x_scaled %*% rotation
  scores_scaled <- sweep(
    scores_raw,
    2,
    reference_representation$pc_center[seq_len(reference_representation$n_pcs)],
    "-"
  )
  scores_scaled <- sweep(
    scores_scaled,
    2,
    reference_representation$pc_scale[seq_len(reference_representation$n_pcs)],
    "/"
  )
  scores_scaled[!is.finite(scores_scaled)] <- 0

  colnames(scores_raw) <- colnames(rotation)
  colnames(scores_scaled) <- colnames(rotation)

  list(scores = scores_scaled, scores_raw = scores_raw)
}

#' Build a cohort-specific scaled-PCA representation
#'
#' Selects variable genes, performs joint gene-wise scaling across query and
#' candidate cell-line samples, fits PCA using all available components, selects
#' the smallest number of PCs reaching the requested cumulative variance
#' threshold, and scales the retained PC scores. All transformation parameters
#' required for projection of future query samples are retained.
#'
#' @param input A `Tumor2LineInput` object from [prepare_tumor2line_data()].
#' @param n_variable_genes Maximum number of variable genes to use.
#' @param pca_variance_threshold Cumulative explained-variance threshold used
#'   to select the number of retained PCs.
#' @param scale_genes Whether to center and scale genes before PCA.
#'
#' @return A list containing scaled PC scores, PCA loadings and all gene- and
#'   PC-scaling parameters required for projection.
#' @export
build_scaled_pca_representation <- function(
    input,
    n_variable_genes = 10000,
    pca_variance_threshold = 0.90,
    scale_genes = TRUE) {

  validate_tumor2line_input(input)
  selected_genes <- select_variable_genes(input, n_variable_genes)
  pca_res <- fit_pca_representation(
    expression_log = input$expression_log,
    selected_genes = selected_genes,
    scale_genes = scale_genes
  )
  n_pcs <- select_pc_by_variance(pca_res, pca_variance_threshold)
  selected_pc_names <- colnames(pca_res$scores_raw)[seq_len(n_pcs)]

  list(
    scores = pca_res$scores_scaled[, seq_len(n_pcs), drop = FALSE],
    scores_raw = pca_res$scores_raw[, seq_len(n_pcs), drop = FALSE],
    selected_genes = selected_genes,
    n_pcs = n_pcs,
    variance = pca_res$variance,
    pca_rotation = pca_res$fit$rotation[, selected_pc_names, drop = FALSE],
    gene_center = pca_res$gene_center,
    gene_scale = pca_res$gene_scale,
    pc_center = pca_res$pc_center[seq_len(n_pcs)],
    pc_scale = pca_res$pc_scale[seq_len(n_pcs)],
    preprocessing = list(
      transformation = "log2(TPM + 1)",
      gene_scaling = scale_genes,
      pc_scaling = TRUE,
      pca_variance_threshold = pca_variance_threshold,
      pc_selection = "smallest number of PCs reaching the requested cumulative variance"
    )
  )
}
