# High-level Tumor2Line workflow ---------------------------------------------

#' Run the Tumor2Line workflow
#'
#' Runs the standard cohort-specific Tumor2Line analysis from a query
#' population and candidate cell-line population. The workflow performs
#' common-gene harmonization, scale-aware conversion to log2(TPM + 1),
#' variable-gene selection, joint gene-wise scaling, PCA, PC selection based on
#' cumulative explained variance, PC-wise scaling and lineage-aware uOT.
#'
#' @param query_expr Expression data for the query population. Samples are rows
#'   and genes are columns.
#' @param cellline_expr Candidate cell-line expression data in the same orientation.
#' @param phenotype Optional sample metadata.
#' @param sample_id_col Sample-ID column name in expression data frames and phenotype.
#' @param pair_col Optional known-pair column name. Set to `NULL` if unavailable.
#' @param lineage_col Optional cancer-lineage column name. Required for
#'   `candidate_scope = "within_lineage"`.
#' @param sample_type_col Optional biological sample-type column name.
#' @param query_type Default biological label assigned to query samples.
#' @param cellline_type Default biological label assigned to candidate cell lines.
#' @param query_input_scale Input scale of query expression: `"raw"` (TPM) or
#'   `"processed"` (log2(TPM + 1)).
#' @param cellline_input_scale Input scale of candidate expression.
#' @param n_variable_genes Maximum number of variable genes.
#' @param pca_variance_threshold Cumulative PCA variance threshold used to
#'   determine the number of retained PCs.
#' @param candidate_scope Candidate pool: `"within_lineage"` or `"global"`.
#' @param distance Distance metric used to construct the uOT cost matrix.
#' @param epsilon Entropic regularization parameter.
#' @param rho_query Query marginal-relaxation parameter.
#' @param rho_cellline Cell-line marginal-relaxation parameter.
#' @param max_iter Maximum Sinkhorn iterations.
#' @param tolerance Sinkhorn convergence tolerance.
#' @param evaluate Whether to evaluate known pairs automatically when available.
#'
#' @return An object of class `Tumor2LineResult`.
#' @export
run_tumor2line <- function(
    query_expr,
    cellline_expr,
    phenotype = NULL,
    sample_id_col = "SampleID",
    pair_col = NULL,
    lineage_col = NULL,
    sample_type_col = NULL,
    query_type = "query",
    cellline_type = "cellline",
    query_input_scale = c("raw", "processed"),
    cellline_input_scale = c("raw", "processed"),
    n_variable_genes = 10000,
    pca_variance_threshold = 0.90,
    candidate_scope = c("within_lineage", "global"),
    distance = c("euclidean", "manhattan", "cosine"),
    epsilon = 0.1,
    rho_query = 1,
    rho_cellline = 1,
    max_iter = 5000,
    tolerance = 1e-8,
    evaluate = TRUE) {

  candidate_scope <- match.arg(candidate_scope)
  distance <- match.arg(distance)
  query_input_scale <- match.arg(query_input_scale)
  cellline_input_scale <- match.arg(cellline_input_scale)

  input <- prepare_tumor2line_data(
    query_expr = query_expr,
    cellline_expr = cellline_expr,
    phenotype = phenotype,
    sample_id_col = sample_id_col,
    pair_col = pair_col,
    lineage_col = lineage_col,
    sample_type_col = sample_type_col,
    query_type = query_type,
    cellline_type = cellline_type,
    query_input_scale = query_input_scale,
    cellline_input_scale = cellline_input_scale
  )

  if (candidate_scope == "within_lineage" && is.null(lineage_col)) {
    stop(
      "lineage_col is required when candidate_scope = \"within_lineage\".",
      call. = FALSE
    )
  }

  representation <- build_scaled_pca_representation(
    input = input,
    n_variable_genes = n_variable_genes,
    pca_variance_threshold = pca_variance_threshold,
    scale_genes = TRUE
  )

  ranking <- run_uot_matching(
    representation = representation$scores,
    metadata = input$metadata,
    candidate_scope = candidate_scope,
    distance = distance,
    epsilon = epsilon,
    rho_query = rho_query,
    rho_cellline = rho_cellline,
    max_iter = max_iter,
    tolerance = tolerance
  )

  top_matches <- ranking |>
    dplyr::filter(.data$rank == 1L) |>
    dplyr::arrange(.data$lineage, .data$Query_ID)

  result <- list(
    ranking = ranking,
    top_matches = top_matches,
    input = input,
    representation = representation,
    parameters = list(
      sample_id_col = sample_id_col,
      pair_col = pair_col,
      lineage_col = lineage_col,
      sample_type_col = sample_type_col,
      query_type = query_type,
      cellline_type = cellline_type,
      query_input_scale = query_input_scale,
      cellline_input_scale = cellline_input_scale,
      n_variable_genes = length(representation$selected_genes),
      selected_pcs = representation$n_pcs,
      pca_variance_threshold = pca_variance_threshold,
      candidate_scope = candidate_scope,
      distance = distance,
      epsilon = epsilon,
      rho_query = rho_query,
      rho_cellline = rho_cellline,
      max_iter = max_iter,
      tolerance = tolerance
    )
  )

  class(result) <- c("Tumor2LineResult", "list")
  if (isTRUE(evaluate) && nrow(input$gold_pairs) > 0) {
    result$evaluation <- evaluate_tumor2line(result, input$gold_pairs)
  } else {
    result$evaluation <- NULL
  }

  result
}

#' Print a Tumor2Line result
#'
#' @param x A `Tumor2LineResult` object.
#' @param ... Additional arguments passed to the print method. Currently unused.
#' @return `x`, invisibly.
#' @export
print.Tumor2LineResult <- function(x, ...) {
  cat("Tumor2Line result\n")
  cat("  Queries:         ", length(x$input$query_ids), "\n", sep = "")
  cat("  Cell lines:      ", length(x$input$cellline_ids), "\n", sep = "")

  lineages <- unique(x$input$metadata$lineage)
  lineages <- lineages[!is.na(lineages) & lineages != ""]
  cat("  Lineages:        ", length(lineages), "\n", sep = "")
  cat("  Selected genes:  ", length(x$representation$selected_genes), "\n", sep = "")
  cat("  Selected PCs:    ", x$representation$n_pcs, "\n", sep = "")
  cat("  PCA threshold:   ", x$parameters$pca_variance_threshold, "\n", sep = "")
  cat("  Candidate scope: ", x$parameters$candidate_scope, "\n", sep = "")
  cat("  Ranking rows:    ", nrow(x$ranking), "\n", sep = "")

  if (!is.null(x$evaluation)) {
    cat("\nKnown-pair evaluation:\n")
    print(x$evaluation$summary, row.names = FALSE)
  }

  invisible(x)
}

#' Get top-ranked cell lines for each query
#'
#' @param result A `Tumor2LineResult`, `Tumor2LineReference` or
#'   `Tumor2LinePrediction` object.
#' @param n Number of top candidates to return per query.
#' @param query_id Optional query sample ID. When supplied, only the requested
#'   query sample or samples are returned.
#' @return A tibble of top-ranked candidate cell lines.
#' @export
get_top_matches <- function(result, n = 5, query_id = NULL) {
  valid_result <- inherits(result, "Tumor2LineResult") ||
    inherits(result, "Tumor2LineReference") ||
    inherits(result, "Tumor2LinePrediction")
  if (!valid_result) {
    stop(
      "result must be a Tumor2LineResult, Tumor2LineReference or Tumor2LinePrediction object.",
      call. = FALSE
    )
  }

  n <- as.integer(n)
  if (length(n) != 1 || is.na(n) || n < 1) {
    stop("n must be an integer >= 1.", call. = FALSE)
  }

  x <- result$ranking
  if (!is.null(query_id)) {
    x <- x[x$Query_ID %in% query_id, , drop = FALSE]
  }

  x |>
    dplyr::filter(.data$rank <= n) |>
    dplyr::arrange(.data$Query_ID, .data$rank)
}
