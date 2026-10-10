# Single-query prediction -----------------------------------------------------

coerce_prediction_expression <- function(x, sample_id_col = "SampleID") {
  if (is.data.frame(x)) {
    if (sample_id_col %in% colnames(x)) {
      sample_ids <- as.character(x[[sample_id_col]])
      gene_data <- x[, setdiff(colnames(x), sample_id_col), drop = FALSE]
    } else {
      default_rownames <- identical(rownames(x), as.character(seq_len(nrow(x))))
      if (is.null(rownames(x)) || default_rownames) {
        stop(
          "new_query_expr must contain a sample-ID column or informative row names.",
          call. = FALSE
        )
      }
      sample_ids <- rownames(x)
      gene_data <- x
    }

    numeric_columns <- vapply(gene_data, is.numeric, logical(1))
    if (!all(numeric_columns)) {
      stop(
        "All gene-expression columns in new_query_expr must be numeric.",
        call. = FALSE
      )
    }

    mat <- as.matrix(gene_data)
  } else if (is.matrix(x)) {
    mat <- x
    if (is.null(rownames(mat))) {
      stop(
        "Matrix new_query_expr must contain sample IDs as row names.",
        call. = FALSE
      )
    }
    sample_ids <- rownames(mat)
    suppressWarnings(storage.mode(mat) <- "double")
  } else {
    stop(
      "new_query_expr must be a data frame or numeric matrix.",
      call. = FALSE
    )
  }

  if (length(sample_ids) != nrow(mat) || any(is.na(sample_ids) | sample_ids == "")) {
    stop("All new query samples must have non-empty sample IDs.", call. = FALSE)
  }
  if (anyDuplicated(sample_ids)) {
    stop("Duplicated sample IDs were found in new_query_expr.", call. = FALSE)
  }
  if (is.null(colnames(mat)) || any(is.na(colnames(mat)) | colnames(mat) == "")) {
    stop("new_query_expr must contain gene names as column names.", call. = FALSE)
  }
  if (anyDuplicated(colnames(mat))) {
    stop("Duplicated gene names were found in new_query_expr.", call. = FALSE)
  }
  if (!is.numeric(mat)) {
    stop("Gene-expression values in new_query_expr must be numeric.", call. = FALSE)
  }
  if (any(!is.finite(mat))) {
    stop("new_query_expr contains non-finite expression values.", call. = FALSE)
  }
  if (any(mat < 0)) {
    stop("new_query_expr must contain non-negative TPM-like expression values.", call. = FALSE)
  }

  rownames(mat) <- sample_ids
  mat
}

prepare_prediction_metadata <- function(reference, sample_ids, new_query_metadata = NULL) {
  sample_id_col <- reference$reference_columns$sample_id_col
  lineage_col <- reference$reference_columns$lineage_col
  sample_type_col <- reference$reference_columns$sample_type_col
  query_type <- reference$reference_columns$query_type
  candidate_scope <- reference$parameters$candidate_scope

  metadata <- data.frame(
    SampleID = sample_ids,
    SampleRole = rep("query", length(sample_ids)),
    SampleType = rep(query_type, length(sample_ids)),
    lineage = rep(NA_character_, length(sample_ids)),
    PairID = rep(NA_character_, length(sample_ids)),
    stringsAsFactors = FALSE
  )

  if (!is.null(new_query_metadata)) {
    new_query_metadata <- as.data.frame(new_query_metadata, check.names = FALSE)
    if (!sample_id_col %in% colnames(new_query_metadata)) {
      stop(
        "new_query_metadata does not contain the sample-ID column: ",
        sample_id_col,
        call. = FALSE
      )
    }

    idx <- match(sample_ids, as.character(new_query_metadata[[sample_id_col]]))
    if (any(is.na(idx))) {
      stop(
        "Some new query samples are missing from new_query_metadata: ",
        paste(sample_ids[is.na(idx)], collapse = ", "),
        call. = FALSE
      )
    }

    if (!is.null(lineage_col)) {
      if (!lineage_col %in% colnames(new_query_metadata)) {
        stop(
          "new_query_metadata does not contain the lineage column: ",
          lineage_col,
          call. = FALSE
        )
      }
      metadata$lineage <- as.character(new_query_metadata[[lineage_col]][idx])
    }

    if (!is.null(sample_type_col) && sample_type_col %in% colnames(new_query_metadata)) {
      sample_type <- as.character(new_query_metadata[[sample_type_col]][idx])
      keep <- !is.na(sample_type) & sample_type != ""
      metadata$SampleType[keep] <- sample_type[keep]
    }
  }

  if (candidate_scope == "within_lineage") {
    if (is.null(lineage_col)) {
      stop(
        "The reference was configured for within-lineage prediction but no lineage column was stored.",
        call. = FALSE
      )
    }
    if (any(is.na(metadata$lineage) | metadata$lineage == "")) {
      stop(
        "Lineage annotations are required for new queries when candidate_scope = \"within_lineage\".",
        call. = FALSE
      )
    }
  } else {
    metadata$lineage[is.na(metadata$lineage) | metadata$lineage == ""] <- "Unspecified"
  }

  metadata
}

project_new_query <- function(reference, new_query_expr) {
  sample_id_col <- reference$reference_columns$sample_id_col
  expression_tpm <- coerce_prediction_expression(
    new_query_expr,
    sample_id_col = sample_id_col
  )

  existing_ids <- reference$input$metadata$SampleID
  overlapping_ids <- intersect(rownames(expression_tpm), existing_ids)
  if (length(overlapping_ids) > 0) {
    stop(
      "New query sample IDs already exist in the reference: ",
      paste(overlapping_ids, collapse = ", "),
      call. = FALSE
    )
  }

  expression_log <- log2(expression_tpm + 1)
  projected <- project_expression_to_reference(
    expression_log = expression_log,
    reference_representation = reference$representation
  )

  list(
    expression_tpm = expression_tpm,
    expression_log = expression_log,
    scores = projected$scores,
    scores_raw = projected$scores_raw
  )
}

#' Fit a reusable Tumor2Line reference cohort
#'
#' Fits the standard Tumor2Line workflow to a reference query cohort and a
#' candidate cell-line population while retaining the gene-scaling parameters,
#' PCA loadings and PC-scaling parameters required to project future query
#' samples into the same molecular space.
#'
#' @param query_expr Expression data for the reference query population.
#' @param cellline_expr Candidate cell-line expression data.
#' @param phenotype Optional sample metadata.
#' @param sample_id_col Sample-ID column name in `phenotype` and data-frame
#'   expression inputs.
#' @param pair_col Optional known-pair column used only for evaluation.
#' @param lineage_col Optional cancer-lineage column name.
#' @param sample_type_col Optional biological sample-type column name.
#' @param query_type Default biological label assigned to query samples.
#' @param cellline_type Default biological label assigned to cell-line samples.
#' @param min_mean_tpm Minimum mean TPM for gene filtering.
#' @param min_detected_fraction Minimum detected fraction for gene filtering.
#' @param n_variable_genes Maximum number of variable genes.
#' @param pca_variance_threshold Cumulative PCA variance threshold.
#' @param max_pcs Optional maximum number of PCs to fit.
#' @param candidate_scope Candidate pool: `"within_lineage"` or `"global"`.
#' @param distance Distance metric used for uOT matching.
#' @param epsilon Entropic regularization parameter.
#' @param rho_query Query-population marginal-relaxation parameter.
#' @param rho_cellline Cell-line marginal-relaxation parameter.
#' @param max_iter Maximum Sinkhorn iterations.
#' @param tolerance Sinkhorn convergence tolerance.
#' @param evaluate Whether to evaluate known pairs when `pair_col` is supplied.
#'
#' @return An object of class `Tumor2LineReference`.
#' @export
fit_tumor2line_reference <- function(
    query_expr,
    cellline_expr,
    phenotype = NULL,
    sample_id_col = "SampleID",
    pair_col = NULL,
    lineage_col = NULL,
    sample_type_col = NULL,
    query_type = "query",
    cellline_type = "cellline",
    min_mean_tpm = 0.1,
    min_detected_fraction = 0.1,
    n_variable_genes = 10000,
    pca_variance_threshold = 0.90,
    max_pcs = NULL,
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

  result <- run_tumor2line(
    query_expr = query_expr,
    cellline_expr = cellline_expr,
    phenotype = phenotype,
    sample_id_col = sample_id_col,
    pair_col = pair_col,
    lineage_col = lineage_col,
    sample_type_col = sample_type_col,
    query_type = query_type,
    cellline_type = cellline_type,
    min_mean_tpm = min_mean_tpm,
    min_detected_fraction = min_detected_fraction,
    n_variable_genes = n_variable_genes,
    pca_variance_threshold = pca_variance_threshold,
    max_pcs = max_pcs,
    candidate_scope = candidate_scope,
    distance = distance,
    epsilon = epsilon,
    rho_query = rho_query,
    rho_cellline = rho_cellline,
    max_iter = max_iter,
    tolerance = tolerance,
    evaluate = evaluate
  )

  required_projection_fields <- c(
    "selected_genes", "gene_center", "gene_scale", "pca_rotation",
    "pc_center", "pc_scale", "n_pcs"
  )
  missing_projection_fields <- setdiff(
    required_projection_fields,
    names(result$representation)
  )
  if (length(missing_projection_fields) > 0) {
    stop(
      "The fitted representation does not contain all projection parameters: ",
      paste(missing_projection_fields, collapse = ", "),
      call. = FALSE
    )
  }

  result$reference_columns <- list(
    sample_id_col = sample_id_col,
    pair_col = pair_col,
    lineage_col = lineage_col,
    sample_type_col = sample_type_col,
    query_type = query_type,
    cellline_type = cellline_type
  )
  result$reference_query_ids <- result$input$query_ids
  result$reference_cellline_ids <- result$input$cellline_ids
  class(result) <- c("Tumor2LineReference", "Tumor2LineResult", "list")
  result
}

#' Print a Tumor2Line reference
#'
#' @param x A `Tumor2LineReference` object.
#' @param ... Additional arguments passed to the print method. Currently unused.
#' @return `x`, invisibly.
#' @keywords internal
#' @export
print.Tumor2LineReference <- function(x, ...) {
  cat("Tumor2Line reference\n")
  cat("  Reference queries: ", length(x$reference_query_ids), "\n", sep = "")
  cat("  Cell lines:        ", length(x$reference_cellline_ids), "\n", sep = "")
  cat("  Selected genes:    ", length(x$representation$selected_genes), "\n", sep = "")
  cat("  Selected PCs:      ", x$representation$n_pcs, "\n", sep = "")
  cat("  Candidate scope:   ", x$parameters$candidate_scope, "\n", sep = "")
  invisible(x)
}

predict_population_mode <- function(reference, projected, new_metadata) {
  candidate_scope <- reference$parameters$candidate_scope
  reference_metadata <- reference$input$metadata

  if (candidate_scope == "within_lineage") {
    affected_lineages <- unique(new_metadata$lineage)

    missing_query_lineages <- setdiff(
      affected_lineages,
      unique(reference_metadata$lineage[reference_metadata$SampleRole == "query"])
    )
    if (length(missing_query_lineages) > 0) {
      stop(
        "No reference queries are available for lineage(s): ",
        paste(missing_query_lineages, collapse = ", "),
        call. = FALSE
      )
    }

    missing_cellline_lineages <- setdiff(
      affected_lineages,
      unique(reference_metadata$lineage[reference_metadata$SampleRole == "cellline"])
    )
    if (length(missing_cellline_lineages) > 0) {
      stop(
        "No candidate cell lines are available for lineage(s): ",
        paste(missing_cellline_lineages, collapse = ", "),
        call. = FALSE
      )
    }

    reference_metadata_use <- reference_metadata[
      reference_metadata$lineage %in% affected_lineages,
      ,
      drop = FALSE
    ]
  } else {
    reference_metadata_use <- reference_metadata
  }

  reference_ids_use <- reference_metadata_use$SampleID
  reference_scores_use <- reference$representation$scores[
    reference_ids_use,
    ,
    drop = FALSE
  ]

  combined_scores <- rbind(
    reference_scores_use,
    projected$scores
  )
  combined_metadata <- rbind(
    reference_metadata_use,
    new_metadata
  )

  ranking_all <- run_uot_matching(
    representation = combined_scores,
    metadata = combined_metadata,
    candidate_scope = candidate_scope,
    distance = reference$parameters$distance,
    epsilon = reference$parameters$epsilon,
    rho_query = reference$parameters$rho_query,
    rho_cellline = reference$parameters$rho_cellline,
    max_iter = reference$parameters$max_iter,
    tolerance = reference$parameters$tolerance
  )

  ranking <- ranking_all[
    ranking_all$Query_ID %in% new_metadata$SampleID,
    ,
    drop = FALSE
  ]
  ranking$prediction_score <- ranking$matching_probability
  ranking$prediction_mode <- "population"
  ranking
}

predict_transfer_mode <- function(
    reference,
    projected,
    new_metadata,
    k = 5,
    neighbor_epsilon = 1e-8) {

  k <- as.integer(k)
  if (length(k) != 1 || is.na(k) || k < 1) {
    stop("k must be an integer >= 1.", call. = FALSE)
  }
  if (!is.finite(neighbor_epsilon) || neighbor_epsilon <= 0) {
    stop("neighbor_epsilon must be > 0.", call. = FALSE)
  }

  candidate_scope <- reference$parameters$candidate_scope
  reference_metadata <- reference$input$metadata
  reference_query_metadata <- reference_metadata[
    reference_metadata$SampleRole == "query",
    ,
    drop = FALSE
  ]

  ranking_out <- vector("list", nrow(new_metadata))
  neighbors_out <- vector("list", nrow(new_metadata))

  for (i in seq_len(nrow(new_metadata))) {
    query_id <- new_metadata$SampleID[i]
    query_lineage <- new_metadata$lineage[i]

    if (candidate_scope == "within_lineage") {
      neighbor_meta <- reference_query_metadata[
        reference_query_metadata$lineage == query_lineage,
        ,
        drop = FALSE
      ]
    } else {
      neighbor_meta <- reference_query_metadata
    }

    if (nrow(neighbor_meta) < 1) {
      stop(
        "No reference queries are available for transfer prediction of ",
        query_id,
        ".",
        call. = FALSE
      )
    }

    k_use <- min(k, nrow(neighbor_meta))
    if (k_use < k) {
      warning(
        "Only ",
        k_use,
        " reference queries are available for ",
        query_id,
        "; using all available reference queries.",
        call. = FALSE
      )
    }

    neighbor_scores <- reference$representation$scores[
      neighbor_meta$SampleID,
      ,
      drop = FALSE
    ]
    distance_matrix <- calculate_cost_matrix(
      query_matrix = projected$scores[query_id, , drop = FALSE],
      cellline_matrix = neighbor_scores,
      distance = reference$parameters$distance,
      normalize = FALSE
    )
    distances <- as.numeric(distance_matrix[1, ])
    names(distances) <- colnames(distance_matrix)

    ord <- order(distances, names(distances))
    selected_ids <- names(distances)[ord[seq_len(k_use)]]
    selected_distances <- distances[selected_ids]
    raw_weights <- 1 / (selected_distances + neighbor_epsilon)
    weights <- raw_weights / sum(raw_weights)

    neighbors <- tibble::tibble(
      Query_ID = query_id,
      Reference_Query_ID = selected_ids,
      lineage = query_lineage,
      distance = as.numeric(selected_distances),
      weight = as.numeric(weights),
      neighbor_rank = seq_len(k_use)
    )
    neighbors_out[[i]] <- neighbors

    reference_ranking <- reference$ranking[
      reference$ranking$Query_ID %in% selected_ids,
      ,
      drop = FALSE
    ]
    if (candidate_scope == "within_lineage") {
      reference_ranking <- reference_ranking[
        reference_ranking$lineage == query_lineage,
        ,
        drop = FALSE
      ]
    }

    if (nrow(reference_ranking) == 0) {
      stop(
        "No reference ranking information is available for transfer prediction of ",
        query_id,
        ".",
        call. = FALSE
      )
    }

    weight_lookup <- stats::setNames(weights, selected_ids)
    reference_ranking$neighbor_weight <- weight_lookup[reference_ranking$Query_ID]
    reference_ranking$weighted_probability <-
      reference_ranking$matching_probability * reference_ranking$neighbor_weight

    cell_ids <- unique(reference_ranking$CellLine_ID)
    transfer_score <- vapply(
      cell_ids,
      function(cell_id) {
        sum(
          reference_ranking$weighted_probability[
            reference_ranking$CellLine_ID == cell_id
          ],
          na.rm = TRUE
        )
      },
      numeric(1)
    )

    ord_cell <- order(-transfer_score, cell_ids)
    ranks <- match(seq_along(cell_ids), ord_cell)

    ranking_out[[i]] <- tibble::tibble(
      Method = "Tumor2Line transfer",
      lineage = query_lineage,
      Query_ID = query_id,
      CellLine_ID = cell_ids,
      transfer_score = as.numeric(transfer_score),
      prediction_score = as.numeric(transfer_score),
      rank = as.integer(ranks),
      n_candidates = length(cell_ids),
      n_reference_neighbors = k_use,
      candidate_scope = candidate_scope,
      prediction_mode = "transfer"
    )
  }

  list(
    ranking = dplyr::bind_rows(ranking_out) |>
      dplyr::arrange(.data$Query_ID, .data$rank),
    neighbors = dplyr::bind_rows(neighbors_out) |>
      dplyr::arrange(.data$Query_ID, .data$neighbor_rank)
  )
}

#' Predict candidate cell lines for a new query
#'
#' Projects a new query into a frozen Tumor2Line reference representation and
#' predicts candidate cell lines using either full population-aware uOT
#' recomputation (`method = "population"`) or nearest-reference ranking transfer
#' (`method = "transfer"`).
#'
#' Population mode adds the projected query to the corresponding reference
#' query population and reruns uOT without refitting variable genes, gene
#' scaling, PCA or PC scaling. Transfer mode identifies molecularly similar
#' reference queries and combines their existing `matching_probability` values
#' using inverse-distance weights.
#'
#' @param reference A `Tumor2LineReference` object.
#' @param new_query_expr Expression data for a single new query sample.
#' @param new_query_metadata Metadata for the new query sample. Required for
#'   within-lineage prediction.
#' @param method Prediction strategy: `"population"` or `"transfer"`.
#' @param k Number of similar reference queries used for transfer mode.
#' @param neighbor_epsilon Small positive constant used for inverse-distance
#'   weighting in transfer mode.
#'
#' @return An object of class `Tumor2LinePrediction` containing the new-query
#'   candidate ranking and, for transfer mode, the selected reference neighbors.
#' @export
predict_tumor2line <- function(
    reference,
    new_query_expr,
    new_query_metadata = NULL,
    method = c("population", "transfer"),
    k = 5,
    neighbor_epsilon = 1e-8) {

  if (!inherits(reference, "Tumor2LineReference")) {
    stop("reference must be a Tumor2LineReference object.", call. = FALSE)
  }

  method <- match.arg(method)
  projected <- project_new_query(
    reference = reference,
    new_query_expr = new_query_expr
  )

  if (nrow(projected$scores) != 1) {
    stop(
      "predict_tumor2line() currently supports exactly one new query sample at a time.",
      call. = FALSE
    )
  }

  new_metadata <- prepare_prediction_metadata(
    reference = reference,
    sample_ids = rownames(projected$scores),
    new_query_metadata = new_query_metadata
  )

  if (method == "population") {
    ranking <- predict_population_mode(
      reference = reference,
      projected = projected,
      new_metadata = new_metadata
    )
    neighbors <- NULL
  } else {
    transfer <- predict_transfer_mode(
      reference = reference,
      projected = projected,
      new_metadata = new_metadata,
      k = k,
      neighbor_epsilon = neighbor_epsilon
    )
    ranking <- transfer$ranking
    neighbors <- transfer$neighbors
  }

  top_matches <- ranking[
    ranking$rank == 1L,
    ,
    drop = FALSE
  ]

  out <- list(
    method = method,
    ranking = tibble::as_tibble(ranking),
    top_matches = tibble::as_tibble(top_matches),
    neighbors = neighbors,
    query_scores = projected$scores,
    query_scores_raw = projected$scores_raw,
    query_metadata = new_metadata,
    parameters = list(
      method = method,
      k = if (method == "transfer") as.integer(k) else NULL,
      neighbor_epsilon = if (method == "transfer") neighbor_epsilon else NULL,
      candidate_scope = reference$parameters$candidate_scope,
      distance = reference$parameters$distance
    )
  )
  class(out) <- c("Tumor2LinePrediction", "list")
  out
}

#' Print a Tumor2Line single-query prediction
#'
#' @param x A `Tumor2LinePrediction` object.
#' @param ... Additional arguments passed to the print method. Currently unused.
#' @return `x`, invisibly.
#' @keywords internal
#' @export
print.Tumor2LinePrediction <- function(x, ...) {
  cat("Tumor2Line single-query prediction\n")
  cat("  Mode:            ", x$method, "\n", sep = "")
  cat("  Query:           ", x$query_metadata$SampleID[1], "\n", sep = "")
  cat("  Lineage:         ", x$query_metadata$lineage[1], "\n", sep = "")
  cat("  Candidates:      ", nrow(x$ranking), "\n", sep = "")
  if (identical(x$method, "transfer")) {
    cat("  Reference k:     ", nrow(x$neighbors), "\n", sep = "")
  }
  invisible(x)
}
