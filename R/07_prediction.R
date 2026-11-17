# Single-query / new-query population prediction ------------------------------

prepare_prediction_metadata <- function(reference, sample_ids, new_query_metadata = NULL) {
  reference_columns <- reference$reference_columns
  sample_id_col <- reference_columns$sample_id_col
  lineage_col <- reference_columns$lineage_col
  sample_type_col <- reference_columns$sample_type_col
  pair_col <- reference_columns$pair_col

  metadata <- data.frame(
    SampleID = sample_ids,
    SampleRole = "query",
    SampleType = reference_columns$query_type,
    lineage = NA_character_,
    PairID = NA_character_,
    stringsAsFactors = FALSE
  )

  if (!is.null(new_query_metadata)) {
    new_query_metadata <- as.data.frame(new_query_metadata, check.names = FALSE)
    assert_columns(new_query_metadata, sample_id_col, "new_query_metadata")
    new_ids <- as.character(new_query_metadata[[sample_id_col]])
    if (anyNA(new_ids) || any(new_ids == "") || anyDuplicated(new_ids)) {
      stop("new_query_metadata sample IDs must be non-missing, non-empty and unique.", call. = FALSE)
    }

    idx <- match(sample_ids, new_ids)
    if (anyNA(idx)) {
      stop(
        "Some new query samples are missing from new_query_metadata: ",
        paste(sample_ids[is.na(idx)], collapse = ", "),
        call. = FALSE
      )
    }

    if (!is.null(lineage_col)) {
      assert_columns(new_query_metadata, lineage_col, "new_query_metadata")
      metadata$lineage <- as.character(new_query_metadata[[lineage_col]][idx])
    }
    if (!is.null(sample_type_col) && sample_type_col %in% colnames(new_query_metadata)) {
      sample_type <- as.character(new_query_metadata[[sample_type_col]][idx])
      use_sample_type <- !is.na(sample_type) & sample_type != ""
      metadata$SampleType[use_sample_type] <- sample_type[use_sample_type]
    }
    if (!is.null(pair_col) && pair_col %in% colnames(new_query_metadata)) {
      metadata$PairID <- as.character(new_query_metadata[[pair_col]][idx])
    }
  }

  if (reference$parameters$candidate_scope == "within_lineage") {
    if (is.null(lineage_col)) {
      stop(
        "The fitted reference has no lineage column but uses within-lineage matching.",
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

project_new_query <- function(
    reference,
    new_query_expr,
    new_query_input_scale = c("raw", "processed")) {

  new_query_input_scale <- match.arg(new_query_input_scale)
  expression_input <- coerce_expression_matrix(
    new_query_expr,
    sample_id_col = reference$reference_columns$sample_id_col,
    object_name = "new_query_expr"
  )

  existing_ids <- reference$input$metadata$SampleID
  overlapping_ids <- intersect(rownames(expression_input), existing_ids)
  if (length(overlapping_ids) > 0) {
    stop(
      "New query sample IDs already exist in the fitted reference: ",
      paste(overlapping_ids, collapse = ", "),
      call. = FALSE
    )
  }

  normalized <- normalize_expression_scale(
    expression_input,
    new_query_input_scale,
    "new_query_expr"
  )
  projected <- project_expression_to_reference(
    expression_log = normalized$expression_log,
    reference_representation = reference$representation
  )

  list(
    expression_log = normalized$expression_log,
    input_scale = new_query_input_scale,
    scores = projected$scores,
    scores_raw = projected$scores_raw
  )
}

#' Fit a reusable Tumor2Line reference cohort
#'
#' Fits cohort-specific Tumor2Line PCA and uOT using an existing query cohort
#' and candidate cell-line population. The fitted gene scaling, PCA loadings
#' and PC scaling are retained so that newly arriving query samples can be
#' projected into the same molecular space. New queries are added to the query
#' population and uOT is rerun; no ranking-transfer mode is used.
#'
#' @inheritParams run_tumor2line
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
    query_input_scale = query_input_scale,
    cellline_input_scale = cellline_input_scale,
    n_variable_genes = n_variable_genes,
    pca_variance_threshold = pca_variance_threshold,
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
  missing_projection_fields <- setdiff(required_projection_fields, names(result$representation))
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

#' Print a Tumor2Line reference cohort
#'
#' @param x A `Tumor2LineReference` object.
#' @param ... Additional arguments passed to the print method. Currently unused.
#' @return `x`, invisibly.
#' @export
print.Tumor2LineReference <- function(x, ...) {
  cat("Tumor2Line reference cohort\n")
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
    reference_query_lineages <- unique(
      reference_metadata$lineage[reference_metadata$SampleRole == "query"]
    )
    reference_cellline_lineages <- unique(
      reference_metadata$lineage[reference_metadata$SampleRole == "cellline"]
    )

    missing_query_lineages <- setdiff(affected_lineages, reference_query_lineages)
    missing_cellline_lineages <- setdiff(affected_lineages, reference_cellline_lineages)
    if (length(missing_query_lineages) > 0 || length(missing_cellline_lineages) > 0) {
      missing <- unique(c(missing_query_lineages, missing_cellline_lineages))
      stop(
        "The fitted reference lacks query or candidate cell-line samples for lineage(s): ",
        paste(missing, collapse = ", "),
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
  combined_scores <- rbind(reference_scores_use, projected$scores)
  combined_metadata <- rbind(reference_metadata_use, new_metadata)

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

  ranking_all[ranking_all$Query_ID %in% new_metadata$SampleID, , drop = FALSE] |>
    dplyr::arrange(.data$Query_ID, .data$rank)
}

#' Predict candidate cell lines for new query samples using population mode
#'
#' Projects new query samples into the frozen molecular representation of a
#' fitted Tumor2Line reference cohort, adds them to the corresponding query
#' population and reruns uOT against the reference candidate cell-line
#' population. The PCA transformation is frozen, but the query population used
#' for uOT is updated to include the new samples.
#'
#' @param reference A `Tumor2LineReference` object created by
#'   [fit_tumor2line_reference()].
#' @param new_query_expr Expression data for one or more new query samples.
#' @param new_query_metadata Optional metadata for new query samples. Lineage
#'   information is required for within-lineage prediction.
#' @param new_query_input_scale Input scale: `"raw"` (TPM) or `"processed"`
#'   (log2(TPM + 1)).
#' @return A `Tumor2LinePrediction` object containing rankings for the new query
#'   samples and optional known-pair evaluation.
#' @export
predict_tumor2line <- function(
    reference,
    new_query_expr,
    new_query_metadata = NULL,
    new_query_input_scale = c("raw", "processed")) {

  if (!inherits(reference, "Tumor2LineReference")) {
    stop("reference must be a Tumor2LineReference object.", call. = FALSE)
  }

  projected <- project_new_query(
    reference = reference,
    new_query_expr = new_query_expr,
    new_query_input_scale = new_query_input_scale
  )
  new_metadata <- prepare_prediction_metadata(
    reference = reference,
    sample_ids = rownames(projected$scores),
    new_query_metadata = new_query_metadata
  )

  ranking <- predict_population_mode(
    reference = reference,
    projected = projected,
    new_metadata = new_metadata
  )
  top_matches <- ranking[ranking$rank == 1L, , drop = FALSE]

  gold_pairs <- data.frame(
    PairID = character(), lineage = character(),
    Query_ID = character(), CellLine_ID = character(),
    stringsAsFactors = FALSE
  )
  if (!is.null(reference$reference_columns$pair_col)) {
    query_pairs <- new_metadata[
      !is.na(new_metadata$PairID) & new_metadata$PairID != "",
      c("PairID", "SampleID", "lineage"),
      drop = FALSE
    ]
    cellline_meta <- reference$input$metadata[
      reference$input$metadata$SampleRole == "cellline" &
        !is.na(reference$input$metadata$PairID) &
        reference$input$metadata$PairID != "",
      c("PairID", "SampleID", "lineage"),
      drop = FALSE
    ]
    if (nrow(query_pairs) > 0 && nrow(cellline_meta) > 0) {
      names(query_pairs) <- c("PairID", "Query_ID", "lineage")
      names(cellline_meta) <- c("PairID", "CellLine_ID", "CellLine_lineage")
      pair_match <- merge(query_pairs, cellline_meta, by = "PairID", all = FALSE)
      if (nrow(pair_match) > 0) {
        gold_pairs <- unique(data.frame(
          PairID = pair_match$PairID,
          lineage = pair_match$lineage,
          Query_ID = pair_match$Query_ID,
          CellLine_ID = pair_match$CellLine_ID,
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  prediction <- list(
    method = "population",
    ranking = tibble::as_tibble(ranking),
    top_matches = tibble::as_tibble(top_matches),
    query_scores = projected$scores,
    query_scores_raw = projected$scores_raw,
    query_metadata = new_metadata,
    gold_pairs = gold_pairs,
    parameters = list(
      method = "population",
      new_query_input_scale = match.arg(new_query_input_scale),
      candidate_scope = reference$parameters$candidate_scope,
      distance = reference$parameters$distance,
      epsilon = reference$parameters$epsilon,
      rho_query = reference$parameters$rho_query,
      rho_cellline = reference$parameters$rho_cellline
    )
  )
  class(prediction) <- c("Tumor2LinePrediction", "list")

  if (nrow(gold_pairs) > 0) {
    prediction$evaluation <- evaluate_tumor2line(prediction, gold_pairs)
  } else {
    prediction$evaluation <- NULL
  }

  prediction
}

#' Print a Tumor2Line population-mode prediction
#'
#' @param x A `Tumor2LinePrediction` object.
#' @param ... Additional arguments passed to the print method. Currently unused.
#' @return `x`, invisibly.
#' @export
print.Tumor2LinePrediction <- function(x, ...) {
  cat("Tumor2Line population-mode prediction\n")
  cat("  New queries:     ", nrow(x$query_metadata), "\n", sep = "")
  cat("  Lineages:        ", length(unique(x$query_metadata$lineage)), "\n", sep = "")
  cat("  Ranking rows:    ", nrow(x$ranking), "\n", sep = "")
  if (!is.null(x$evaluation)) {
    cat("\nKnown-pair evaluation:\n")
    print(x$evaluation$summary, row.names = FALSE)
  }
  invisible(x)
}
