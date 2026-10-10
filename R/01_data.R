# Data preparation ------------------------------------------------------------

#' Prepare transcriptomic data for Tumor2Line
#'
#' Standardizes the two expression populations used by Tumor2Line: a query
#' population and a candidate cell-line population. The query population can
#' contain patient tumors, PDOs, PDXs, organoids, or other tumor-derived
#' samples. Known query-cell-line pairs are not required for Tumor2Line
#' matching and are used only for optional downstream evaluation.
#'
#' Expression matrices are harmonized to shared genes, transformed as
#' log2(TPM + 1), and summarized for variable-gene selection. Sample roles are
#' inferred directly from the two expression matrices rather than phenotype
#' labels.
#'
#' @param query_expr Expression data for the query population. Samples must be
#'   rows and genes must be columns. A data frame may contain a sample-ID column;
#'   a matrix may use sample IDs as row names.
#' @param cellline_expr Candidate cell-line expression data in the same
#'   orientation as `query_expr`.
#' @param phenotype Optional data frame containing sample annotations.
#' @param sample_id_col Name of the sample-ID column in `phenotype`.
#' @param pair_col Optional name of a known pairing column. Set to `NULL` when
#'   known query-cell-line pairs are unavailable.
#' @param lineage_col Optional name of the cancer-lineage column. Set to `NULL`
#'   when lineage information is unavailable.
#' @param sample_type_col Optional name of a column describing biological sample
#'   type, such as tumor, PDO, PDX, organoid, or cell line.
#' @param query_type Default biological label assigned to query samples when
#'   `sample_type_col` is unavailable.
#' @param cellline_type Default biological label assigned to candidate cell-line
#'   samples when `sample_type_col` is unavailable.
#' @param min_mean_tpm Minimum mean TPM required for a gene to be eligible for
#'   variable-gene selection.
#' @param min_detected_fraction Minimum fraction of samples with TPM > 0.
#'
#' @return An object of class `Tumor2LineInput` containing standardized metadata,
#'   TPM and log-expression matrices, gene statistics, sample IDs and optional
#'   known query-cell-line pairs.
#' @export
prepare_tumor2line_data <- function(
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
    min_detected_fraction = 0.1) {

  query_mat <- coerce_expression_matrix(
    query_expr,
    sample_id_col,
    "query_expr"
  )
  cellline_mat <- coerce_expression_matrix(
    cellline_expr,
    sample_id_col,
    "cellline_expr"
  )

  query_ids <- rownames(query_mat)
  cellline_ids <- rownames(cellline_mat)

  if (length(intersect(query_ids, cellline_ids)) > 0) {
    stop("Query and cell-line sample-ID sets overlap.", call. = FALSE)
  }
  if (length(query_ids) == 0) {
    stop("query_expr contains no samples.", call. = FALSE)
  }
  if (length(cellline_ids) == 0) {
    stop("cellline_expr contains no samples.", call. = FALSE)
  }

  common_genes <- intersect(
    colnames(query_mat),
    colnames(cellline_mat)
  )
  if (length(common_genes) < 2) {
    stop(
      "Fewer than two common genes were found between query and cell-line data.",
      call. = FALSE
    )
  }

  query_mat <- query_mat[, common_genes, drop = FALSE]
  cellline_mat <- cellline_mat[, common_genes, drop = FALSE]
  expression_tpm <- rbind(query_mat, cellline_mat)
  expression_log <- log2(expression_tpm + 1)

  all_ids <- c(query_ids, cellline_ids)
  metadata <- data.frame(
    SampleID = all_ids,
    SampleRole = c(
      rep("query", length(query_ids)),
      rep("cellline", length(cellline_ids))
    ),
    SampleType = c(
      rep(as.character(query_type)[1], length(query_ids)),
      rep(as.character(cellline_type)[1], length(cellline_ids))
    ),
    lineage = rep(NA_character_, length(all_ids)),
    PairID = rep(NA_character_, length(all_ids)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (is.null(phenotype)) {
    if (!is.null(pair_col)) {
      stop("pair_col cannot be used when phenotype is NULL.", call. = FALSE)
    }
    if (!is.null(lineage_col)) {
      stop("lineage_col cannot be used when phenotype is NULL.", call. = FALSE)
    }
    if (!is.null(sample_type_col)) {
      stop("sample_type_col cannot be used when phenotype is NULL.", call. = FALSE)
    }
  } else {
    phenotype <- as.data.frame(phenotype, check.names = FALSE)
    assert_columns(phenotype, sample_id_col, "phenotype")
    phenotype[[sample_id_col]] <- as.character(phenotype[[sample_id_col]])

    if (anyNA(phenotype[[sample_id_col]]) || any(phenotype[[sample_id_col]] == "")) {
      stop("phenotype contains missing or empty sample IDs.", call. = FALSE)
    }
    if (anyDuplicated(phenotype[[sample_id_col]])) {
      stop("Duplicated sample IDs were found in phenotype.", call. = FALSE)
    }

    phenotype_index <- match(metadata$SampleID, phenotype[[sample_id_col]])

    if (!is.null(lineage_col)) {
      assert_columns(phenotype, lineage_col, "phenotype")
      metadata$lineage <- as.character(
        phenotype[[lineage_col]][phenotype_index]
      )
    }

    if (!is.null(pair_col)) {
      assert_columns(phenotype, pair_col, "phenotype")
      metadata$PairID <- as.character(
        phenotype[[pair_col]][phenotype_index]
      )
    }

    if (!is.null(sample_type_col)) {
      assert_columns(phenotype, sample_type_col, "phenotype")
      sample_type <- as.character(
        phenotype[[sample_type_col]][phenotype_index]
      )
      use_sample_type <- !is.na(sample_type) & sample_type != ""
      metadata$SampleType[use_sample_type] <- sample_type[use_sample_type]
    }
  }

  expression_tpm <- expression_tpm[
    metadata$SampleID,
    ,
    drop = FALSE
  ]
  expression_log <- expression_log[
    metadata$SampleID,
    ,
    drop = FALSE
  ]

  gene_stats <- data.frame(
    gene = colnames(expression_log),
    mean_tpm = colMeans(expression_tpm),
    detected_fraction = colMeans(expression_tpm > 0),
    log_sd = apply(expression_log, 2, stats::sd),
    stringsAsFactors = FALSE
  )
  gene_stats <- gene_stats[
    is.finite(gene_stats$mean_tpm) &
      is.finite(gene_stats$detected_fraction) &
      is.finite(gene_stats$log_sd) &
      gene_stats$mean_tpm >= min_mean_tpm &
      gene_stats$detected_fraction >= min_detected_fraction &
      gene_stats$log_sd > 0,
    ,
    drop = FALSE
  ]
  gene_stats <- gene_stats[
    order(-gene_stats$log_sd, gene_stats$gene),
    ,
    drop = FALSE
  ]
  rownames(gene_stats) <- NULL

  gold_pairs <- data.frame(
    PairID = character(),
    lineage = character(),
    Query_ID = character(),
    CellLine_ID = character(),
    stringsAsFactors = FALSE
  )

  if (!is.null(pair_col) && any(!is.na(metadata$PairID) & metadata$PairID != "")) {
    query_meta <- metadata[
      metadata$SampleRole == "query" &
        !is.na(metadata$PairID) &
        metadata$PairID != "",
      c("PairID", "SampleID", "lineage"),
      drop = FALSE
    ]
    cellline_meta <- metadata[
      metadata$SampleRole == "cellline" &
        !is.na(metadata$PairID) &
        metadata$PairID != "",
      c("PairID", "SampleID", "lineage"),
      drop = FALSE
    ]

    names(query_meta) <- c("PairID", "Query_ID", "Query_lineage")
    names(cellline_meta) <- c("PairID", "CellLine_ID", "CellLine_lineage")

    pair_table <- merge(
      query_meta,
      cellline_meta,
      by = "PairID",
      all = FALSE
    )

    if (nrow(pair_table) > 0) {
      gold_pairs <- data.frame(
        PairID = pair_table$PairID,
        lineage = pair_table$Query_lineage,
        Query_ID = pair_table$Query_ID,
        CellLine_ID = pair_table$CellLine_ID,
        stringsAsFactors = FALSE
      )
      gold_pairs <- unique(gold_pairs)
      rownames(gold_pairs) <- NULL
    }
  }

  metadata_summary <- metadata
  metadata_summary$lineage[
    is.na(metadata_summary$lineage) |
      metadata_summary$lineage == ""
  ] <- "Unspecified"

  cohort_summary <- tibble::as_tibble(metadata_summary) |>
    dplyr::group_by(.data$lineage) |>
    dplyr::summarise(
      Query_number = sum(.data$SampleRole == "query"),
      Cell_line_number = sum(.data$SampleRole == "cellline"),
      Total_samples = dplyr::n(),
      .groups = "drop"
    )

  if (nrow(gold_pairs) > 0) {
    pair_summary <- tibble::as_tibble(gold_pairs) |>
      dplyr::mutate(
        lineage = dplyr::if_else(
          is.na(.data$lineage) | .data$lineage == "",
          "Unspecified",
          .data$lineage
        )
      ) |>
      dplyr::group_by(.data$lineage) |>
      dplyr::summarise(
        Paired_number = dplyr::n_distinct(.data$PairID),
        .groups = "drop"
      )

    cohort_summary <- cohort_summary |>
      dplyr::left_join(pair_summary, by = "lineage") |>
      dplyr::mutate(
        Paired_number = dplyr::coalesce(.data$Paired_number, 0L)
      )
  } else {
    cohort_summary$Paired_number <- 0L
  }

  cohort_summary <- cohort_summary |>
    dplyr::select(
      .data$lineage,
      .data$Query_number,
      .data$Cell_line_number,
      .data$Paired_number,
      .data$Total_samples
    ) |>
    dplyr::arrange(
      dplyr::desc(.data$Total_samples),
      .data$lineage
    )

  out <- list(
    metadata = metadata,
    phenotype = phenotype,
    expression_tpm = expression_tpm,
    expression_log = expression_log,
    gene_stats = gene_stats,
    common_genes = common_genes,
    gold_pairs = gold_pairs,
    cohort_summary = cohort_summary,
    query_ids = query_ids,
    cellline_ids = cellline_ids,
    input_columns = list(
      sample_id_col = sample_id_col,
      pair_col = pair_col,
      lineage_col = lineage_col,
      sample_type_col = sample_type_col
    )
  )

  class(out) <- c("Tumor2LineInput", "list")
  return(out)
}
