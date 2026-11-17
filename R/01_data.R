# Data preparation ------------------------------------------------------------

#' Prepare query and cell-line transcriptomic data
#'
#' Standardizes query expression, cell-line expression and optional phenotype
#' metadata. Samples without matching phenotype annotations are removed when
#' phenotype is supplied. Expression matrices are harmonized to common genes,
#' transformed to log2(TPM + 1), and summarized for variable-gene selection.
#' Known query-cell-line pairs are constructed when pairing information
#' is available.
#'
#' @param query_expr Query expression data. Samples must be rows and genes
#'   must be columns. A data frame may contain a sample-ID column; a matrix
#'   may use sample IDs as row names.
#' @param cellline_expr Cell-line expression data in the same orientation.
#' @param phenotype Optional data frame containing sample annotations.
#' @param sample_id_col Name of the sample-ID column.
#' @param pair_col Name of the known pairing column. Set to NULL if unavailable.
#' @param lineage_col Name of the cancer-lineage column. Set to NULL if unavailable.
#' @param sample_type_col Name of the biological sample-type column.
#' @param query_type Biological type assigned to query samples.
#' @param cellline_type Biological type assigned to candidate cell lines.
#' @param query_input_scale Expression scale of query data: "raw" (TPM)
#'   or "processed" (log2(TPM + 1)).
#' @param cellline_input_scale Expression scale of cell-line data: "raw" (TPM)
#'   or "processed" (log2(TPM + 1)).
#'
#' @return An object of class Tumor2LineInput containing standardized
#'   metadata, log-expression matrix, gene statistics, sample IDs and
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
    query_input_scale = c("raw", "processed"),
    cellline_input_scale = c("raw", "processed")) {
  
  query_input_scale <- match.arg(query_input_scale)
  cellline_input_scale <- match.arg(cellline_input_scale)
  
  query_mat <- coerce_expression_matrix(
    query_expr, sample_id_col, "query_expr"
  )
  
  cellline_mat <- coerce_expression_matrix(
    cellline_expr, sample_id_col, "cellline_expr"
  )
  
  query_ids <- rownames(query_mat)
  cellline_ids <- rownames(cellline_mat)
  
  if (length(query_ids) == 0 || length(cellline_ids) == 0) {
    stop(
      "Both query_expr and cellline_expr must contain at least one sample.",
      call. = FALSE
    )
  }
  
  if (length(intersect(query_ids, cellline_ids)) > 0) {
    stop(
      "Query and cell-line sample-ID sets overlap.",
      call. = FALSE
    )
  }
  
  if (is.null(phenotype)) {
    
    if (!is.null(pair_col) ||
        !is.null(lineage_col) ||
        !is.null(sample_type_col)) {
      stop(
        "pair_col, lineage_col and sample_type_col require phenotype.",
        call. = FALSE
      )
    }
    
  } else {
    
    phenotype <- as.data.frame(
      phenotype,
      check.names = FALSE
    )
    
    assert_columns(
      phenotype,
      sample_id_col,
      "phenotype"
    )
    
    phenotype[[sample_id_col]] <- as.character(
      phenotype[[sample_id_col]]
    )
    
    phenotype_ids <- phenotype[[sample_id_col]]
    
    if (anyNA(phenotype_ids) ||
        any(phenotype_ids == "") ||
        anyDuplicated(phenotype_ids)) {
      stop(
        paste0(
          "phenotype sample IDs must be non-missing, ",
          "non-empty and unique."
        ),
        call. = FALSE
      )
    }
    
    # Combined expression sample IDs
    expression_ids <- c(query_ids, cellline_ids)
    
    # Common sample IDs between expression and phenotype
    common_ids <- expression_ids[
      expression_ids %in% phenotype_ids
    ]
    
    # Identify samples without phenotype annotations
    dropped_query <- sum(!query_ids %in% common_ids)
    dropped_cellline <- sum(!cellline_ids %in% common_ids)
    
    # Retain only samples with matching phenotype annotations
    query_keep <- query_ids %in% common_ids
    cellline_keep <- cellline_ids %in% common_ids
    
    query_mat <- query_mat[query_keep, , drop = FALSE]
    cellline_mat <- cellline_mat[cellline_keep, , drop = FALSE]
    
    query_ids <- rownames(query_mat)
    cellline_ids <- rownames(cellline_mat)
    
    if (dropped_query > 0 || dropped_cellline > 0) {
      warning(
        sprintf(
          paste0(
            "Samples without matching phenotype annotations were removed: ",
            "%d query sample(s) and %d cell-line sample(s). ",
            "Retained: %d query sample(s) and %d cell-line sample(s)."
          ),
          dropped_query,
          dropped_cellline,
          length(query_ids),
          length(cellline_ids)
        ),
        call. = FALSE
      )
    }
    
    if (length(query_ids) == 0 || length(cellline_ids) == 0) {
      stop(
        paste0(
          "No valid samples remain for at least one sample group ",
          "after matching expression data with phenotype."
        ),
        call. = FALSE
      )
    }
    
    # Reorder phenotype according to the retained expression samples
    all_ids <- c(query_ids, cellline_ids)
    
    phenotype <- phenotype[
      match(all_ids, phenotype_ids),
      ,
      drop = FALSE
    ]
    
    rownames(phenotype) <- NULL
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
  
  query_norm <- normalize_expression_scale(
    query_mat,
    query_input_scale,
    "query_expr"
  )
  
  cellline_norm <- normalize_expression_scale(
    cellline_mat,
    cellline_input_scale,
    "cellline_expr"
  )
  
  expression_log <- rbind(
    query_norm$expression_log,
    cellline_norm$expression_log
  )
  
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
    lineage = NA_character_,
    PairID = NA_character_,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  if (!is.null(phenotype)) {
    
    phenotype_index <- match(
      metadata$SampleID,
      phenotype[[sample_id_col]]
    )
    
    if (anyNA(phenotype_index)) {
      stop(
        "Internal error: phenotype and expression sample IDs are inconsistent.",
        call. = FALSE
      )
    }
    
    if (!is.null(lineage_col)) {
      assert_columns(
        phenotype,
        lineage_col,
        "phenotype"
      )
      
      metadata$lineage <- as.character(
        phenotype[[lineage_col]][phenotype_index]
      )
    }
    
    if (!is.null(pair_col)) {
      assert_columns(
        phenotype,
        pair_col,
        "phenotype"
      )
      
      metadata$PairID <- as.character(
        phenotype[[pair_col]][phenotype_index]
      )
    }
    
    if (!is.null(sample_type_col)) {
      assert_columns(
        phenotype,
        sample_type_col,
        "phenotype"
      )
      
      sample_type <- as.character(
        phenotype[[sample_type_col]][phenotype_index]
      )
      
      use_sample_type <- !is.na(sample_type) & sample_type != ""
      
      metadata$SampleType[use_sample_type] <-
        sample_type[use_sample_type]
    }
  }
  
  gene_stats <- data.frame(
    gene = colnames(expression_log),
    log_sd = apply(expression_log, 2, stats::sd),
    stringsAsFactors = FALSE
  )
  
  gene_stats <- gene_stats[
    is.finite(gene_stats$log_sd) & gene_stats$log_sd > 0,
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
  
  if (!is.null(pair_col) &&
      any(!is.na(metadata$PairID) & metadata$PairID != "")) {
    
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
    
    names(query_meta) <- c(
      "PairID", "Query_ID", "Query_lineage"
    )
    
    names(cellline_meta) <- c(
      "PairID", "CellLine_ID", "CellLine_lineage"
    )
    
    pair_table <- merge(
      query_meta,
      cellline_meta,
      by = "PairID",
      all = FALSE
    )
    
    if (nrow(pair_table) > 0) {
      
      gold_pairs <- unique(data.frame(
        PairID = pair_table$PairID,
        lineage = pair_table$Query_lineage,
        Query_ID = pair_table$Query_ID,
        CellLine_ID = pair_table$CellLine_ID,
        stringsAsFactors = FALSE
      ))
      
      rownames(gold_pairs) <- NULL
    }
  }
  
  metadata_summary <- metadata
  
  metadata_summary$lineage[
    is.na(metadata_summary$lineage) |
      metadata_summary$lineage == ""
  ] <- "Unspecified"
  
  cohort_summary <- tibble::as_tibble(metadata_summary) |>
    dplyr::group_by(lineage) |>
    dplyr::summarise(
      Query_number = sum(SampleRole == "query"),
      Cell_line_number = sum(SampleRole == "cellline"),
      Total_samples = dplyr::n(),
      .groups = "drop"
    )
  
  if (nrow(gold_pairs) > 0) {
    
    pair_summary <- tibble::as_tibble(gold_pairs) |>
      dplyr::mutate(
        lineage = dplyr::if_else(
          is.na(lineage) | lineage == "",
          "Unspecified", lineage)) |>
      dplyr::group_by(lineage) |>
      dplyr::summarise(
        Paired_number = dplyr::n_distinct(PairID),
        .groups = "drop"
      )
    
    cohort_summary <- cohort_summary |>
      dplyr::left_join(pair_summary, by = "lineage") |>
      dplyr::mutate(
        Paired_number = dplyr::coalesce(Paired_number, 0L)
      )
    
  } else {
    
    cohort_summary$Paired_number <- 0L
  }
  
  cohort_summary <- cohort_summary |>
    dplyr::select(
      dplyr::all_of(c(
        "lineage",
        "Query_number",
        "Cell_line_number",
        "Paired_number",
        "Total_samples"
      ))
    ) |>
    dplyr::arrange(
      dplyr::desc(Total_samples),
      lineage
    )

  out <- list(
    metadata = metadata,
    phenotype = phenotype,
    expression_log = expression_log,
    gene_stats = gene_stats,
    common_genes = common_genes,
    gold_pairs = gold_pairs,
    cohort_summary = cohort_summary,
    query_ids = query_ids,
    cellline_ids = cellline_ids,
    input_scales = list(
      query = query_input_scale,
      cellline = cellline_input_scale
    ),
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
