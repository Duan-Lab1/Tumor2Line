# Internal utilities ----------------------------------------------------------

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

assert_columns <- function(x, columns, object_name = "object") {
  columns <- columns[!is.na(columns) & nzchar(columns)]
  missing <- setdiff(columns, colnames(x))
  if (length(missing) > 0) {
    stop(
      object_name, " is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

coerce_expression_matrix <- function(
    x,
    sample_id_col = "SampleID",
    object_name = "expression data") {

  if (is.matrix(x)) {
    out <- x
    if (is.null(rownames(out)) || any(is.na(rownames(out))) || any(rownames(out) == "")) {
      stop(object_name, " matrix must have sample IDs as row names.", call. = FALSE)
    }
    if (is.null(colnames(out)) || any(is.na(colnames(out))) || any(colnames(out) == "")) {
      stop(object_name, " matrix must have gene names as column names.", call. = FALSE)
    }
    if (anyDuplicated(rownames(out))) {
      stop(object_name, " contains duplicated sample IDs.", call. = FALSE)
    }
    if (anyDuplicated(colnames(out))) {
      stop(object_name, " contains duplicated gene names.", call. = FALSE)
    }
    suppressWarnings(storage.mode(out) <- "double")
    if (any(!is.finite(out))) {
      stop(object_name, " contains non-finite or non-numeric expression values.", call. = FALSE)
    }
    if (any(out < 0)) {
      stop(object_name, " must contain non-negative expression values.", call. = FALSE)
    }
    return(out)
  }

  x <- as.data.frame(x, check.names = FALSE)
  if (sample_id_col %in% colnames(x)) {
    sample_ids <- as.character(x[[sample_id_col]])
    gene_df <- x[, setdiff(colnames(x), sample_id_col), drop = FALSE]
  } else if (!is.null(rownames(x)) && all(!is.na(rownames(x))) && all(rownames(x) != "")) {
    sample_ids <- rownames(x)
    gene_df <- x
  } else {
    stop(
      object_name, " must contain '", sample_id_col,
      "' or have sample IDs as row names.",
      call. = FALSE
    )
  }

  if (anyNA(sample_ids) || any(sample_ids == "")) {
    stop(object_name, " contains missing or empty sample IDs.", call. = FALSE)
  }
  if (anyDuplicated(sample_ids)) {
    stop(object_name, " contains duplicated sample IDs.", call. = FALSE)
  }
  if (ncol(gene_df) == 0) {
    stop(object_name, " contains no gene columns.", call. = FALSE)
  }
  if (anyDuplicated(colnames(gene_df))) {
    stop(object_name, " contains duplicated gene names.", call. = FALSE)
  }

  out <- as.matrix(data.frame(
    lapply(gene_df, function(v) suppressWarnings(as.numeric(as.character(v)))),
    check.names = FALSE
  ))
  rownames(out) <- sample_ids
  colnames(out) <- colnames(gene_df)
  storage.mode(out) <- "double"

  if (any(!is.finite(out))) {
    stop(object_name, " contains non-finite or non-numeric expression values.", call. = FALSE)
  }
  if (any(out < 0)) {
    stop(object_name, " must contain non-negative expression values.", call. = FALSE)
  }

  out
}

validate_tumor2line_input <- function(x) {
  required <- c(
    "metadata", "expression_tpm", "expression_log", "gene_stats",
    "common_genes", "gold_pairs", "query_ids", "cellline_ids"
  )
  if (!is.list(x) || !all(required %in% names(x))) {
    stop("Invalid Tumor2Line prepared-data object.", call. = FALSE)
  }
  invisible(TRUE)
}
