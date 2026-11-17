# Internal utilities ----------------------------------------------------------

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

assert_columns <- function(x, columns, object_name = "object") {
  columns <- unique(columns[!is.na(columns) & nzchar(columns)])
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
    if (is.null(rownames(out)) || anyNA(rownames(out)) || any(rownames(out) == "")) {
      stop(object_name, " matrix must have sample IDs as row names.", call. = FALSE)
    }
    if (is.null(colnames(out)) || anyNA(colnames(out)) || any(colnames(out) == "")) {
      stop(object_name, " matrix must have gene names as column names.", call. = FALSE)
    }
    if (anyDuplicated(rownames(out))) {
      stop(object_name, " contains duplicated sample IDs.", call. = FALSE)
    }
    if (anyDuplicated(colnames(out))) {
      stop(object_name, " contains duplicated gene names.", call. = FALSE)
    }
    suppressWarnings(storage.mode(out) <- "double")
  } else {
    x <- as.data.frame(x, check.names = FALSE)
    if (sample_id_col %in% colnames(x)) {
      sample_ids <- as.character(x[[sample_id_col]])
      gene_df <- x[, setdiff(colnames(x), sample_id_col), drop = FALSE]
    } else if (!is.null(rownames(x)) &&
               all(!is.na(rownames(x))) &&
               all(rownames(x) != "") &&
               !identical(rownames(x), as.character(seq_len(nrow(x))))) {
      sample_ids <- rownames(x)
      gene_df <- x
    } else {
      stop(
        object_name, " must contain '", sample_id_col,
        "' or have informative sample IDs as row names.",
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

    numeric_columns <- vapply(gene_df, is.numeric, logical(1))
    if (all(numeric_columns)) {
      out <- as.matrix(gene_df)
    } else {
      converted <- lapply(gene_df, function(v) {
        if (is.numeric(v)) return(as.numeric(v))
        suppressWarnings(as.numeric(as.character(v)))
      })
      out <- as.matrix(as.data.frame(converted, check.names = FALSE))
    }
    rownames(out) <- sample_ids
    colnames(out) <- colnames(gene_df)
    storage.mode(out) <- "double"
  }

  if (any(!is.finite(out))) {
    stop(object_name, " contains non-finite or non-numeric expression values.", call. = FALSE)
  }
  if (any(out < 0)) {
    stop(object_name, " must contain non-negative expression values.", call. = FALSE)
  }
  out
}

# Expression scale handling ---------------------------------------------------

#' Normalize expression input to log2(TPM + 1)
#'
#' `raw` means untransformed TPM, not raw RNA-seq counts. Raw counts must be
#' normalized to an appropriate expression unit such as TPM before they are
#' supplied to this function. `processed` means expression already transformed
#' as log2(TPM + 1).
#'
#' @param x Numeric expression matrix.
#' @param input_scale Either `"raw"` (TPM) or `"processed"` (log2(TPM + 1)).
#' @param object_name Name used in error messages.
#'
#' @return A list containing log-expression and the declared input scale.
#' @keywords internal
normalize_expression_scale <- function(
    x,
    input_scale = c("raw", "processed"),
    object_name = "expression data") {

  input_scale <- match.arg(input_scale)
  x <- as.matrix(x)

  if (any(!is.finite(x)) || any(x < 0)) {
    stop(
      object_name,
      " must contain finite, non-negative values. `processed` inputs must be log2(TPM + 1).",
      call. = FALSE
    )
  }

  expression_log <- if (input_scale == "raw") {
    log2(x + 1)
  } else {
    x
  }

  if (any(!is.finite(expression_log))) {
    stop(
      object_name,
      " produced non-finite values during expression-scale conversion.",
      call. = FALSE
    )
  }

  list(
    expression_log = expression_log,
    input_scale = input_scale
  )
}

validate_tumor2line_input <- function(x) {
  required <- c(
    "metadata", "expression_log", "gene_stats", "common_genes",
    "gold_pairs", "query_ids", "cellline_ids"
  )
  if (!is.list(x) || !all(required %in% names(x))) {
    stop("Invalid Tumor2Line prepared-data object.", call. = FALSE)
  }
  invisible(TRUE)
}
