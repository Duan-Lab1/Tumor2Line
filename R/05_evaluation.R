# Evaluation -----------------------------------------------------------------

#' Evaluate Tumor2Line using known query-cell-line pairs
#'
#' Computes query-level true-pair ranks and summary metrics including Top1
#' accuracy, Top3 accuracy and mean reciprocal rank (MRR). Ground-truth pairs
#' are used only after matching and are not required by the matching algorithm.
#'
#' @param result A `Tumor2LineResult`, `Tumor2LineReference`,
#'   `Tumor2LinePrediction` object or ranking data frame.
#' @param true_pairs Optional data frame containing `Query_ID` and `CellLine_ID`.
#'   If omitted for a fitted result/reference, stored known pairs are used.
#' @return A `Tumor2LineEvaluation` object containing query-level true-pair
#'   ranks, overall performance and lineage-specific performance.
#' @export
evaluate_tumor2line <- function(result, true_pairs = NULL) {
  if (inherits(result, "Tumor2LineResult") || inherits(result, "Tumor2LineReference")) {
    ranking <- result$ranking
    if (is.null(true_pairs)) true_pairs <- result$input$gold_pairs
  } else if (inherits(result, "Tumor2LinePrediction")) {
    ranking <- result$ranking
    if (is.null(true_pairs)) true_pairs <- result$gold_pairs
  } else {
    ranking <- as.data.frame(result)
  }

  assert_columns(ranking, c("Query_ID", "CellLine_ID", "rank", "n_candidates"), "ranking")
  if (is.null(true_pairs)) {
    stop("true_pairs is required when known pairs are not stored in result.", call. = FALSE)
  }

  true_pairs <- as.data.frame(true_pairs, check.names = FALSE)
  assert_columns(true_pairs, c("Query_ID", "CellLine_ID"), "true_pairs")
  pair_cols <- intersect(c("Query_ID", "CellLine_ID", "PairID", "lineage"), colnames(true_pairs))
  pair_ref <- unique(true_pairs[, pair_cols, drop = FALSE])

  matched <- merge(
    ranking,
    pair_ref[, c("Query_ID", "CellLine_ID"), drop = FALSE],
    by = c("Query_ID", "CellLine_ID"),
    all = FALSE
  )
  if (nrow(matched) == 0) {
    stop("No known query-cell-line pairs were found in the ranking output.", call. = FALSE)
  }

  matched <- matched[order(matched$Query_ID, matched$rank, matched$CellLine_ID), , drop = FALSE]
  pair_ranks <- matched[!duplicated(matched$Query_ID), , drop = FALSE]
  pair_ranks$True_CellLine_ID <- pair_ranks$CellLine_ID
  pair_ranks$true_rank <- as.integer(pair_ranks$rank)
  pair_ranks$Top1 <- pair_ranks$true_rank == 1L
  pair_ranks$Top3 <- pair_ranks$true_rank <= 3L
  pair_ranks$reciprocal_rank <- 1 / pair_ranks$true_rank
  pair_ranks$normalized_rank <- pair_ranks$true_rank / pair_ranks$n_candidates
  pair_ranks$rank_percentile <- ifelse(
    pair_ranks$n_candidates > 1,
    (pair_ranks$true_rank - 1) / (pair_ranks$n_candidates - 1),
    0
  )

  keep <- c(
    "lineage", "Query_ID", "True_CellLine_ID", "true_rank", "Top1", "Top3",
    "reciprocal_rank", "normalized_rank", "rank_percentile", "n_candidates",
    "transport_mass", "matching_probability", "query_total_mass",
    "cellline_total_mass", "cost"
  )
  pair_ranks <- pair_ranks[, intersect(keep, colnames(pair_ranks)), drop = FALSE]

  summary <- tibble::tibble(
    n_query = length(unique(pair_ranks$Query_ID)),
    Top1_accuracy = mean(pair_ranks$true_rank == 1),
    Top3_accuracy = mean(pair_ranks$true_rank <= 3),
    MRR = mean(1 / pair_ranks$true_rank),
    median_true_rank = stats::median(pair_ranks$true_rank),
    mean_true_rank = mean(pair_ranks$true_rank)
  )

  if ("lineage" %in% colnames(pair_ranks)) {
    lineage <- tibble::as_tibble(pair_ranks) |>
      dplyr::group_by(.data$lineage) |>
      dplyr::summarise(
        n_query = dplyr::n_distinct(.data$Query_ID),
        n_candidates = as.integer(stats::median(.data$n_candidates)),
        Top1_accuracy = mean(.data$true_rank == 1),
        Top3_accuracy = mean(.data$true_rank <= 3),
        MRR = mean(1 / .data$true_rank),
        median_true_rank = stats::median(.data$true_rank),
        .groups = "drop"
      )
  } else {
    lineage <- tibble::tibble()
  }

  out <- list(
    pair_ranks = tibble::as_tibble(pair_ranks),
    summary = summary,
    lineage = lineage
  )
  class(out) <- c("Tumor2LineEvaluation", "list")
  out
}

#' Print a Tumor2Line evaluation
#'
#' @param x A `Tumor2LineEvaluation` object.
#' @param ... Additional arguments passed to the print method. Currently unused.
#' @return `x`, invisibly.
#' @export
print.Tumor2LineEvaluation <- function(x, ...) {
  cat("Tumor2Line known-pair evaluation\n")
  print(x$summary, row.names = FALSE)
  invisible(x)
}
