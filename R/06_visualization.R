# Visualization ---------------------------------------------------------------

resolve_evaluation <- function(x) {
  if (inherits(x, "Tumor2LineEvaluation")) return(x)
  if ((inherits(x, "Tumor2LineResult") || inherits(x, "Tumor2LineReference")) &&
      !is.null(x$evaluation)) {
    return(x$evaluation)
  }
  if (inherits(x, "Tumor2LinePrediction") && !is.null(x$evaluation)) {
    return(x$evaluation)
  }
  stop(
    "Provide a Tumor2LineEvaluation or a result with known-pair evaluation.",
    call. = FALSE
  )
}

#' Plot Top1, Top3 and MRR performance
#'
#' @param x A `Tumor2LineEvaluation` or evaluated Tumor2Line result.
#' @return A ggplot object.
#' @export
plot_tumor2line_accuracy <- function(x) {
  ev <- resolve_evaluation(x)
  dat <- tibble::tibble(
    Metric = factor(
      c("Top1 accuracy", "Top3 accuracy", "MRR"),
      levels = c("Top1 accuracy", "Top3 accuracy", "MRR")
    ),
    Value = c(
      ev$summary$Top1_accuracy,
      ev$summary$Top3_accuracy,
      ev$summary$MRR
    )
  )

  ggplot2::ggplot(dat, ggplot2::aes(x = .data$Metric, y = .data$Value)) +
    ggplot2::geom_col(width = 0.68) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.3f", .data$Value)),
      vjust = -0.4,
      size = 3.8
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      expand = ggplot2::expansion(mult = c(0, 0.08))
    ) +
    ggplot2::labs(x = NULL, y = "Performance", title = "Tumor2Line matching performance") +
    ggplot2::theme_classic(base_size = 12)
}

#' Plot the distribution of true-pair ranks
#'
#' @param x A `Tumor2LineEvaluation` or evaluated Tumor2Line result.
#' @return A ggplot object.
#' @export
plot_true_pair_ranks <- function(x) {
  ev <- resolve_evaluation(x)
  dat <- ev$pair_ranks

  ggplot2::ggplot(dat, ggplot2::aes(x = .data$true_rank)) +
    ggplot2::geom_histogram(binwidth = 1, boundary = 0.5, closed = "left") +
    ggplot2::scale_x_continuous(breaks = tumor2line_pretty_breaks(dat$true_rank)) +
    ggplot2::labs(
      x = "Rank of known paired cell line",
      y = "Number of queries",
      title = "True-pair rank distribution"
    ) +
    ggplot2::theme_classic(base_size = 12)
}

tumor2line_pretty_breaks <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NULL)
  unique(pretty(range(x), n = min(10, length(unique(x)))))
}

#' Plot top candidate cell lines for one query
#'
#' @param result A `Tumor2LineResult`, `Tumor2LineReference` or
#'   `Tumor2LinePrediction` object.
#' @param query_id Query sample ID.
#' @param n Number of top candidates to show.
#' @return A ggplot object.
#' @export
plot_top_matches <- function(result, query_id, n = 10) {
  dat <- get_top_matches(result, n = n, query_id = query_id)
  if (nrow(dat) == 0) {
    stop("query_id was not found in the ranking output.", call. = FALSE)
  }
  if (length(unique(dat$Query_ID)) != 1) {
    stop("plot_top_matches() requires exactly one query_id.", call. = FALSE)
  }
  if (!"matching_probability" %in% colnames(dat)) {
    stop("matching_probability is not available in the ranking output.", call. = FALSE)
  }

  dat$CellLine_ID <- factor(
    dat$CellLine_ID,
    levels = rev(dat$CellLine_ID[order(dat$rank)])
  )

  ggplot2::ggplot(
    dat,
    ggplot2::aes(x = .data$matching_probability, y = .data$CellLine_ID)
  ) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::labs(
      x = "Row-normalized transport mass (matching probability)",
      y = "Cell line",
      title = paste0("Top cell-line matches for ", query_id)
    ) +
    ggplot2::theme_classic(base_size = 12)
}
