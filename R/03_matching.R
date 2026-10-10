# Matching algorithms ---------------------------------------------------------

calculate_cost_matrix <- function(
    query_matrix,
    cellline_matrix,
    distance = c("euclidean", "manhattan", "cosine"),
    normalize = TRUE) {
  
  distance <- match.arg(distance)
  if (ncol(query_matrix) != ncol(cellline_matrix)) {
    stop("Query and cell-line representations must have the same dimensions.", call. = FALSE)
  }
  
  if (distance == "euclidean") {
    cost <- outer(rowSums(query_matrix^2), rowSums(cellline_matrix^2), "+") -
      2 * tcrossprod(query_matrix, cellline_matrix)
    cost <- pmax(cost / ncol(query_matrix), 0)
  } else if (distance == "manhattan") {
    cost <- vapply(seq_len(nrow(cellline_matrix)), function(j) {
      rowMeans(abs(sweep(query_matrix, 2, cellline_matrix[j, ], "-")))
    }, numeric(nrow(query_matrix)))
    if (is.null(dim(cost))) {
      cost <- matrix(
        cost,
        nrow = nrow(query_matrix),
        ncol = nrow(cellline_matrix)
      )
    }
  } else {
    query_norm <- sqrt(rowSums(query_matrix^2))
    cell_norm <- sqrt(rowSums(cellline_matrix^2))
    query_unit <- query_matrix / pmax(query_norm, 1e-12)
    cell_unit <- cellline_matrix / pmax(cell_norm, 1e-12)
    cost <- pmax(1 - tcrossprod(query_unit, cell_unit), 0)
  }
  
  rownames(cost) <- rownames(query_matrix)
  colnames(cost) <- rownames(cellline_matrix)
  
  if (normalize) {
    positive <- cost[is.finite(cost) & cost > 0]
    scale_value <- if (length(positive)) stats::median(positive) else 1
    if (!is.finite(scale_value) || scale_value <= 0) scale_value <- 1
    cost <- cost / scale_value
  }
  
  cost
}

run_unbalanced_sinkhorn <- function(
    cost_matrix,
    epsilon = 0.1,
    rho_query = 1,
    rho_cellline = 1,
    max_iter = 5000,
    tolerance = 1e-8) {
  
  if (!is.finite(epsilon) || epsilon <= 0) {
    stop("epsilon must be > 0.", call. = FALSE)
  }
  if (!is.finite(rho_query) || rho_query <= 0) {
    stop("rho_query must be > 0.", call. = FALSE)
  }
  if (!is.finite(rho_cellline) || rho_cellline <= 0) {
    stop("rho_cellline must be > 0.", call. = FALSE)
  }
  
  n_query <- nrow(cost_matrix)
  n_cellline <- ncol(cost_matrix)
  if (n_query < 1 || n_cellline < 1) {
    stop("The cost matrix must contain at least one query and one cell line.", call. = FALSE)
  }
  
  a <- rep(1 / n_query, n_query)
  b <- rep(1 / n_cellline, n_cellline)
  kernel <- pmax(exp(-cost_matrix / epsilon), 1e-300)
  tau_a <- rho_query / (rho_query + epsilon)
  tau_b <- rho_cellline / (rho_cellline + epsilon)
  
  u <- rep(1, n_query)
  v <- rep(1, n_cellline)
  delta <- Inf
  
  for (iteration in seq_len(max_iter)) {
    u_new <- (a / pmax(as.vector(kernel %*% v), 1e-300))^tau_a
    v_new <- (b / pmax(as.vector(crossprod(kernel, u_new)), 1e-300))^tau_b
    u_new <- pmin(pmax(u_new, 1e-100), 1e100)
    v_new <- pmin(pmax(v_new, 1e-100), 1e100)
    
    delta <- max(
      abs(log(u_new) - log(u)),
      abs(log(v_new) - log(v))
    )
    u <- u_new
    v <- v_new
    if (is.finite(delta) && delta < tolerance) break
  }
  
  plan <- sweep(kernel, 1, u, "*")
  plan <- sweep(plan, 2, v, "*")
  rownames(plan) <- rownames(cost_matrix)
  colnames(plan) <- colnames(cost_matrix)
  
  list(
    transport_plan = plan,
    converged = is.finite(delta) && delta < tolerance,
    iterations = iteration,
    final_delta = delta
  )
}

#' Run population-aware unbalanced optimal transport matching
#'
#' Matches a query population and candidate cell-line population jointly within
#' cancer lineage or across a global candidate pool using an entropically
#' regularized unbalanced Sinkhorn algorithm.
#'
#' The raw uOT transport plan is retained as `transport_mass`. For each query,
#' transport masses are additionally normalized across its candidate cell lines
#' to obtain `matching_probability`. Query-specific candidate ranks are based on
#' these normalized matching probabilities.
#'
#' @param representation Numeric matrix with samples in rows and molecular
#'   features in columns.
#' @param metadata Standardized metadata containing `SampleID`, `SampleRole`
#'   and `lineage`.
#' @param candidate_scope Either `"within_lineage"` or `"global"`.
#' @param distance Distance metric: `"euclidean"`, `"manhattan"` or `"cosine"`.
#' @param epsilon Entropic regularization parameter.
#' @param rho_query Query-population marginal-relaxation parameter.
#' @param rho_cellline Cell-line marginal-relaxation parameter.
#' @param max_iter Maximum Sinkhorn iterations.
#' @param tolerance Convergence tolerance.
#'
#' @return A tibble containing all query-cell-line candidate pairs, raw
#'   transport masses, query-normalized matching probabilities, query and
#'   cell-line total transported masses, molecular costs, query-specific ranks
#'   and uOT convergence information.
#' @export
run_uot_matching <- function(
    representation,
    metadata,
    candidate_scope = c("within_lineage", "global"),
    distance = c("euclidean", "manhattan", "cosine"),
    epsilon = 0.1,
    rho_query = 1,
    rho_cellline = 1,
    max_iter = 5000,
    tolerance = 1e-8) {
  
  candidate_scope <- match.arg(candidate_scope)
  distance <- match.arg(distance)
  representation <- as.matrix(representation)
  metadata <- as.data.frame(metadata, check.names = FALSE)
  assert_columns(metadata, c("SampleID", "SampleRole", "lineage"), "metadata")
  
  if (is.null(rownames(representation))) {
    stop("representation must contain sample IDs as row names.", call. = FALSE)
  }
  if (anyDuplicated(rownames(representation))) {
    stop("Duplicated sample IDs were found in representation.", call. = FALSE)
  }
  
  missing_samples <- setdiff(metadata$SampleID, rownames(representation))
  if (length(missing_samples) > 0) {
    stop(
      "Some metadata samples are missing from representation: ",
      paste(utils::head(missing_samples, 10), collapse = ", "),
      call. = FALSE
    )
  }
  
  if (candidate_scope == "within_lineage") {
    if (any(is.na(metadata$lineage) | metadata$lineage == "")) {
      stop(
        paste0(
          "Lineage annotations are required for candidate_scope = \"within_lineage\". ",
          "Provide lineage metadata or use candidate_scope = \"global\"."
        ),
        call. = FALSE
      )
    }
    groups <- split(metadata, metadata$lineage)
  } else {
    groups <- list(ALL = metadata)
  }
  
  out <- lapply(names(groups), function(group_name) {
    group_meta <- groups[[group_name]]
    query_ids <- intersect(
      group_meta$SampleID[group_meta$SampleRole == "query"],
      rownames(representation)
    )
    cell_ids <- intersect(
      group_meta$SampleID[group_meta$SampleRole == "cellline"],
      rownames(representation)
    )
    if (length(query_ids) < 1 || length(cell_ids) < 1) return(NULL)
    
    query <- representation[query_ids, , drop = FALSE]
    cell <- representation[cell_ids, , drop = FALSE]
    cost <- calculate_cost_matrix(query, cell, distance = distance)
    fit <- run_unbalanced_sinkhorn(
      cost,
      epsilon = epsilon,
      rho_query = rho_query,
      rho_cellline = rho_cellline,
      max_iter = max_iter,
      tolerance = tolerance
    )
    
    transport_plan <- fit$transport_plan
    query_total_mass <- rowSums(transport_plan)
    cellline_total_mass <- colSums(transport_plan)
    matching_probability <- transport_plan / pmax(query_total_mass, 1e-300)
    
    dplyr::bind_rows(lapply(seq_along(query_ids), function(i) {
      transport_mass_i <- transport_plan[i, , drop = TRUE]
      matching_probability_i <- matching_probability[i, , drop = TRUE]
      cost_i <- cost[i, , drop = TRUE]
      
      ord <- order(-matching_probability_i, cost_i, cell_ids)
      ranks <- match(seq_along(cell_ids), ord)
      
      query_lineage <- metadata$lineage[match(query_ids[i], metadata$SampleID)]
      if (candidate_scope == "global" && (is.na(query_lineage) || query_lineage == "")) {
        query_lineage <- "Unspecified"
      }
      
      tibble::tibble(
        Method = "Tumor2Line",
        lineage = query_lineage,
        Query_ID = query_ids[i],
        CellLine_ID = cell_ids,
        transport_mass = as.numeric(transport_mass_i),
        matching_probability = as.numeric(matching_probability_i),
        query_total_mass = as.numeric(query_total_mass[i]),
        cellline_total_mass = as.numeric(cellline_total_mass),
        cost = as.numeric(cost_i),
        rank = as.integer(ranks),
        n_candidates = length(cell_ids),
        candidate_scope = candidate_scope,
        OT_converged = fit$converged,
        OT_iterations = fit$iterations,
        OT_final_delta = fit$final_delta
      )
    }))
  })
  
  result <- dplyr::bind_rows(out)
  if (nrow(result) == 0) {
    stop("No query-cell-line candidate groups could be matched.", call. = FALSE)
  }
  
  ranked_queries <- unique(result$Query_ID)
  all_queries <- metadata$SampleID[metadata$SampleRole == "query"]
  unmatched_queries <- setdiff(all_queries, ranked_queries)
  if (length(unmatched_queries) > 0) {
    stop(
      "No candidate cell lines were available for query sample(s): ",
      paste(unmatched_queries, collapse = ", "),
      call. = FALSE
    )
  }
  
  result |>
    dplyr::arrange(.data$lineage, .data$Query_ID, .data$rank)
}