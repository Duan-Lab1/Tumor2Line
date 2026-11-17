test_that("run_tumor2line returns transport mass and query-specific rankings", {
  demo <- read_bundled_pdmr_example()
  result <- run_tumor2line(
    query_expr = demo$query_expr,
    cellline_expr = demo$cellline_expr,
    phenotype = demo$phenotype,
    pair_col = "PairID",
    lineage_col = "lineage",
    sample_type_col = "SampleClass",
    query_type = "tumor",
    query_input_scale = "processed",
    cellline_input_scale = "processed",
    n_variable_genes = 200,
    pca_variance_threshold = 0.90,
    candidate_scope = "within_lineage"
  )

  expect_s3_class(result, "Tumor2LineResult")
  expect_gt(nrow(result$ranking), 0)
  expect_true(all(c(
    "Query_ID", "CellLine_ID", "transport_mass",
    "matching_probability", "cost", "rank"
  ) %in% names(result$ranking)))
  expect_false("score" %in% names(result$ranking))
  expect_true(all(result$top_matches$rank == 1L))
  expect_false(is.null(result$evaluation))

  probability_sums <- stats::aggregate(
    matching_probability ~ Query_ID,
    data = result$ranking,
    FUN = sum
  )$matching_probability
  expect_equal(probability_sums, rep(1, length(probability_sums)), tolerance = 1e-8)
})

test_that("run_tumor2line supports unpaired global matching", {
  demo <- read_bundled_pdmr_example()
  result <- run_tumor2line(
    query_expr = demo$query_expr,
    cellline_expr = demo$cellline_expr,
    phenotype = NULL,
    query_input_scale = "processed",
    cellline_input_scale = "processed",
    n_variable_genes = 200,
    candidate_scope = "global"
  )
  expect_s3_class(result, "Tumor2LineResult")
  expect_null(result$evaluation)
  expect_equal(length(unique(result$ranking$Query_ID)), nrow(demo$query_expr))
})

test_that("default Euclidean cost is mean squared L2 distance and is median-normalized", {
  query <- matrix(c(0, 0), nrow = 1, dimnames = list("Q1", c("G1", "G2")))
  cellline <- rbind(c(3, 4), c(0, 0))
  rownames(cellline) <- c("C1", "C2")

  raw_cost <- calculate_cost_matrix(query, cellline, normalize = FALSE)
  expect_equal(as.numeric(raw_cost[1, "C1"]), 12.5, tolerance = 1e-10)
  expect_equal(as.numeric(raw_cost[1, "C2"]), 0, tolerance = 1e-10)

  normalized_cost <- calculate_cost_matrix(query, cellline, normalize = TRUE)
  expect_equal(as.numeric(normalized_cost[1, "C1"]), 1, tolerance = 1e-10)
})
