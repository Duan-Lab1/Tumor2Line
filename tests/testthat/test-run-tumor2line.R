test_that("run_tumor2line returns query-specific candidate rankings", {
  demo <- read_bundled_pdmr_example()

  result <- run_tumor2line(
    query_expr = demo$query_expr,
    cellline_expr = demo$cellline_expr,
    phenotype = demo$phenotype,
    pair_col = "PairID",
    lineage_col = "lineage",
    sample_type_col = "SampleClass",
    query_type = "tumor",
    n_variable_genes = 200,
    pca_variance_threshold = 0.90,
    max_pcs = NULL,
    candidate_scope = "within_lineage"
  )

  expect_s3_class(result, "Tumor2LineResult")
  expect_true(nrow(result$ranking) > 0)
  expect_true(all(c("Query_ID", "CellLine_ID", "matching_probability", "rank") %in% names(result$ranking)))
  expect_true(all(result$top_matches$rank == 1L))
  expect_false(is.null(result$evaluation))
  expect_true(result$representation$n_pcs >= 1)
})

test_that("run_tumor2line supports unpaired global matching", {
  demo <- read_bundled_pdmr_example()

  result <- run_tumor2line(
    query_expr = demo$query_expr,
    cellline_expr = demo$cellline_expr,
    phenotype = NULL,
    n_variable_genes = 200,
    pca_variance_threshold = 0.90,
    max_pcs = NULL,
    candidate_scope = "global"
  )

  expect_s3_class(result, "Tumor2LineResult")
  expect_null(result$evaluation)
  expect_equal(length(unique(result$ranking$Query_ID)), nrow(demo$query_expr))
})
