test_that("known-pair evaluation returns Top1, Top3 and MRR", {
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
    max_pcs = NULL,
    candidate_scope = "within_lineage"
  )

  ev <- evaluate_tumor2line(result)

  expect_s3_class(ev, "Tumor2LineEvaluation")
  expect_true(all(c("Top1_accuracy", "Top3_accuracy", "MRR") %in% names(ev$summary)))
  expect_true(ev$summary$Top1_accuracy >= 0 && ev$summary$Top1_accuracy <= 1)
  expect_true(ev$summary$Top3_accuracy >= 0 && ev$summary$Top3_accuracy <= 1)
  expect_true(ev$summary$MRR > 0 && ev$summary$MRR <= 1)
  expect_true(all(c("Query_ID", "True_CellLine_ID", "true_rank") %in% names(ev$pair_ranks)))
})
