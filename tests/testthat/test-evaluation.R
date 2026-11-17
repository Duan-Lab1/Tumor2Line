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
    query_input_scale = "processed",
    cellline_input_scale = "processed",
    n_variable_genes = 200,
    candidate_scope = "within_lineage"
  )
  evaluation <- evaluate_tumor2line(result)
  expect_s3_class(evaluation, "Tumor2LineEvaluation")
  expect_true(all(c("Top1_accuracy", "Top3_accuracy", "MRR") %in% names(evaluation$summary)))
  expect_true(all(c("Query_ID", "True_CellLine_ID", "true_rank") %in% names(evaluation$pair_ranks)))
})
