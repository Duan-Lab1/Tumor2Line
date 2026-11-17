test_that("single-query prediction adds the query and reruns population uOT", {
  demo <- read_bundled_pdmr_example()
  heldout_id <- demo$query_expr$SampleID[1]
  heldout_query <- demo$query_expr[demo$query_expr$SampleID == heldout_id, , drop = FALSE]
  reference_query <- demo$query_expr[demo$query_expr$SampleID != heldout_id, , drop = FALSE]
  reference_ids <- c(reference_query$SampleID, demo$cellline_expr$SampleID)
  reference_phenotype <- demo$phenotype[
    demo$phenotype$SampleID %in% reference_ids,
    ,
    drop = FALSE
  ]
  new_query_metadata <- demo$phenotype[
    demo$phenotype$SampleID == heldout_id,
    ,
    drop = FALSE
  ]

  reference <- fit_tumor2line_reference(
    query_expr = reference_query,
    cellline_expr = demo$cellline_expr,
    phenotype = reference_phenotype,
    sample_id_col = "SampleID",
    pair_col = "PairID",
    lineage_col = "lineage",
    sample_type_col = "SampleClass",
    query_input_scale = "processed",
    cellline_input_scale = "processed",
    n_variable_genes = 200,
    pca_variance_threshold = 0.90,
    candidate_scope = "within_lineage"
  )

  prediction <- predict_tumor2line(
    reference = reference,
    new_query_expr = heldout_query,
    new_query_metadata = new_query_metadata,
    new_query_input_scale = "processed"
  )

  expect_s3_class(prediction, "Tumor2LinePrediction")
  expect_equal(unique(prediction$ranking$Query_ID), heldout_id)
  expect_true(all(c("transport_mass", "matching_probability", "rank") %in% names(prediction$ranking)))
  expect_false("method" %in% names(formals(predict_tumor2line)))
  expect_false("k" %in% names(formals(predict_tumor2line)))
  expect_equal(nrow(prediction$gold_pairs), 1)
  expect_false(is.null(prediction$evaluation))
})
