test_that("bundled PDMR-format data can be prepared", {
  demo <- read_bundled_pdmr_example()

  x <- prepare_tumor2line_data(
    query_expr = demo$query_expr,
    cellline_expr = demo$cellline_expr,
    phenotype = demo$phenotype,
    sample_id_col = "SampleID",
    pair_col = "PairID",
    lineage_col = "lineage",
    sample_type_col = "SampleClass",
    query_type = "tumor"
  )

  expect_s3_class(x, "Tumor2LineInput")
  expect_true(length(x$query_ids) > 0)
  expect_true(length(x$cellline_ids) > 0)
  expect_true(ncol(x$expression_tpm) >= 100)
  expect_true(all(c("SampleID", "SampleRole", "SampleType", "lineage", "PairID") %in% names(x$metadata)))
  expect_setequal(unique(x$metadata$SampleRole), c("query", "cellline"))
})

test_that("pairing and phenotype metadata are optional", {
  demo <- read_bundled_pdmr_example()

  x <- prepare_tumor2line_data(
    query_expr = demo$query_expr,
    cellline_expr = demo$cellline_expr,
    phenotype = NULL
  )

  expect_s3_class(x, "Tumor2LineInput")
  expect_equal(nrow(x$gold_pairs), 0)
  expect_true(all(is.na(x$metadata$PairID)))
  expect_true(all(is.na(x$metadata$lineage)))
})
