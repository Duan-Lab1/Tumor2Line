test_that("bundled PDMR-format demonstration data can be prepared", {
  demo <- read_bundled_pdmr_example()
  x <- prepare_tumor2line_data(
    query_expr = demo$query_expr,
    cellline_expr = demo$cellline_expr,
    phenotype = demo$phenotype,
    sample_id_col = "SampleID",
    pair_col = "PairID",
    lineage_col = "lineage",
    sample_type_col = "SampleClass",
    query_type = "tumor",
    query_input_scale = "processed",
    cellline_input_scale = "processed"
  )

  expect_s3_class(x, "Tumor2LineInput")
  expect_gt(length(x$query_ids), 0)
  expect_gt(length(x$cellline_ids), 0)
  expect_gt(nrow(x$gene_stats), 0)
  expect_named(x$gene_stats, c("gene", "log_sd"))
  expect_false("expression_tpm" %in% names(x))
  expect_equal(nrow(x$gold_pairs), length(x$query_ids))
  expect_setequal(unique(x$metadata$SampleRole), c("query", "cellline"))
})

test_that("raw TPM and processed log2(TPM + 1) are handled exactly once", {
  tpm <- matrix(
    c(0, 1, 9, 3, 15, 31, 7, 4, 2, 8, 11, 5, 0, 0, 0, 0.01),
    nrow = 4,
    dimnames = list(paste0("S", 1:4), c(paste0("G", 1:3), "LOW"))
  )
  log_expr <- log2(tpm + 1)

  raw_input <- prepare_tumor2line_data(
    tpm[1:2, , drop = FALSE],
    tpm[3:4, , drop = FALSE],
    query_input_scale = "raw",
    cellline_input_scale = "raw"
  )
  processed_input <- prepare_tumor2line_data(
    log_expr[1:2, , drop = FALSE],
    log_expr[3:4, , drop = FALSE],
    query_input_scale = "processed",
    cellline_input_scale = "processed"
  )

  expect_equal(raw_input$expression_log, processed_input$expression_log, tolerance = 1e-10)
  expect_true("LOW" %in% raw_input$gene_stats$gene)
  expect_false(any(c("mean_tpm", "detected_fraction") %in% names(raw_input$gene_stats)))
  expect_false("expression_tpm" %in% names(raw_input))
})

test_that("pairing and phenotype metadata are optional", {
  demo <- read_bundled_pdmr_example()
  x <- prepare_tumor2line_data(
    query_expr = demo$query_expr,
    cellline_expr = demo$cellline_expr,
    phenotype = NULL,
    query_input_scale = "processed",
    cellline_input_scale = "processed"
  )
  expect_equal(nrow(x$gold_pairs), 0)
  expect_true(all(is.na(x$metadata$PairID)))
  expect_true(all(is.na(x$metadata$lineage)))
})

test_that("only raw TPM and processed log2(TPM + 1) scales are accepted", {
  demo <- read_bundled_pdmr_example()
  expect_error(
    prepare_tumor2line_data(
      query_expr = demo$query_expr,
      cellline_expr = demo$cellline_expr,
      query_input_scale = "counts",
      cellline_input_scale = "processed"
    ),
    "should be one of"
  )
})
