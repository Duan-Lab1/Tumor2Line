test_that("expression samples are matched with phenotype IDs", {
  
  query_expr <- data.frame(
    SampleID = c("Q1", "Q2", "Q3"),
    G1 = c(1, 2, 3),
    G2 = c(2, 3, 4)
  )
  
  cellline_expr <- data.frame(
    SampleID = c("C1", "C2", "C3"),
    G1 = c(4, 5, 6),
    G2 = c(5, 6, 7)
  )
  
  phenotype <- data.frame(
    SampleID = c("Q1", "Q2", "C1", "C2", "C4"),
    lineage = "COAD"
  )
  
  expect_warning(
    result <- prepare_tumor2line_data(
      query_expr = query_expr,
      cellline_expr = cellline_expr,
      phenotype = phenotype,
      lineage_col = "lineage"
    ),
    "without matching phenotype"
  )
  
  expect_equal(result$query_ids, c("Q1", "Q2"))
  expect_equal(result$cellline_ids, c("C1", "C2"))
  
  expect_equal(
    result$metadata$SampleID,
    c("Q1", "Q2", "C1", "C2")
  )
  
  expect_equal(nrow(result$phenotype), 4)
})
