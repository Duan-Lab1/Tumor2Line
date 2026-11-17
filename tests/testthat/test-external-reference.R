test_that("external reference PCA projects raw and processed query data", {
  set.seed(11)
  n_gene <- 50
  tumor <- matrix(
    stats::rexp(30 * n_gene, rate = 0.2),
    nrow = 30,
    dimnames = list(paste0("RT", seq_len(30)), paste0("G", seq_len(n_gene)))
  )
  cellline <- matrix(
    stats::rexp(20 * n_gene, rate = 0.2),
    nrow = 20,
    dimnames = list(paste0("RC", seq_len(20)), paste0("G", seq_len(n_gene)))
  )
  reference_phenotype <- data.frame(
    sampleID = c(rownames(tumor), rownames(cellline)),
    type = c(rep("tumor", nrow(tumor)), rep("CL", nrow(cellline))),
    lineage = rep(c("colon", "lung"), length.out = nrow(tumor) + nrow(cellline)),
    stringsAsFactors = FALSE
  )

  reference <- build_reference_pca(
    tumor_expr = log2(tumor + 1),
    cellline_expr = log2(cellline + 1),
    phenotype = reference_phenotype,
    sample_id_col = "sampleID",
    type_col = "type",
    lineage_col = "lineage",
    input_scale = "processed",
    n_variable_genes = 40,
    pca_variance_threshold = 0.80
  )
  expect_s3_class(reference, "Tumor2LineReferencePCA")
  expect_equal(nrow(reference$cellline_scores), nrow(cellline))
  expect_true(reference$representation$n_pcs >= 1)
  reference_log <- log2(rbind(tumor, cellline) + 1)
  reference_sd <- apply(reference_log, 2, stats::sd)
  expected_genes <- names(reference_sd)[order(-reference_sd, names(reference_sd))][1:40]
  expect_identical(reference$representation$selected_genes, expected_genes)

  query <- matrix(
    stats::rexp(4 * n_gene, rate = 0.2),
    nrow = 4,
    dimnames = list(paste0("Q", 1:4), paste0("G", 1:n_gene))
  )
  projected <- project_to_reference_pca(query, reference, input_scale = "raw")
  expect_equal(nrow(projected$scores), nrow(query))
  expect_equal(ncol(projected$scores), reference$representation$n_pcs)

  query_phenotype <- data.frame(
    SampleID = rownames(query),
    Celligner_lineage = c("colon", "colon", "lung", "lung")
  )
  matched <- run_tumor2line_with_reference(
    query_expr = query,
    reference = reference,
    query_phenotype = query_phenotype,
    query_lineage_col = "Celligner_lineage",
    query_input_scale = "raw",
    candidate_scope = "within_lineage"
  )
  expect_s3_class(matched, "Tumor2LineResult")
  expect_true(all(matched$ranking$rank >= 1))
  expect_true("transport_mass" %in% names(matched$ranking))
})

test_that("lineage crosswalk rejects unmapped labels", {
  phenotype <- data.frame(TCGA_lineage = c("COADREAD", "PAAD"))
  crosswalk <- data.frame(
    TCGA_lineage = c("COADREAD", "PAAD"),
    Celligner_lineage = c("colorectal", "pancreas")
  )
  out <- apply_lineage_crosswalk(phenotype, crosswalk)
  expect_equal(out$Celligner_lineage, c("colorectal", "pancreas"))
  expect_error(
    apply_lineage_crosswalk(data.frame(TCGA_lineage = "SKCM"), crosswalk),
    "does not cover"
  )
})
