# Tumor2Line v0.1.3: build the frozen Celligner reference PCA and check package
# Run from the Tumor2Line package root after placing full input matrices in
# data-raw/raw/. package.R belongs in R/; this script belongs in data-raw/.

# build package
roxygen2::roxygenise(getwd())
devtools::check(document = FALSE)
devtools::build(binary = FALSE, manual = TRUE, quiet = FALSE)

# build pkgdown
pkgdown::build_site()

library(data.table)
library(tidyverse)
library(Tumor2Line)

raw_dir <- file.path("data-raw", "raw")
required_expression_files <- c(
  "Celligner_tumor_gene_expression.tsv",
  "Celligner_cellline_gene_expression.tsv",
  "Celligner_phenotype.tsv",
  "PDMR_tumor_gene_expression.tsv",
  "PDMR_cellline_gene_expression.tsv",
  "PDMR_phenotype.tsv"
)

celligner_tumor <- data.table::fread(
  file.path(raw_dir, "Celligner_tumor_gene_expression.tsv"),
  data.table = FALSE,
  check.names = FALSE,
  colClasses = list(character = "sampleID")
)

celligner_cellline <- data.table::fread(
  file.path(raw_dir, "Celligner_cellline_gene_expression.tsv"),
  data.table = FALSE,
  check.names = FALSE,
  colClasses = list(character = "sampleID")
)

celligner_phenotype <- data.table::fread(
  file.path(raw_dir, "Celligner_phenotype.tsv"),
  data.table = FALSE,
  check.names = FALSE,
  colClasses = list(character = "sampleID")
)

celligner_reference <- build_reference_pca(
  tumor_expr = celligner_tumor,
  cellline_expr = celligner_cellline,
  phenotype = celligner_phenotype,
  sample_id_col = "sampleID",
  type_col = "type",
  lineage_col = "lineage",
  tumor_label = "tumor",
  cellline_label = "CL",
  input_scale = "processed",
  n_variable_genes = 10000,
  pca_variance_threshold = 0.90,
  scale_genes = TRUE
)

saveRDS(
  celligner_reference,
  file.path(raw_dir, "Celligner_reference_pca.rds"),
  compress = "xz"
)

message(
  "Celligner reference PCA saved to inst/extdata/Celligner_reference_pca.rds. ",
  "Selected genes: ", length(celligner_reference$representation$selected_genes),
  "; selected PCs: ", celligner_reference$representation$n_pcs,
  "; reference samples: ", celligner_reference$parameters$n_reference_samples,
  "."
)