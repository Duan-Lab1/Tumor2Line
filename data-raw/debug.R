# build package
roxygen2::roxygenise(getwd())
devtools::check(document = FALSE)
devtools::build(binary = FALSE, manual = TRUE, quiet = FALSE)

# build pkgdown
pkgdown::build_site()

# PDMR example-data preparation ----------------------------------------------
library(Tumor2Line)

raw_data_dir <- "data-raw/raw"
extdata_dir <- "inst/extdata"

dir.create(extdata_dir, recursive = TRUE, showWarnings = FALSE)

query_file <- file.path(
  raw_data_dir,
  "M3_paired_tumor_gene_expression.tsv"
)
cellline_file <- file.path(
  raw_data_dir,
  "M3_paired_cellline_gene_expression.tsv"
)
phenotype_file <- file.path(
  raw_data_dir,
  "M3_paired_phenotype.tsv"
)

stopifnot(file.exists(query_file))
stopifnot(file.exists(cellline_file))
stopifnot(file.exists(phenotype_file))

phenotype <- utils::read.delim(
  phenotype_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
) |>
  dplyr::filter(model_lineage %in% c("COAD", "HNSC", "OV", "PAAD", "SKCM", "UCEC"))

table(phenotype$model_lineage, phenotype$SampleType)
phenotype$lineage <- phenotype$model_lineage

query_expr <- utils::read.delim(
  query_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
) |>
  dplyr::filter(SampleID %in% phenotype$SampleID)
cellline_expr <- utils::read.delim(
  cellline_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
) |>
  dplyr::filter(SampleID %in% phenotype$SampleID)


cat("\nQuery expression:\n")
print(dim(query_expr))

cat("\nCell-line expression:\n")
print(dim(cellline_expr))

cat("\nPhenotype:\n")
print(dim(phenotype))
print(colnames(phenotype))

prepared <- Tumor2Line::prepare_tumor2line_data(
  query_expr = query_expr,
  cellline_expr = cellline_expr,
  phenotype = phenotype,
  sample_id_col = "SampleID",
  pair_col = "PairID",
  lineage_col = "lineage",
  sample_type_col = "SampleClass",
  query_type = "tumor",
  min_mean_tpm = 0.1,
  min_detected_fraction = 0.1
)

PDMR_tumor_expr <- query_expr
PDMR_cellline_expr <- cellline_expr
PDMR_phenotype <- phenotype
PDMR_pairs <- prepared$gold_pairs

utils::write.table(
  PDMR_tumor_expr,
  file.path(extdata_dir, "PDMR_tumor_expr.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  PDMR_cellline_expr,
  file.path(extdata_dir, "PDMR_cellline_expr.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  PDMR_phenotype,
  file.path(extdata_dir, "PDMR_phenotype.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  PDMR_pairs,
  file.path(extdata_dir, "PDMR_pairs.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n============================================================\n")
cat("Tumor2Line PDMR example data generated successfully\n")
cat("============================================================\n")
cat("Query samples:    ", length(prepared$query_ids), "\n", sep = "")
cat("Cell lines:       ", length(prepared$cellline_ids), "\n", sep = "")
cat("Common genes:     ", length(prepared$common_genes), "\n", sep = "")
cat("Known pairs:      ", nrow(prepared$gold_pairs), "\n", sep = "")
cat("Output directory: ", extdata_dir, "\n", sep = "")
cat("============================================================\n")
