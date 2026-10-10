read_bundled_pdmr_example <- function() {
  query_path <- system.file(
    "extdata",
    "PDMR_tumor_expr.tsv",
    package = "Tumor2Line"
  )
  cellline_path <- system.file(
    "extdata",
    "PDMR_cellline_expr.tsv",
    package = "Tumor2Line"
  )
  phenotype_path <- system.file(
    "extdata",
    "PDMR_phenotype.tsv",
    package = "Tumor2Line"
  )

  stopifnot(nzchar(query_path), nzchar(cellline_path), nzchar(phenotype_path))

  list(
    query_expr = utils::read.delim(
      query_path,
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    cellline_expr = utils::read.delim(
      cellline_path,
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    phenotype = utils::read.delim(
      phenotype_path,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  )
}
