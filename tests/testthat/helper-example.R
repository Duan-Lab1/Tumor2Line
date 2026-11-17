read_bundled_pdmr_example <- function() {
  ext <- system.file("extdata", package = "Tumor2Line")
  stopifnot(nzchar(ext))
  list(
    query_expr = utils::read.delim(
      file.path(ext, "demo_tumor_expr.tsv"),
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    cellline_expr = utils::read.delim(
      file.path(ext, "demo_cellline_expr.tsv"),
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    phenotype = utils::read.delim(
      file.path(ext, "demo_phenotype.tsv"),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  )
}
