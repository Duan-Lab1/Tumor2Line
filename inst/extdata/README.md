# Tumor2Line bundled example data

This directory contains a small PDMR-format demonstration dataset used by
`vignettes/Tumor2Line-tutorial.Rmd` and the package tests.

Files:

- `PDMR_tumor_expr.tsv`: query expression table for tumor samples in the demo.
- `PDMR_cellline_expr.tsv`: candidate cell-line expression table.
- `PDMR_phenotype.tsv`: sample annotations including `SampleID`, `PairID`,
  `SampleClass`, and `lineage`.
- `PDMR_pairs.tsv`: known query-cell-line pairs for demonstration/evaluation.

The files are intentionally small and schema-compatible with the development
PDMR workflow. They are demonstration data, not the complete study matrices.
Tumor2Line itself is not restricted to paired tumor/cell-line data: the query
population may also contain PDO, PDX, organoid, or other tumor-derived samples,
and `PairID` is optional.
