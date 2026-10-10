# Raw PDMR development files

Place the original PDMR development inputs here before running
`data-raw/debug.R`:

- `M3_paired_tumor_gene_expression.tsv`
- `M3_paired_cellline_gene_expression.tsv`
- `M3_paired_phenotype.tsv`

These source files are development inputs and are not required by package users.
`data-raw/debug.R` validates them with the public Tumor2Line preprocessing API
and writes standardized example resources to `inst/extdata/`.
