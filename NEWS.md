# Tumor2Line 0.1.3

* Added a reusable external reference PCA workflow through `build_reference_pca()`.
* Added `project_to_reference_pca()` and `run_tumor2line_with_reference()` for projecting new query and candidate samples into a frozen reference representation.

# Tumor2Line 0.1.2

* Generalized the core API from a tumor-only first population to a `query` population, allowing tumors, PDOs, PDXs, organoids, and other tumor-derived samples.
* Made known pairing information optional (`pair_col = NULL` by default); pairs are used only for evaluation.
* Added `SampleRole` (`query`/`cellline`) independently of biological `SampleType`.
* Renamed algorithm-facing identifiers from `Tumor_ID` to `Query_ID` and `rho_tumor` to `rho_query`.
* Changed `max_pcs` to `NULL` by default so the retained PC number is determined by the cumulative PCA variance threshold unless the user explicitly supplies an upper limit.
* Updated the tutorial to read bundled PDMR-format TSV files directly from `inst/extdata` with `system.file()`.
* Removed the former I/O helper script and tests that depended on it; package tests now read the bundled example files directly.

# Tumor2Line 0.1.1

* Added a complete user-facing workflow centered on `run_tumor2line()`.
* Added standardized data preparation for the three PDMR input tables.
* Added scaled-PCA representation and population-aware unbalanced optimal transport matching.
* Added evaluation using known tumor-cell line pairs (Top1, Top3, MRR and true-pair ranks).
* Added visualization helpers, package tests, pkgdown configuration and a complete vignette.