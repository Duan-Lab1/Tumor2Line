# External reference PCA ------------------------------------------------------

#' Build a reusable external reference PCA space
#'
#' Fits a joint scaled-PCA representation using an external tumor population
#' and an external cancer cell-line population. Shared genes are ranked by
#' standard deviation across the combined log-expression matrix, and the top
#' requested genes are retained. The selected genes, gene-wise scaling
#' parameters, PCA loadings, retained PC number and PC-wise scaling parameters
#' are stored so future query and candidate samples can be projected without
#' refitting the reference transformation.
#'
#' `input_scale = "raw"` means untransformed TPM; `input_scale = "processed"`
#' means data already transformed as log2(TPM + 1). Raw RNA-seq counts are not
#' accepted as TPM and must be normalized before using this function.
#'
#' @param tumor_expr External reference tumor expression data, samples in rows.
#' @param cellline_expr External reference cell-line expression data, samples in rows.
#' @param phenotype Reference phenotype table containing sample type and lineage.
#' @param sample_id_col Sample-ID column in expression data frames and phenotype.
#' @param type_col Phenotype column distinguishing tumors from cell lines.
#' @param lineage_col Phenotype column containing lineage annotations.
#' @param tumor_label Value in `type_col` identifying tumors.
#' @param cellline_label Value in `type_col` identifying cell lines.
#' @param input_scale Either `"raw"` (TPM) or `"processed"` (log2(TPM + 1)).
#' @param n_variable_genes Number of most variable genes used for PCA.
#' @param pca_variance_threshold Cumulative variance threshold used to select
#'   the retained number of PCs.
#' @param scale_genes Whether to center and scale genes before PCA.
#'
#' @return A `Tumor2LineReferencePCA` object containing the frozen PCA
#'   transformation, reference scores, reference metadata and cell-line scores.
#' @export
build_reference_pca <- function(
    tumor_expr,
    cellline_expr,
    phenotype,
    sample_id_col = "sampleID",
    type_col = "type",
    lineage_col = "lineage",
    tumor_label = "tumor",
    cellline_label = "CL",
    input_scale = c("processed", "raw"),
    n_variable_genes = 10000,
    pca_variance_threshold = 0.90,
    scale_genes = TRUE) {

  input_scale <- match.arg(input_scale)
  tumor_mat <- coerce_expression_matrix(tumor_expr, sample_id_col, "tumor_expr")
  cellline_mat <- coerce_expression_matrix(cellline_expr, sample_id_col, "cellline_expr")

  if (length(intersect(rownames(tumor_mat), rownames(cellline_mat))) > 0) {
    stop("Reference tumor and cell-line sample IDs overlap.", call. = FALSE)
  }
  common_genes <- intersect(colnames(tumor_mat), colnames(cellline_mat))
  if (length(common_genes) < 2) {
    stop("Fewer than two common genes were found between reference matrices.", call. = FALSE)
  }

  expression_input <- rbind(
    tumor_mat[, common_genes, drop = FALSE],
    cellline_mat[, common_genes, drop = FALSE]
  )
  normalized <- normalize_expression_scale(
    expression_input,
    input_scale = input_scale,
    object_name = "reference expression"
  )
  expression_log <- normalized$expression_log

  phenotype <- as.data.frame(phenotype, check.names = FALSE)
  assert_columns(phenotype, c(sample_id_col, type_col, lineage_col), "phenotype")
  phenotype[[sample_id_col]] <- as.character(phenotype[[sample_id_col]])
  if (anyNA(phenotype[[sample_id_col]]) ||
      any(phenotype[[sample_id_col]] == "") ||
      anyDuplicated(phenotype[[sample_id_col]])) {
    stop("Reference phenotype sample IDs must be non-missing, non-empty and unique.", call. = FALSE)
  }

  all_ids <- rownames(expression_log)
  phenotype_index <- match(all_ids, phenotype[[sample_id_col]])
  if (anyNA(phenotype_index)) {
    stop(
      "Reference phenotype is missing expression sample IDs: ",
      paste(utils::head(all_ids[is.na(phenotype_index)], 10), collapse = ", "),
      call. = FALSE
    )
  }

  phenotype_use <- phenotype[phenotype_index, , drop = FALSE]
  sample_type <- as.character(phenotype_use[[type_col]])
  lineage <- as.character(phenotype_use[[lineage_col]])
  is_tumor <- !is.na(sample_type) & sample_type == tumor_label
  is_cellline <- !is.na(sample_type) & sample_type == cellline_label
  unknown_type <- !(is_tumor | is_cellline)

  if (any(unknown_type)) {
    stop(
      "Reference phenotype contains unsupported sample type value(s): ",
      paste(unique(sample_type[unknown_type]), collapse = ", "),
      call. = FALSE
    )
  }
  if (!any(is_tumor) || !any(is_cellline)) {
    stop("Reference data must contain both tumor and cell-line samples.", call. = FALSE)
  }
  if (any(is.na(lineage[is_cellline]) | lineage[is_cellline] == "")) {
    stop("Reference cell-line samples have missing lineage annotations.", call. = FALSE)
  }

  n_variable_genes <- as.integer(n_variable_genes)
  if (length(n_variable_genes) != 1 || is.na(n_variable_genes) || n_variable_genes < 2) {
    stop("n_variable_genes must be an integer >= 2.", call. = FALSE)
  }
  gene_stats <- data.frame(
    gene = colnames(expression_log),
    log_sd = apply(expression_log, 2, stats::sd),
    stringsAsFactors = FALSE
  )
  gene_stats <- gene_stats[
    is.finite(gene_stats$log_sd) & gene_stats$log_sd > 0,
    ,
    drop = FALSE
  ]
  gene_stats <- gene_stats[order(-gene_stats$log_sd, gene_stats$gene), , drop = FALSE]
  if (nrow(gene_stats) < 2) {
    stop("Fewer than two variable genes are available in the reference data.", call. = FALSE)
  }
  selected_genes <- utils::head(
    gene_stats$gene,
    min(n_variable_genes, nrow(gene_stats))
  )

  pca_res <- fit_pca_representation(
    expression_log = expression_log,
    selected_genes = selected_genes,
    scale_genes = scale_genes
  )
  n_pcs <- select_pc_by_variance(pca_res, pca_variance_threshold)
  selected_pc_names <- colnames(pca_res$scores_raw)[seq_len(n_pcs)]

  representation <- list(
    selected_genes = selected_genes,
    n_pcs = n_pcs,
    variance = pca_res$variance,
    pca_rotation = pca_res$fit$rotation[, selected_pc_names, drop = FALSE],
    gene_center = pca_res$gene_center,
    gene_scale = pca_res$gene_scale,
    pc_center = pca_res$pc_center[seq_len(n_pcs)],
    pc_scale = pca_res$pc_scale[seq_len(n_pcs)],
    preprocessing = list(
      reference_dataset = "External reference",
      input_scale = input_scale,
      transformation = if (input_scale == "raw") "log2(TPM + 1)" else "none; input already log2(TPM + 1)",
      n_variable_genes = length(selected_genes),
      pca_variance_threshold = pca_variance_threshold,
      gene_scaling = scale_genes,
      pc_scaling = TRUE,
      pc_selection = "smallest number of PCs reaching the requested cumulative variance"
    )
  )

  scores <- pca_res$scores_scaled[, seq_len(n_pcs), drop = FALSE]
  colnames(scores) <- selected_pc_names

  metadata <- data.frame(
    SampleID = all_ids,
    SampleRole = ifelse(is_cellline, "cellline", "reference_tumor"),
    SampleType = sample_type,
    lineage = lineage,
    PairID = NA_character_,
    stringsAsFactors = FALSE
  )

  cellline_index <- which(is_cellline)
  tumor_index <- which(is_tumor)
  out <- list(
    representation = representation,
    metadata = metadata,
    cellline_scores = scores[cellline_index, , drop = FALSE],
    cellline_metadata = metadata[cellline_index, , drop = FALSE],
    tumor_scores = scores[tumor_index, , drop = FALSE],
    tumor_metadata = metadata[tumor_index, , drop = FALSE],
    parameters = list(
      sample_id_col = sample_id_col,
      type_col = type_col,
      lineage_col = lineage_col,
      tumor_label = tumor_label,
      cellline_label = cellline_label,
      input_scale = input_scale,
      n_variable_genes = length(selected_genes),
      pca_variance_threshold = pca_variance_threshold,
      selected_pcs = n_pcs,
      common_gene_count = length(common_genes),
      n_reference_samples = nrow(expression_log),
      n_reference_tumors = sum(is_tumor),
      n_reference_celllines = sum(is_cellline)
    )
  )
  class(out) <- c("Tumor2LineReferencePCA", "list")
  out
}

#' Project expression profiles into a frozen reference PCA space
#'
#' Applies the reference gene set, gene-wise scaling, PCA loadings and reference
#' PC-score scaling without refitting any transformation.
#'
#' @param expression Expression data with samples in rows and genes in columns.
#' @param reference A `Tumor2LineReferencePCA` object.
#' @param input_scale Either `"raw"` (TPM) or `"processed"` (log2(TPM + 1)).
#' @param sample_id_col Sample-ID column in data-frame expression input.
#' @return A list containing projected scaled and raw PC scores, normalized
#'   log-expression and the declared input scale.
#' @export
project_to_reference_pca <- function(
    expression,
    reference,
    input_scale = c("raw", "processed"),
    sample_id_col = "SampleID") {

  if (!inherits(reference, "Tumor2LineReferencePCA")) {
    stop("reference must be a Tumor2LineReferencePCA object.", call. = FALSE)
  }
  input_scale <- match.arg(input_scale)
  expression_mat <- coerce_expression_matrix(expression, sample_id_col, "expression")
  normalized <- normalize_expression_scale(expression_mat, input_scale, "expression")
  projected <- project_expression_to_reference(
    expression_log = normalized$expression_log,
    reference_representation = reference$representation
  )

  list(
    scores = projected$scores,
    scores_raw = projected$scores_raw,
    expression_log = normalized$expression_log,
    input_scale = input_scale
  )
}

#' Apply an explicit source-to-reference lineage crosswalk
#'
#' Adds a harmonized lineage column using a user-supplied mapping. This avoids
#' assuming that source-specific lineage codes are interchangeable with the
#' descriptive lineage labels used by an external reference dataset.
#'
#' @param phenotype Data frame containing source lineage annotations.
#' @param crosswalk Data frame with one row per source-to-reference mapping.
#' @param source_lineage_col Source lineage column in `phenotype`.
#' @param crosswalk_source_col Source lineage column in `crosswalk`.
#' @param crosswalk_target_col Target reference lineage column in `crosswalk`.
#' @param target_col Name of the output harmonized lineage column.
#' @return `phenotype` with an additional harmonized lineage column.
#' @export
apply_lineage_crosswalk <- function(
    phenotype,
    crosswalk,
    source_lineage_col = "TCGA_lineage",
    crosswalk_source_col = "TCGA_lineage",
    crosswalk_target_col = "Celligner_lineage",
    target_col = "Celligner_lineage") {

  phenotype <- as.data.frame(phenotype, check.names = FALSE)
  crosswalk <- as.data.frame(crosswalk, check.names = FALSE)
  assert_columns(phenotype, source_lineage_col, "phenotype")
  assert_columns(crosswalk, c(crosswalk_source_col, crosswalk_target_col), "crosswalk")

  source <- as.character(crosswalk[[crosswalk_source_col]])
  target <- as.character(crosswalk[[crosswalk_target_col]])
  if (anyNA(source) || any(source == "") || anyDuplicated(source)) {
    stop("Crosswalk source lineage values must be non-missing, non-empty and unique.", call. = FALSE)
  }
  if (anyNA(target) || any(target == "")) {
    stop("Crosswalk target lineage values must be non-missing and non-empty.", call. = FALSE)
  }

  phenotype_source <- as.character(phenotype[[source_lineage_col]])
  index <- match(phenotype_source, source)
  unmapped <- unique(phenotype_source[
    !is.na(phenotype_source) & phenotype_source != "" & is.na(index)
  ])
  if (length(unmapped) > 0) {
    stop(
      "The lineage crosswalk does not cover: ",
      paste(unmapped, collapse = ", "),
      call. = FALSE
    )
  }
  phenotype[[target_col]] <- target[index]
  phenotype
}

#' Run Tumor2Line matching in a frozen external reference PCA space
#'
#' Projects query samples and candidate cell lines into a frozen external
#' reference representation. Expression genes are harmonized across the
#' reference PCA, query and candidate cell-line datasets before projection.
#'
#' Missing reference genes are treated as reference-mean expression values,
#' corresponding to zero after reference gene-wise standardization. The
#' original PCA axes and PC scaling parameters are retained.
#'
#' @param query_expr Query expression matrix.
#' @param reference Frozen reference PCA object from [build_reference_pca()].
#' @param query_phenotype Query phenotype data frame.
#' @param cellline_expr Optional candidate cell-line expression matrix. If
#'   `NULL`, reference cell-line scores are used.
#' @param cellline_phenotype Optional phenotype for `cellline_expr`.
#' @param query_sample_id_col Sample-ID column in query data and phenotype.
#' @param cellline_sample_id_col Sample-ID column in candidate data and phenotype.
#' @param query_lineage_col Harmonized lineage column in query phenotype.
#' @param cellline_lineage_col Lineage column in candidate phenotype.
#' @param pair_col Optional known-pair column present in both phenotypes.
#' @param query_type Biological label assigned to query samples.
#' @param cellline_type Biological label assigned to candidate cell lines.
#' @param query_input_scale Input scale of query expression: `"raw"` or `"processed"`.
#' @param cellline_input_scale Input scale of candidate expression.
#' @param candidate_scope `"within_lineage"` or `"global"`.
#' @param distance Distance metric for matching.
#' @param epsilon Entropic regularization parameter.
#' @param rho_query Query marginal-relaxation parameter.
#' @param rho_cellline Cell-line marginal-relaxation parameter.
#' @param max_iter Maximum Sinkhorn iterations.
#' @param tolerance Sinkhorn convergence tolerance.
#' @param evaluate Whether to evaluate known pairs when supplied.
#'
#' @return A `Tumor2LineResult` object.
#' @export
run_tumor2line_with_reference <- function(
    query_expr,
    reference,
    query_phenotype,
    cellline_expr = NULL,
    cellline_phenotype = NULL,
    query_sample_id_col = "SampleID",
    cellline_sample_id_col = "SampleID",
    query_lineage_col = "lineage",
    cellline_lineage_col = "lineage",
    pair_col = NULL,
    query_type = "tumor",
    cellline_type = "cellline",
    query_input_scale = c("raw", "processed"),
    cellline_input_scale = c("raw", "processed"),
    candidate_scope = c("within_lineage", "global"),
    distance = c("euclidean", "manhattan", "cosine"),
    epsilon = 0.1,
    rho_query = 1,
    rho_cellline = 1,
    max_iter = 5000,
    tolerance = 1e-8,
    evaluate = TRUE) {
  
  if (!inherits(reference, "Tumor2LineReferencePCA")) {
    stop(
      "reference must be created by build_reference_pca().",
      call. = FALSE
    )
  }
  
  candidate_scope <- match.arg(candidate_scope)
  distance <- match.arg(distance)
  query_input_scale <- match.arg(query_input_scale)
  cellline_input_scale <- match.arg(cellline_input_scale)
  
  query_mat <- coerce_expression_matrix(
    query_expr,
    query_sample_id_col,
    "query_expr"
  )
  
  query_phenotype <- as.data.frame(
    query_phenotype,
    check.names = FALSE
  )
  
  assert_columns(
    query_phenotype,
    c(query_sample_id_col, query_lineage_col),
    "query_phenotype"
  )
  
  query_phenotype[[query_sample_id_col]] <- as.character(
    query_phenotype[[query_sample_id_col]]
  )
  
  if (anyNA(query_phenotype[[query_sample_id_col]]) ||
      any(query_phenotype[[query_sample_id_col]] == "") ||
      anyDuplicated(query_phenotype[[query_sample_id_col]])) {
    stop(
      "query_phenotype sample IDs must be non-missing, non-empty and unique.",
      call. = FALSE
    )
  }
  
  # Match query expression samples with phenotype IDs
  query_keep <- rownames(query_mat) %in% query_phenotype[[query_sample_id_col]]
  
  if (any(!query_keep)) {
    warning(
      sum(!query_keep),
      " query sample(s) were removed because they have no matching phenotype.",
      call. = FALSE
    )
  }
  
  query_mat <- query_mat[query_keep, , drop = FALSE]
  
  if (nrow(query_mat) == 0) {
    stop(
      "No query samples remain after matching expression and phenotype IDs.",
      call. = FALSE
    )
  }
  
  query_ids <- rownames(query_mat)
  
  query_idx <- match(
    query_ids,
    query_phenotype[[query_sample_id_col]]
  )
  
  # Read and validate candidate cell-line expression ------------------------
  
  cellline_mat <- NULL
  
  if (!is.null(cellline_expr)) {
    
    cellline_mat <- coerce_expression_matrix(
      cellline_expr,
      cellline_sample_id_col,
      "cellline_expr"
    )
    
    if (is.null(cellline_phenotype)) {
      cellline_phenotype <- query_phenotype
    }
    
    cellline_phenotype <- as.data.frame(
      cellline_phenotype,
      check.names = FALSE
    )
    
    assert_columns(
      cellline_phenotype,
      c(cellline_sample_id_col, cellline_lineage_col),
      "cellline_phenotype"
    )
    
    cellline_phenotype[[cellline_sample_id_col]] <- as.character(
      cellline_phenotype[[cellline_sample_id_col]]
    )
    
    if (anyNA(cellline_phenotype[[cellline_sample_id_col]]) ||
        any(cellline_phenotype[[cellline_sample_id_col]] == "") ||
        anyDuplicated(cellline_phenotype[[cellline_sample_id_col]])) {
      stop(
        "cellline_phenotype sample IDs must be non-missing, non-empty and unique.",
        call. = FALSE
      )
    }
    
    # Match cell-line expression samples with phenotype IDs
    cellline_keep <- rownames(cellline_mat) %in%
      cellline_phenotype[[cellline_sample_id_col]]
    
    if (any(!cellline_keep)) {
      warning(
        sum(!cellline_keep),
        " cell-line sample(s) were removed because they have no matching phenotype.",
        call. = FALSE
      )
    }
    
    cellline_mat <- cellline_mat[cellline_keep, , drop = FALSE]
    
    if (nrow(cellline_mat) == 0) {
      stop(
        "No cell-line samples remain after matching expression and phenotype IDs.",
        call. = FALSE
      )
    }
    
    cellline_ids <- rownames(cellline_mat)
    
    cellline_idx <- match(
      cellline_ids,
      cellline_phenotype[[cellline_sample_id_col]]
    )
    
  } else {
    
    # Use the cell-line scores already stored in the reference
    cellline_scores <- as.matrix(reference$cellline_scores)
    
    if (is.null(rownames(cellline_scores))) {
      stop(
        "Reference cell-line scores must have sample IDs as row names.",
        call. = FALSE
      )
    }
    
    cellline_ids <- rownames(cellline_scores)
    
    reference_cellline_metadata <- as.data.frame(
      reference$cellline_metadata,
      check.names = FALSE
    )
    
    if ("SampleID" %in% colnames(reference_cellline_metadata)) {
      
      reference_cellline_idx <- match(
        cellline_ids,
        as.character(reference_cellline_metadata$SampleID)
      )
      
    } else if (all(cellline_ids %in% rownames(reference_cellline_metadata))) {
      
      reference_cellline_idx <- match(
        cellline_ids,
        rownames(reference_cellline_metadata)
      )
      
    } else if (nrow(reference_cellline_metadata) == length(cellline_ids)) {
      
      reference_cellline_idx <- seq_along(cellline_ids)
      
    } else {
      
      stop(
        "Reference cell-line metadata cannot be aligned with cell-line scores.",
        call. = FALSE
      )
    }
    
    if (anyNA(reference_cellline_idx)) {
      stop(
        "Some reference cell lines are missing from reference metadata.",
        call. = FALSE
      )
    }
    
    cellline_meta <- data.frame(
      SampleID = cellline_ids,
      SampleRole = "cellline",
      SampleType = cellline_type,
      lineage = if ("lineage" %in% colnames(reference_cellline_metadata)) {
        as.character(
          reference_cellline_metadata$lineage[reference_cellline_idx]
        )
      } else {
        rep(NA_character_, length(cellline_ids))
      },
      PairID = if ("PairID" %in% colnames(reference_cellline_metadata)) {
        as.character(
          reference_cellline_metadata$PairID[reference_cellline_idx]
        )
      } else {
        rep(NA_character_, length(cellline_ids))
      },
      stringsAsFactors = FALSE
    )
  }
  
  reference_genes <- as.character(
    reference$representation$selected_genes
  )
  
  if (length(reference_genes) == 0 ||
      anyNA(reference_genes) ||
      any(reference_genes == "") ||
      anyDuplicated(reference_genes)) {
    stop(
      "Reference PCA must contain valid, unique selected genes.",
      call. = FALSE
    )
  }
  
  # Genes available in both query and candidate cell-line data
  available_genes <- colnames(query_mat)
  
  if (!is.null(cellline_mat)) {
    available_genes <- intersect(
      available_genes,
      colnames(cellline_mat)
    )
  }
  
  # Preserve the original reference PCA gene order
  common_genes <- reference_genes[
    reference_genes %in% available_genes
  ]
  
  if (length(common_genes) < 2) {
    stop(
      "Fewer than two common genes are available between the input data and reference PCA.",
      call. = FALSE
    )
  }
  
  missing_reference_genes <- setdiff(
    reference_genes,
    common_genes
  )
  
  if (length(missing_reference_genes) > 0) {
    warning(
      length(missing_reference_genes),
      " of ",
      length(reference_genes),
      " reference PCA genes are unavailable in the input data (",
      round(
        100 * length(missing_reference_genes) / length(reference_genes),
        2
      ),
      "%). Missing genes will contribute zero after reference gene-wise standardization.",
      call. = FALSE
    )
  }
  
  # Subset expression matrices to the common genes
  query_mat <- query_mat[, common_genes, drop = FALSE]
  
  if (!is.null(cellline_mat)) {
    cellline_mat <- cellline_mat[, common_genes, drop = FALSE]
  }
  
  # Subset the frozen PCA projection parameters -----------------------------
  
  projection_reference <- reference
  
  pca_rotation <- as.matrix(
    reference$representation$pca_rotation
  )
  
  if (is.null(rownames(pca_rotation))) {
    
    if (nrow(pca_rotation) != length(reference_genes)) {
      stop(
        "PCA rotation rows cannot be aligned with reference selected genes.",
        call. = FALSE
      )
    }
    
    rownames(pca_rotation) <- reference_genes
  }
  
  if (anyDuplicated(rownames(pca_rotation)) ||
      !all(common_genes %in% rownames(pca_rotation))) {
    stop(
      "Reference PCA rotation cannot be aligned with common genes.",
      call. = FALSE
    )
  }
  
  gene_center <- reference$representation$gene_center
  gene_scale <- reference$representation$gene_scale
  
  if (is.null(names(gene_center))) {
    if (length(gene_center) != length(reference_genes)) {
      stop("Reference gene_center cannot be aligned with selected genes.")
    }
    names(gene_center) <- reference_genes
  }
  
  if (is.null(names(gene_scale))) {
    if (length(gene_scale) != length(reference_genes)) {
      stop("Reference gene_scale cannot be aligned with selected genes.")
    }
    names(gene_scale) <- reference_genes
  }
  
  if (!all(common_genes %in% names(gene_center)) ||
      !all(common_genes %in% names(gene_scale))) {
    stop(
      "Reference gene scaling parameters are missing common genes.",
      call. = FALSE
    )
  }
  
  if (any(!is.finite(gene_center[common_genes])) ||
      any(!is.finite(gene_scale[common_genes])) ||
      any(gene_scale[common_genes] <= 0)) {
    stop(
      "Reference gene scaling parameters contain invalid values.",
      call. = FALSE
    )
  }
  
  projection_reference$representation$selected_genes <- common_genes
  
  projection_reference$representation$gene_center <-
    gene_center[common_genes]
  
  projection_reference$representation$gene_scale <-
    gene_scale[common_genes]
  
  projection_reference$representation$pca_rotation <-
    pca_rotation[common_genes, , drop = FALSE]
  
  query_projection <- project_to_reference_pca(
    query_mat,
    projection_reference,
    input_scale = query_input_scale,
    sample_id_col = query_sample_id_col
  )
  
  query_ids <- rownames(query_projection$scores)
  
  query_meta <- data.frame(
    SampleID = query_ids,
    SampleRole = "query",
    SampleType = query_type,
    lineage = if (is.null(query_lineage_col)) {
      rep(NA_character_, length(query_ids))
    } else {
      as.character(query_phenotype[[query_lineage_col]][query_idx])
    },
    PairID = NA_character_,
    stringsAsFactors = FALSE
  )
  
  if (!is.null(pair_col)) {
    assert_columns(
      query_phenotype,
      pair_col,
      "query_phenotype"
    )
    
    query_meta$PairID <- as.character(
      query_phenotype[[pair_col]][query_idx]
    )
  }
  
  if (!is.null(cellline_mat)) {
    
    cellline_projection <- project_to_reference_pca(
      cellline_mat,
      projection_reference,
      input_scale = cellline_input_scale,
      sample_id_col = cellline_sample_id_col
    )
    
    cellline_scores <- cellline_projection$scores
    cellline_ids <- rownames(cellline_scores)
    
    cellline_meta <- data.frame(
      SampleID = cellline_ids,
      SampleRole = "cellline",
      SampleType = cellline_type,
      lineage = if (is.null(cellline_lineage_col)) {
        rep(NA_character_, length(cellline_ids))
      } else {
        as.character(
          cellline_phenotype[[cellline_lineage_col]][cellline_idx]
        )
      },
      PairID = NA_character_,
      stringsAsFactors = FALSE
    )
    
    if (!is.null(pair_col) &&
        pair_col %in% colnames(cellline_phenotype)) {
      cellline_meta$PairID <- as.character(
        cellline_phenotype[[pair_col]][cellline_idx]
      )
    }
  }
  
  if (length(intersect(query_ids, cellline_ids)) > 0) {
    stop(
      "Query and candidate cell-line sample IDs overlap.",
      call. = FALSE
    )
  }
  
  if (candidate_scope == "within_lineage") {
    
    if (any(is.na(query_meta$lineage) | query_meta$lineage == "")) {
      stop(
        "Query lineage annotations are required for within-lineage matching.",
        call. = FALSE
      )
    }
    
    if (any(is.na(cellline_meta$lineage) | cellline_meta$lineage == "")) {
      stop(
        "Candidate cell-line lineage annotations are required for within-lineage matching.",
        call. = FALSE
      )
    }
  }
  
  all_scores <- rbind(
    query_projection$scores,
    cellline_scores
  )
  
  all_meta <- rbind(
    query_meta,
    cellline_meta
  )
  
  ranking <- run_uot_matching(
    representation = all_scores,
    metadata = all_meta,
    candidate_scope = candidate_scope,
    distance = distance,
    epsilon = epsilon,
    rho_query = rho_query,
    rho_cellline = rho_cellline,
    max_iter = max_iter,
    tolerance = tolerance
  )
  
  top_matches <- ranking |>
    dplyr::filter(.data$rank == 1L) |>
    dplyr::arrange(.data$lineage, .data$Query_ID)
  
  gold_pairs <- data.frame(
    PairID = character(),
    lineage = character(),
    Query_ID = character(),
    CellLine_ID = character(),
    stringsAsFactors = FALSE
  )
  
  if (!is.null(pair_col)) {
    
    query_pairs <- query_meta[
      !is.na(query_meta$PairID) & query_meta$PairID != "",
      c("PairID", "SampleID", "lineage"),
      drop = FALSE
    ]
    
    cellline_pairs <- cellline_meta[
      !is.na(cellline_meta$PairID) & cellline_meta$PairID != "",
      c("PairID", "SampleID", "lineage"),
      drop = FALSE
    ]
    
    if (nrow(query_pairs) > 0 && nrow(cellline_pairs) > 0) {
      
      names(query_pairs) <- c(
        "PairID", "Query_ID", "lineage"
      )
      
      names(cellline_pairs) <- c(
        "PairID", "CellLine_ID", "CellLine_lineage"
      )
      
      pair_match <- merge(
        query_pairs,
        cellline_pairs,
        by = "PairID",
        all = FALSE
      )
      
      if (nrow(pair_match) > 0) {
        gold_pairs <- unique(data.frame(
          PairID = pair_match$PairID,
          lineage = pair_match$lineage,
          Query_ID = pair_match$Query_ID,
          CellLine_ID = pair_match$CellLine_ID,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  cohort_summary <- tibble::as_tibble(all_meta) |>
    dplyr::mutate(
      lineage = dplyr::if_else(
        is.na(.data$lineage) | .data$lineage == "",
        "Unspecified",
        .data$lineage
      )
    ) |>
    dplyr::group_by(.data$lineage) |>
    dplyr::summarise(
      Query_number = sum(.data$SampleRole == "query"),
      Cell_line_number = sum(.data$SampleRole == "cellline"),
      Total_samples = dplyr::n(),
      .groups = "drop"
    )
  
  input <- list(
    metadata = all_meta,
    gold_pairs = gold_pairs,
    query_ids = query_ids,
    cellline_ids = cellline_ids,
    cohort_summary = cohort_summary,
    input_scales = list(
      query = query_input_scale,
      cellline = if (is.null(cellline_expr)) {
        "reference_scores"
      } else {
        cellline_input_scale
      }
    )
  )
  
  representation <- list(
    reference_parameters = reference$parameters,
    selected_genes = common_genes,
    n_pcs = reference$representation$n_pcs,
    scores = all_scores,
    scores_raw = NULL,
    pca_rotation = projection_reference$representation$pca_rotation,
    gene_center = projection_reference$representation$gene_center,
    gene_scale = projection_reference$representation$gene_scale,
    pc_center = reference$representation$pc_center,
    pc_scale = reference$representation$pc_scale,
    variance = reference$representation$variance
  )
  
  result <- list(
    ranking = ranking,
    top_matches = top_matches,
    input = input,
    representation = representation,
    parameters = list(
      pca_reference = "External reference PCA",
      pca_variance_threshold = reference$parameters$pca_variance_threshold,
      selected_pcs = reference$representation$n_pcs,
      reference_gene_count = length(reference_genes),
      projection_gene_count = length(common_genes),
      omitted_reference_gene_count = length(missing_reference_genes),
      projection_gene_coverage = length(common_genes) / length(reference_genes),
      candidate_scope = candidate_scope,
      distance = distance,
      epsilon = epsilon,
      rho_query = rho_query,
      rho_cellline = rho_cellline,
      max_iter = max_iter,
      tolerance = tolerance,
      query_input_scale = query_input_scale,
      cellline_input_scale = if (is.null(cellline_expr)) {
        "reference_scores"
      } else {
        cellline_input_scale
      }
    )
  )
  
  class(result) <- c("Tumor2LineResult", "list")
  
  if (isTRUE(evaluate) && nrow(gold_pairs) > 0) {
    result$evaluation <- evaluate_tumor2line(
      result,
      gold_pairs
    )
  } else {
    result$evaluation <- NULL
  }
  
  result
}
