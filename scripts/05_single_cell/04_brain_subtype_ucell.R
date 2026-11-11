source("R/runtime.R")

suppressPackageStartupMessages({
  library(BiocParallel)
  library(data.table)
  library(hdf5r)
  library(Matrix)
  library(UCell)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !args[[1]] %chin% c("broad", "subtype")) {
  stop("Usage: Rscript <script> broad|subtype")
}

MODE <- args[[1]]
ROOT <- "."
WORK <- file.path(ROOT, "results/single_cell/results/brain_work")
MODE_DIR <- file.path(WORK, paste0(MODE, "_ucell"))
dir.create(MODE_DIR, recursive = TRUE, showWarnings = FALSE)

if (MODE == "broad") {
  H5AD <- file.path(
    ROOT, "data", "raw", "single_cell", "sea_ad_pfc_a9",
    "SEAAD_DFC_RNAseq_final-nuclei.2026-06-22.h5ad"
  )
  META_FILE <- file.path(WORK, "broad_metadata.tsv.gz")
  MATRIX_GROUP <- "X"
  DATASET <- "SEA-AD_DFC_final-nuclei_2026-06-22"
} else {
  H5AD <- file.path(
    ROOT, "data", "raw", "single_cell", "sea_ad_pfc_a9",
    "SEAAD_Immune_HIP_MEC_LEC_ITG_MTG_FI_STG_DFC_AnG_V1C_RNAseq_final-nuclei.2026-06-22.h5ad"
  )
  META_FILE <- file.path(WORK, "subtype_metadata.tsv.gz")
  MATRIX_GROUP <- "X"
  DATASET <- "SEA-AD_Immune_final-nuclei_2026-06-22"
}

GENE_FILE <- file.path(WORK, "genes.txt")
AXIS_FILE <- file.path(WORK, "axis_measured.txt")
COMMON_FILE <- file.path(WORK, "axis_measured.txt")
UNIVERSE_FILE <- file.path(WORK, "universe_measured.txt")
MARKER_FILE <- file.path(WORK, "independent_markers.tsv")

RMAX <- params$ucell_max_rank
OUTER_CHUNK <- 5000L
UCELL_CHUNK <- 200L
NCORES <- 4L

CELL_SCORE_FILE <- file.path(MODE_DIR, "cell_scores.tsv")
RANK_AGG_FILE <- file.path(MODE_DIR, "profile_gene_rank_aggregate.rds")
GENE_PROPERTY_FILE <- file.path(MODE_DIR, "gene_properties.tsv")
MARKER_PROPERTY_FILE <- file.path(MODE_DIR, "marker_profile_properties.tsv")
VALIDATION_FILE <- file.path(MODE_DIR, "ucell_exact_validation.tsv")
SESSION_FILE <- file.path(MODE_DIR, "sessionInfo.txt")
COMPLETION_FILE <- file.path(MODE_DIR, "COMPLETED.txt")
CHECKPOINT_FILE <- file.path(MODE_DIR, "stream_checkpoint.rds")
CHECKPOINT_EVERY_OBJECT_ROWS <- 100000L

if (file.exists(COMPLETION_FILE)) {
  stop("Streamed UCell run is already complete: ", COMPLETION_FILE)
}
resuming <- file.exists(CHECKPOINT_FILE)
outputs <- c(RANK_AGG_FILE, GENE_PROPERTY_FILE, MARKER_PROPERTY_FILE, SESSION_FILE)
if (any(file.exists(outputs))) {
  stop(
    "Refusing to overwrite an existing final streamed-output file: ",
    paste(outputs[file.exists(outputs)], collapse = "; ")
  )
}
if (!resuming && any(file.exists(c(CELL_SCORE_FILE, VALIDATION_FILE)))) {
  stop("Partial score/validation files exist without a resumable checkpoint")
}

genes <- readLines(GENE_FILE)
axis_genes <- readLines(AXIS_FILE)
common_genes <- readLines(COMMON_FILE)
universe_genes <- readLines(UNIVERSE_FILE)
marker_dt <- fread(MARKER_FILE)[measured == TRUE]
marker_genes <- unique(marker_dt$gene)
meta <- fread(META_FILE)
setorder(meta, row_index)

stopifnot(
  length(axis_genes) == 136L,
  length(common_genes) == 136L,
  uniqueN(axis_genes) == 136L,
  uniqueN(universe_genes) == length(universe_genes),
  all(axis_genes %chin% genes),
  all(common_genes %chin% genes),
  all(universe_genes %chin% genes),
  all(marker_genes %chin% genes),
  !any(marker_genes %chin% axis_genes)
)

profile_dt <- unique(meta[, .(donor, annotation_level, annotation)])
setorder(profile_dt, donor, annotation)
profile_dt[, profile_index := .I]
meta <- merge(
  meta, profile_dt,
  by = c("donor", "annotation_level", "annotation"),
  all.x = TRUE,
  sort = FALSE
)
setorder(meta, row_index)

h5 <- H5File$new(H5AD, mode = "r")
on.exit(h5$close_all(), add = TRUE)
group <- h5[[MATRIX_GROUP]]
shape <- as.integer(group$attr_open("shape")$read())
n_object_cells <- shape[[1]]
n_genes <- shape[[2]]
if (n_genes != length(genes)) stop("H5AD feature dimension drift")
if (max(meta$row_index) >= n_object_cells) stop("Metadata row index exceeds H5AD")

profile_by_object_row <- integer(n_object_cells)
meta_by_object_row <- integer(n_object_cells)
profile_by_object_row[meta$row_index + 1L] <- meta$profile_index
meta_by_object_row[meta$row_index + 1L] <- seq_len(nrow(meta))

universe_idx <- match(universe_genes, genes)
marker_idx <- match(marker_genes, genes)
n_profiles <- nrow(profile_dt)
n_universe <- length(universe_genes)
n_markers <- length(marker_genes)

rank_sum_nonzero <- matrix(0, nrow = n_universe, ncol = n_profiles)
rank_count_nonzero <- matrix(0, nrow = n_universe, ncol = n_profiles)
profile_n_cells <- numeric(n_profiles)
gene_sum <- numeric(n_universe)
gene_sum_sq <- numeric(n_universe)
gene_detected <- numeric(n_universe)
marker_sum <- matrix(0, nrow = n_markers, ncol = n_profiles)
marker_detected <- matrix(0, nrow = n_markers, ncol = n_profiles)

indptr <- as.numeric(group[["indptr"]][])
indices_ds <- group[["indices"]]
data_ds <- group[["data"]]
if (length(indptr) != n_object_cells + 1L) stop("CSR indptr length drift")

signatures <- list(
  AxisA_full_measured = axis_genes,
  AxisA_common129 = common_genes
)

score_from_stored_rankings <- function(rankings, signatures, max_rank) {
  out <- vapply(
    signatures,
    function(signature) {
      signature_idx <- match(signature, rownames(rankings))
      if (anyNA(signature_idx)) stop("Measured signature gene missing from stored rankings")
      signature_ranks <- as.matrix(rankings[signature_idx, , drop = FALSE])
      signature_ranks[signature_ranks == 0] <- max_rank
      signature_size <- length(signature_idx)
      rank_sum_min <- signature_size * (signature_size + 1) / 2
      1 - (colSums(signature_ranks) - rank_sum_min) /
        (signature_size * max_rank - rank_sum_min)
    },
    numeric(ncol(rankings))
  )
  colnames(out) <- paste0(names(signatures), "_UCell")
  rownames(out) <- colnames(rankings)
  out
}

validation_done <- FALSE
started <- Sys.time()
selected_processed <- 0L
resume_row_start <- 1L

if (resuming) {
  checkpoint <- readRDS(CHECKPOINT_FILE)
  stopifnot(
    identical(checkpoint$mode, MODE),
    checkpoint$n_object_cells == n_object_cells,
    checkpoint$n_profiles == n_profiles,
    checkpoint$n_universe == n_universe,
    file.exists(CELL_SCORE_FILE)
  )
  rank_sum_nonzero <- checkpoint$rank_sum_nonzero
  rank_count_nonzero <- checkpoint$rank_count_nonzero
  profile_n_cells <- checkpoint$profile_n_cells
  gene_sum <- checkpoint$gene_sum
  gene_sum_sq <- checkpoint$gene_sum_sq
  gene_detected <- checkpoint$gene_detected
  marker_sum <- checkpoint$marker_sum
  marker_detected <- checkpoint$marker_detected
  validation_done <- checkpoint$validation_done
  started <- checkpoint$started
  selected_processed <- checkpoint$selected_processed
  resume_row_start <- checkpoint$next_row_start

  prior_scores <- fread(CELL_SCORE_FILE)
  if (nrow(prior_scores) < selected_processed) {
    stop("Cell-score file has fewer rows than the committed checkpoint")
  }
  if (nrow(prior_scores) > selected_processed) {
    prior_scores <- prior_scores[seq_len(selected_processed)]
    score_tmp <- paste0(CELL_SCORE_FILE, ".resume_tmp")
    fwrite(prior_scores, score_tmp, sep = "\t")
    if (!file.rename(score_tmp, CELL_SCORE_FILE)) stop("Could not trim score file to checkpoint")
  }
  rm(prior_scores, checkpoint)
  message(sprintf(
    "%s: resuming at object row %d with %d selected cells committed",
    MODE, resume_row_start, selected_processed
  ))
}

for (row_start in seq.int(resume_row_start, n_object_cells, by = OUTER_CHUNK)) {
  row_end <- min(n_object_cells, row_start + OUTER_CHUNK - 1L)
  selected_local <- which(profile_by_object_row[row_start:row_end] > 0L)
  if (!length(selected_local)) next

  ptr <- indptr[row_start:(row_end + 1L)]
  nnz_start <- ptr[[1]] + 1
  nnz_end <- ptr[[length(ptr)]]
  if (nnz_end < nnz_start) stop("Unexpected empty CSR chunk")
  indices <- as.integer(indices_ds[nnz_start:nnz_end]) + 1L
  values <- as.numeric(data_ds[nnz_start:nnz_end])
  counts_per_row <- diff(ptr)
  cell_j <- rep.int(seq_along(counts_per_row), counts_per_row)
  matrix_chunk <- sparseMatrix(
    i = indices,
    j = cell_j,
    x = values,
    dims = c(n_genes, length(counts_per_row)),
    dimnames = list(genes, as.character((row_start - 1L):(row_end - 1L)))
  )
  if (length(selected_local) != ncol(matrix_chunk)) {
    matrix_chunk <- matrix_chunk[, selected_local, drop = FALSE]
  }
  global_rows <- (row_start - 1L) + selected_local
  profile_idx <- profile_by_object_row[global_rows]
  meta_idx <- meta_by_object_row[global_rows]
  if (any(profile_idx <= 0L) || any(meta_idx <= 0L)) stop("Selected-row mapping failure")

  rankings <- StoreRankings_UCell(
    matrix = matrix_chunk,
    maxRank = RMAX,
    chunk.size = UCELL_CHUNK,
    BPPARAM = MulticoreParam(workers = NCORES),
    ncores = NCORES,
    ties.method = "average",
    force.gc = TRUE
  )
  scores <- score_from_stored_rankings(rankings, signatures, RMAX)

  if (!validation_done) {
    n_validate <- min(200L, ncol(matrix_chunk))
    direct <- ScoreSignatures_UCell(
      matrix = matrix_chunk[, seq_len(n_validate), drop = FALSE],
      features = signatures,
      maxRank = RMAX,
      missing_genes = "skip",
      BPPARAM = SerialParam(),
      ncores = 1L,
      chunk.size = 100L,
      ties.method = "average"
    )
    delta <- abs(
      direct[, "AxisA_full_measured_UCell"] -
        scores[seq_len(n_validate), "AxisA_full_measured_UCell"]
    )
    validation <- data.table(
      dataset = DATASET,
      matrix = MATRIX_GROUP,
      UCell_version = as.character(packageVersion("UCell")),
      rmax = RMAX,
      n_validation_cells = n_validate,
      max_abs_delta_direct_vs_precalculated = max(delta),
      exact = max(delta) == 0
    )
    fwrite(validation, VALIDATION_FILE, sep = "\t")
    if (!validation$exact) stop("Formal UCell direct/precalculated validation failed")
    validation_done <- TRUE
  }

  score_dt <- data.table(
    dataset = DATASET,
    region = "DFC",
    annotation_level = meta$annotation_level[meta_idx],
    row_index = global_rows - 1L,
    cell_id = meta$cell_id[meta_idx],
    donor = meta$donor[meta_idx],
    annotation = meta$annotation[meta_idx],
    axisA_full_measured_ucell = as.numeric(scores[, "AxisA_full_measured_UCell"]),
    axisA_common129_ucell = as.numeric(scores[, "AxisA_common129_UCell"])
  )
  fwrite(
    score_dt,
    CELL_SCORE_FILE,
    sep = "\t",
    append = file.exists(CELL_SCORE_FILE),
    col.names = !file.exists(CELL_SCORE_FILE)
  )

  active_profiles <- sort(unique(profile_idx))
  membership <- sparseMatrix(
    i = seq_along(profile_idx),
    j = match(profile_idx, active_profiles),
    x = 1,
    dims = c(length(profile_idx), length(active_profiles))
  )
  profile_n_cells[active_profiles] <- profile_n_cells[active_profiles] +
    as.numeric(colSums(membership))

  rank_sub <- rankings[universe_idx, , drop = FALSE]
  rank_sum_nonzero[, active_profiles] <- rank_sum_nonzero[, active_profiles, drop = FALSE] +
    as.matrix(rank_sub %*% membership)
  rank_count_nonzero[, active_profiles] <- rank_count_nonzero[, active_profiles, drop = FALSE] +
    as.matrix((rank_sub != 0) %*% membership)

  x_universe <- matrix_chunk[universe_idx, , drop = FALSE]
  gene_sum <- gene_sum + as.numeric(rowSums(x_universe))
  gene_sum_sq <- gene_sum_sq + as.numeric(rowSums(x_universe^2))
  gene_detected <- gene_detected + as.numeric(rowSums(x_universe != 0))

  if (MODE == "subtype" && n_markers > 0L) {
    x_marker <- matrix_chunk[marker_idx, , drop = FALSE]
    marker_sum[, active_profiles] <- marker_sum[, active_profiles, drop = FALSE] +
      as.matrix(x_marker %*% membership)
    marker_detected[, active_profiles] <- marker_detected[, active_profiles, drop = FALSE] +
      as.matrix((x_marker != 0) %*% membership)
  }

  selected_processed <- selected_processed + length(selected_local)
  if (row_end %% CHECKPOINT_EVERY_OBJECT_ROWS == 0L || row_end == n_object_cells) {
    checkpoint <- list(
      mode = MODE,
      n_object_cells = n_object_cells,
      n_profiles = n_profiles,
      n_universe = n_universe,
      next_row_start = row_end + 1L,
      selected_processed = selected_processed,
      validation_done = validation_done,
      started = started,
      rank_sum_nonzero = rank_sum_nonzero,
      rank_count_nonzero = rank_count_nonzero,
      profile_n_cells = profile_n_cells,
      gene_sum = gene_sum,
      gene_sum_sq = gene_sum_sq,
      gene_detected = gene_detected,
      marker_sum = marker_sum,
      marker_detected = marker_detected
    )
    checkpoint_tmp <- paste0(CHECKPOINT_FILE, ".tmp")
    saveRDS(checkpoint, checkpoint_tmp, compress = FALSE)
    if (!file.rename(checkpoint_tmp, CHECKPOINT_FILE)) stop("Could not commit checkpoint")
    rm(checkpoint)
  }
  if ((row_start - 1L) %% (OUTER_CHUNK * 10L) == 0L || row_end == n_object_cells) {
    message(sprintf(
      "%s: processed object rows %d/%d; selected cells %d/%d",
      MODE, row_end, n_object_cells, selected_processed, nrow(meta)
    ))
    gc()
  }
  rm(
    matrix_chunk, rankings, scores, score_dt, membership, rank_sub,
    x_universe, indices, values, cell_j
  )
}

if (selected_processed != nrow(meta)) {
  stop(sprintf("Processed %d selected cells but expected %d", selected_processed, nrow(meta)))
}
if (any(profile_n_cells <= 0)) stop("At least one observed profile was not processed")

mean_truncated_rank <- (
  rank_sum_nonzero +
    (matrix(profile_n_cells, nrow = n_universe, ncol = n_profiles, byrow = TRUE) -
      rank_count_nonzero) * RMAX
) / matrix(profile_n_cells, nrow = n_universe, ncol = n_profiles, byrow = TRUE)
rownames(mean_truncated_rank) <- universe_genes
colnames(mean_truncated_rank) <- profile_dt$profile_index

aggregate <- list(
  dataset = DATASET,
  mode = MODE,
  rmax = RMAX,
  axis_genes = axis_genes,
  common_genes = common_genes,
  universe_genes = universe_genes,
  profiles = profile_dt,
  profile_n_cells = profile_n_cells,
  mean_truncated_rank = mean_truncated_rank,
  implementation = "UCell::StoreRankings_UCell 2.14.0; exact u-statistic aggregation"
)
saveRDS(aggregate, RANK_AGG_FILE, compress = "xz")

n_selected <- nrow(meta)
gene_properties <- data.table(
  gene = universe_genes,
  mean_expression = gene_sum / n_selected,
  expression_variance = pmax(0, gene_sum_sq / n_selected - (gene_sum / n_selected)^2),
  detection_rate = gene_detected / n_selected,
  n_cells = n_selected,
  matrix = MATRIX_GROUP
)
fwrite(gene_properties, GENE_PROPERTY_FILE, sep = "\t")

if (MODE == "subtype" && n_markers > 0L) {
  marker_rows <- vector("list", n_profiles)
  for (j in seq_len(n_profiles)) {
    marker_rows[[j]] <- data.table(
      dataset = DATASET,
      region = "DFC",
      donor = profile_dt$donor[[j]],
      Supertype = profile_dt$annotation[[j]],
      gene = marker_genes,
      mean_expression = marker_sum[, j] / profile_n_cells[[j]],
      detection_rate = marker_detected[, j] / profile_n_cells[[j]],
      n_nuclei = profile_n_cells[[j]]
    )
  }
  marker_out <- rbindlist(marker_rows)
  marker_out <- merge(marker_out, marker_dt[, .(marker_family, gene)], by = "gene", all.x = TRUE)
  setcolorder(
    marker_out,
    c(
      "dataset", "region", "donor", "Supertype", "marker_family", "gene",
      "mean_expression", "detection_rate", "n_nuclei"
    )
  )
  fwrite(marker_out, MARKER_PROPERTY_FILE, sep = "\t")
} else {
  fwrite(
    data.table(note = "not_applicable_to_broad_Subclass_run"),
    MARKER_PROPERTY_FILE,
    sep = "\t"
  )
}

capture.output(sessionInfo(), file = SESSION_FILE)
writeLines(
  c(
    paste0("dataset=", DATASET),
    paste0("matrix=", MATRIX_GROUP),
    paste0("selected_cells=", selected_processed),
    paste0("profiles=", n_profiles),
    paste0("rmax=", RMAX),
    paste0("started=", format(started, "%Y-%m-%dT%H:%M:%S%z")),
    paste0("completed=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
  ),
  COMPLETION_FILE
)
if (file.exists(CHECKPOINT_FILE) && !file.remove(CHECKPOINT_FILE)) {
  warning("Run completed, but checkpoint file could not be removed: ", CHECKPOINT_FILE)
}
