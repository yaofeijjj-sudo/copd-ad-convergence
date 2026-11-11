source("R/runtime.R")
.libPaths(c("./env/library", .libPaths()))
suppressPackageStartupMessages({
  library(UCell)
  library(jsonlite)
})
stopifnot(as.character(packageVersion("UCell")) == "2.14.0")
task_dir <- commandArgs(trailingOnly = TRUE)[1]
genes <- readLines(file.path(task_dir, "genes.txt"))
cells <- readLines(file.path(task_dir, "cells.txt"))
axis <- readLines(file.path(task_dir, "axis.txt"))
read_matrix <- function(path) {
  con <- if (grepl("\\.pipe$", path)) fifo(path, "rb", blocking = TRUE) else file(path, "rb")
  on.exit(close(con))
  total <- length(genes) * length(cells)
  x <- numeric(total)
  pos <- 1L
  while (pos <= total) {
    chunk <- readBin(con, what = "double", n = min(1048576, total - pos + 1), size = 8, endian = "little")
    if (!length(chunk)) stop("Expression transport ended before complete matrix")
    x[pos:(pos + length(chunk) - 1)] <- chunk
    pos <- pos + length(chunk)
  }
  stopifnot(length(x) == length(genes) * length(cells))
  matrix(x, nrow = length(genes), dimnames = list(genes, cells))
}
baseline <- read_matrix(file.path(task_dir, "baseline.bin"))
score <- function(m, tf) {
  keep <- genes != tf
  feature <- setdiff(axis, tf)
  stopifnot(sum(keep) >= 1500, length(intersect(feature, genes[keep])) > 0)
  ix <- split(seq_along(cells), ceiling(seq_along(cells) / 8000))
  unlist(lapply(ix, function(ii) {
    s <- ScoreSignatures_UCell(
      matrix = m[keep, ii, drop = FALSE],
      features = list(AxisA = feature), maxRank = params$ucell_max_rank, missing_genes = "skip",
      ties.method = "average", chunk.size = params$celloracle_ucell_chunk, ncores = params$celloracle_ucell_cores
    )
    stopifnot(identical(rownames(s), cells[ii]), all(is.finite(s)), all(s >= 0 & s <= 1))
    as.numeric(s[, 1])
  }), use.names = FALSE)
}
cin <- file("stdin", "r")
cat("READY\n")
flush.console()
repeat {
  line <- readLines(cin, n = 1, warn = FALSE)
  if (!length(line) || line == "QUIT") break
  parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
  tf <- parts[1]
  mode <- parts[2]
  t0 <- proc.time()[["elapsed"]]
  baseline_seconds <- 0
  if (mode == "main") {
    tf <- tf
    bs <- proc.time()[["elapsed"]]
    base_score <- score(baseline, tf)
    baseline_seconds <- proc.time()[["elapsed"]] - bs
    con <- file(file.path(task_dir, "baseline_score.bin"), "wb")
    writeBin(base_score, con, size = 8, endian = "little")
    close(con)
  } else {
    stopifnot(identical(tf, tf))
  }
  transport_start <- proc.time()[["elapsed"]]
  predicted <- read_matrix(file.path(task_dir, "simulated.pipe"))
  transport_seconds <- proc.time()[["elapsed"]] - transport_start
  ss <- proc.time()[["elapsed"]]
  simulated_score <- score(predicted, tf)
  simulated_seconds <- proc.time()[["elapsed"]] - ss
  rm(predicted)
  gc(verbose = FALSE)
  con <- file(file.path(task_dir, "simulated_score.bin"), "wb")
  writeBin(simulated_score, con, size = 8, endian = "little")
  close(con)
  write_json(list(baseline_score_seconds = baseline_seconds, matrix_read_seconds = transport_seconds, simulated_score_seconds = simulated_seconds, total_R_seconds = proc.time()[["elapsed"]] - t0), file.path(task_dir, paste0(tf, "_", mode, "_timing.json")), auto_unbox = TRUE)
  cat("DONE\n")
  flush.console()
}
close(cin)
