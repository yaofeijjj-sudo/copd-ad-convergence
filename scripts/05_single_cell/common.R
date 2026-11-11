source("R/runtime.R")
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})
ROOT <- "."
OUT <- file.path(ROOT, "results/single_cell/results/lung")
LUNG <- file.path(ROOT, "results/single_cell/preprocessing/lung")
for (d in c("tables", "qa", "objects", "figures", "docs", "logs")) dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)
RMAX <- params$ucell_max_rank
AXIS <- readLines(file.path(ROOT, "results/single_cell/results/lung_axis_measured.txt"))
UNIVERSE <- readLines(file.path(ROOT, "results/single_cell/results/lung_universe_measured.txt"))
stopifnot(length(AXIS) == 136, length(UNIVERSE) == 3608)
write_tsv <- function(x, name) fwrite(x, file.path(OUT, "tables", name), sep = "\t", na = "NA")
read_counts <- function(path) {
  f <- hdf5r::H5File$new(path, "r")
  on.exit(f$close_all())
  raw <- new("dgCMatrix",
    i = as.integer(f[["matrix/indices"]][]), p = as.integer(f[["matrix/indptr"]][]),
    x = as.numeric(f[["matrix/data"]][]), Dim = as.integer(f[["matrix/shape"]][])
  )
  symbols <- as.character(f[["matrix/features/id"]][])
  barcodes <- as.character(f[["matrix/barcodes"]][])
  keep <- !is.na(symbols) & nzchar(symbols)
  raw <- raw[keep, , drop = FALSE]
  symbols <- symbols[keep]
  u <- unique(symbols)
  m <- sparseMatrix(i = match(symbols, u), j = seq_along(symbols), x = 1, dims = c(length(u), length(symbols)))
  counts <- as(m %*% raw, "dgCMatrix")
  rownames(counts) <- u
  colnames(counts) <- barcodes
  counts
}
fit_model <- function(profile, types) {
  library(clubSandwich)
  d <- copy(profile[n_cells >= 20 & annotation %chin% types])
  setorder(d, donor, annotation)
  d[, type := factor(annotation, levels = types)]
  d[, donor_factor := factor(donor)]
  fit <- lm(mean_ucell ~ 0 + type + donor_factor, data = d)
  stopifnot(fit$rank == length(coef(fit)), !anyNA(coef(fit)))
  V <- vcovCR(fit, cluster = d$donor, type = "CR2")
  list(data = d, fit = fit, V = V, types = types, X = model.matrix(fit))
}
contrast <- function(model, weights, cohort, label, role) {
  c <- rep(0, length(coef(model$fit)))
  c[seq_along(model$types)] <- weights[model$types]
  stopifnot(!anyNA(c), abs(sum(c)) < 1e-10)
  estimate <- sum(c * coef(model$fit))
  se <- sqrt(drop(t(c) %*% model$V %*% c))
  wt <- clubSandwich::Wald_test(model$fit, constraints = matrix(c, nrow = 1), vcov = model$V, test = "HTZ")
  df <- wt$df_denom
  p <- wt$p_val
  old_se <- sqrt(drop(t(c) %*% vcov(model$fit) %*% c))
  list(
    row = data.table(
      cohort = cohort, contrast = label, role = role, effect = estimate, se_CR2 = se, df_Satterthwaite = df,
      ci_low = estimate - qt(.975, df) * se, ci_high = estimate + qt(.975, df) * se, p = p,
      n_donors = uniqueN(model$data$donor), n_profiles = nrow(model$data),
      ordinary_se = old_se, ordinary_residual_df = df.residual(model$fit)
    ),
    cvec = c, profile_weights = drop(model$X %*% solve(crossprod(model$X), c))
  )
}
target_contrasts <- function(model, target, cohort, prefix) {
  others <- setdiff(model$types, target)
  w <- setNames(rep(-1 / length(others), length(model$types)), model$types)
  w[target] <- 1
  main <- contrast(model, w, cohort, paste0(target, " minus mean(other types)"), "primary")
  pairs <- rbindlist(lapply(others, function(other) {
    w <- setNames(rep(0, length(model$types)), model$types)
    w[target] <- 1
    w[other] <- -1
    contrast(model, w, cohort, paste0(target, " minus ", other), "pairwise")$row
  }))
  pairs[, BH_FDR := p.adjust(p, "BH")]
  write_tsv(main$row, paste0(prefix, "_primary_contrast.tsv"))
  write_tsv(pairs, paste0(prefix, "_pairwise_contrasts.tsv"))
  pw <- model$data[, .(donor, annotation, n_cells, mean_ucell)]
  pw[, contrast_weight := main$profile_weights]
  write_tsv(pw, paste0(prefix, "_model_profile_weights.tsv"))
  list(main = main, pairs = pairs)
}
