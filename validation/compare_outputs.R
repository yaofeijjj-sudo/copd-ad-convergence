source("R/runtime.R")
suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})
spec <- fromJSON("validation/comparisons.json", simplifyVector = FALSE)
dir.create("results/validation", recursive = TRUE, showWarnings = FALSE)
rows <- list()
check <- function(ref, out, keys, columns = NULL) {
  a <- fread(ref)
  b <- fread(out)
  for (k in keys) {
    a[[k]] <- as.character(a[[k]])
    b[[k]] <- as.character(b[[k]])
  }
  stopifnot(!anyDuplicated(a[, ..keys]), !anyDuplicated(b[, ..keys]))
  if (basename(ref) == "lung_module_traits.csv") b <- b[a[, ..keys], on = keys]
  setorderv(a, keys)
  setorderv(b, keys)
  if (!identical(a[, ..keys], b[, ..keys])) {
    return(data.table(column = "keys", pass = FALSE, max_absolute_error = NA_real_))
  }
  if (is.null(columns)) columns <- setdiff(names(a), c(keys, "model", "db_source", "theme", "term_name", "display", "fontface", "text_colour"))
  rbindlist(lapply(columns, function(k) {
    if (!k %in% names(b)) {
      return(data.table(column = k, pass = FALSE, max_absolute_error = NA_real_))
    }
    x <- a[[k]]
    y <- b[[k]]
    if (is.numeric(x) && is.numeric(y)) {
      na_ok <- identical(is.na(x), is.na(y))
      finite <- is.finite(x) & is.finite(y)
      err <- abs(x[finite] - y[finite])
      # Zeroes, missing values and gene membership are exact; P values use relative tolerance.
      tol <- if (grepl("(^p$|p_value|p_two|fgsea_p|fdr|FDR|holm|Holm)", k)) 1e-10 * pmax(abs(x[finite]), .Machine$double.xmin) else 1e-12 + 1e-10 * abs(x[finite])
      # MAGMA gene-level P values of the strongest genes depend on the random seed, so MAGMA results get a 1% tolerance.
      if (basename(ref) == "magma.csv") tol <- 1e-12 + 1e-2 * abs(x[finite])
      ok <- na_ok && all(err <= tol) && all(x[!finite] == y[!finite], na.rm = TRUE)
      data.table(column = k, pass = ok, max_absolute_error = if (length(err)) max(err) else 0)
    } else {
      data.table(column = k, pass = identical(as.character(x), as.character(y)), max_absolute_error = NA_real_)
    }
  }))
}
for (s in spec) {
  ref <- file.path("validation/reference", s$reference)
  out <- file.path("results", s$output)
  if (!file.exists(out)) {
    rows[[length(rows) + 1L]] <- data.table(reference = s$reference, column = "output", pass = FALSE, status = "missing", max_absolute_error = NA_real_)
    next
  }
  d <- check(ref, out, unlist(s$keys), if (is.null(s$columns)) NULL else unlist(s$columns))
  d[, `:=`(reference = s$reference, status = ifelse(pass, "pass", "different"))]
  rows[[length(rows) + 1L]] <- d
}
for (axis in c("A", "B")) {
  ref <- fread(paste0("resources/Axis", axis, "_confirmed_genes.csv"), colClasses = "character")
  out <- fread(paste0("results/module_overlap/results/cross/Axis", axis, "_confirmed_genes.csv"), colClasses = "character")
  ok <- identical(sort(ref$entrez_id), sort(out$entrez_id))
  rows[[length(rows) + 1L]] <- data.table(reference = paste0("Axis", axis), column = "exact_Entrez_membership", pass = ok, status = ifelse(ok, "pass", "different"), max_absolute_error = 0)
}
d <- rbindlist(rows, fill = TRUE)
fwrite(d, "results/validation/comparison_details.tsv", sep = "\t")
fwrite(d[, .(columns = .N, passed = sum(pass), status = if (all(pass)) "pass" else paste(unique(status[!pass]), collapse = ";")), by = reference], "results/validation/comparison_summary.tsv", sep = "\t")
print(d[pass == FALSE])
cat(sum(d$pass), "/", nrow(d), "comparisons passed\n")
if (any(!d$pass)) quit(status = 1)
