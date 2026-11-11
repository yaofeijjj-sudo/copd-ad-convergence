source("R/runtime.R")

suppressPackageStartupMessages({
  library(data.table)
  library(hdf5r)
  library(digest)
})

project_root <- "."
result_dir <- file.path(
  project_root,
  "results/single_cell/preprocessing",
  "lung"
)
soft_file <- file.path(project_root, "data/raw/single_cell/copd_lung_gse302339/GSE302339_family.soft.gz")
h5_dir <- file.path(project_root, "data/raw/lung_h5")
axis_file <- file.path(project_root, "resources/AxisA_confirmed_genes.csv")

dir.create(file.path(result_dir, "inputs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(result_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(result_dir, "qa"), recursive = TRUE, showWarnings = FALSE)

read_soft_samples <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  starts <- grep("^\\^SAMPLE = ", lines)
  ends <- c(starts[-1L] - 1L, length(lines))
  rows <- vector("list", length(starts))
  for (i in seq_along(starts)) {
    block <- lines[starts[i]:ends[i]]
    gsm <- sub("^\\^SAMPLE = ", "", block[1L])
    get_one <- function(prefix) {
      hit <- grep(paste0("^", prefix), block, value = TRUE)
      if (!length(hit)) {
        return(NA_character_)
      }
      sub(paste0("^", prefix), "", hit[1L])
    }
    chars <- grep("^!Sample_characteristics_ch1 = ", block, value = TRUE)
    chars <- sub("^!Sample_characteristics_ch1 = ", "", chars)
    char_names <- trimws(sub(":.*$", "", chars))
    char_values <- trimws(sub("^[^:]+:\\s*", "", chars))
    char_map <- setNames(char_values, tolower(char_names))
    get_char <- function(name) {
      out <- unname(char_map[[tolower(name)]])
      if (is.null(out)) NA_character_ else out
    }
    supp <- grep("^!Sample_supplementary_file_[0-9]+ = ", block, value = TRUE)
    supp <- if (length(supp)) sub("^!Sample_supplementary_file_[0-9]+ = ", "", supp[1L]) else NA_character_
    rows[[i]] <- data.table(
      gsm = gsm,
      title = get_one("!Sample_title = "),
      subject_group = get_char("subject group"),
      subject_id_raw = get_char("subject identity"),
      lab_id = get_char("lab id"),
      steroidstype = get_char("steroidstype"),
      smoke = get_char("smoke"),
      gold = get_char("gold"),
      age = suppressWarnings(as.numeric(get_char("age"))),
      sex = toupper(get_char("sex")),
      copd = get_char("copd"),
      fev1 = suppressWarnings(as.numeric(get_char("fev1"))),
      fev1_fvc = suppressWarnings(as.numeric(get_char("fev1/fvc"))),
      quitdate = get_char("quitdate"),
      steroids = get_char("steroids"),
      lobe = toupper(get_char("lobe")),
      lobe_laa_percent = suppressWarnings(as.numeric(get_char("lobe emphysema_severity_laa%"))),
      lobe_emphysema = get_char("lobe emphysema"),
      packyear = suppressWarnings(as.numeric(get_char("packyear"))),
      supplementary_url = supp
    )
  }
  rbindlist(rows, fill = TRUE)
}

normalize_subject_id <- function(x) {
  x <- tolower(trimws(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

inspect_h5 <- function(path) {
  f <- H5File$new(path, mode = "r")
  on.exit(f$close_all(), add = TRUE)
  shape <- as.integer(f[["matrix/shape"]][])
  genes <- as.character(f[["matrix/features/name"]][])
  feature_types <- as.character(f[["matrix/features/feature_type"]][])
  barcodes <- as.character(f[["matrix/barcodes"]][])
  indptr <- as.numeric(f[["matrix/indptr"]][])
  data_head <- as.numeric(f[["matrix/data"]][seq_len(min(10000L, f[["matrix/data"]]$dims))])
  data.table(
    matrix_file = basename(path),
    gsm = sub("_.*$", "", basename(path)),
    n_features = shape[1L],
    n_cells = shape[2L],
    n_nonzero = tail(indptr, 1L),
    n_duplicate_gene_symbols = sum(duplicated(genes)),
    n_duplicate_barcodes = sum(duplicated(barcodes)),
    feature_type_values = paste(sort(unique(feature_types)), collapse = ";"),
    sampled_counts_nonnegative = all(data_head >= 0),
    sampled_counts_integer = all(data_head == floor(data_head))
  )
}

ledger <- read_soft_samples(soft_file)
ledger[, subject_id := normalize_subject_id(subject_id_raw)]
ledger[, disease_group := fifelse(tolower(copd) == "copd", "COPD", "Control")]
ledger[, emphysema_group := fifelse(
  grepl("nonemphysematous", tolower(subject_group)), "NE-COPD",
  fifelse(grepl("emphysematous", tolower(subject_group)), "E-COPD", "Control")
)]

h5_files <- sort(list.files(h5_dir, pattern = "\\.h5$", full.names = TRUE))
inventory <- rbindlist(lapply(h5_files, inspect_h5))
ledger <- merge(ledger, inventory, by = "gsm", all.x = TRUE, sort = FALSE)
ledger[, matrix_available := !is.na(matrix_file)]

conflict_fields <- c("disease_group", "emphysema_group", "age", "sex", "smoke", "gold", "fev1", "fev1_fvc", "packyear")
conflicts <- rbindlist(lapply(conflict_fields, function(field) {
  ledger[matrix_available == TRUE,
    {
      vals <- unique(na.omit(as.character(get(field))))
      .(field = field, n_distinct = length(vals), values = paste(sort(vals), collapse = ";"))
    },
    by = subject_id
  ][n_distinct > 1L]
}), fill = TRUE)

axis <- fread(axis_file)
axis_genes <- unique(na.omit(as.character(axis$current_symbol)))
axis_genes <- axis_genes[nzchar(axis_genes)]
fwrite(data.table(gene = axis_genes), file.path(result_dir, "inputs", "axisA_genes.csv"))

checks <- data.table(
  check = c(
    "soft_sample_records", "h5_matrices", "soft_samples_with_matrix", "soft_samples_without_matrix",
    "unique_donors_with_matrix", "all_h5_join_to_soft", "all_shapes_positive",
    "no_duplicate_barcodes_within_matrix", "sampled_counts_nonnegative_integer",
    "donor_level_core_phenotype_conflicts", "axis_gene_count", "axis_source_sha256"
  ),
  observed = c(
    nrow(ledger), nrow(inventory), sum(ledger$matrix_available), sum(!ledger$matrix_available),
    uniqueN(ledger[matrix_available == TRUE]$subject_id),
    all(inventory$gsm %in% ledger$gsm),
    all(inventory$n_features > 0 & inventory$n_cells > 0),
    all(inventory$n_duplicate_barcodes == 0),
    all(inventory$sampled_counts_nonnegative & inventory$sampled_counts_integer),
    nrow(conflicts[field %in% c("disease_group", "emphysema_group", "age", "sex")]),
    length(axis_genes), digest(file = axis_file, algo = "sha256")
  ),
  expected = c("67", "65", "65", "2", "reported", "TRUE", "TRUE", "TRUE", "TRUE", "0", "137", digest(file = axis_file, algo = "sha256"))
)

fwrite(ledger, file.path(result_dir, "tables", "sample_ledger.tsv"), sep = "\t", na = "NA")
fwrite(inventory, file.path(result_dir, "tables", "matrix_inventory.tsv"), sep = "\t")
fwrite(conflicts, file.path(result_dir, "qa", "donor_characteristics_conflicts.tsv"), sep = "\t", na = "NA")
fwrite(checks, file.path(result_dir, "qa", "input_checks.tsv"), sep = "\t")

cat("Sample ledger complete\n")
print(checks)
cat("\nSamples without count matrices:\n")
print(ledger[matrix_available == FALSE, .(gsm, title, subject_id_raw, supplementary_url)])
cat("\nRepeated donors:\n")
print(ledger[matrix_available == TRUE, .N, by = subject_id][N > 1L][order(-N, subject_id)])
cat("\nDonor characteristic conflicts:\n")
print(conflicts)
