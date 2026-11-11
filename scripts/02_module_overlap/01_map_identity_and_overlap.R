source("R/runtime.R")
PACK <- file.path(getwd(), "results/module_overlap")
ROOT <- getwd()
options(stringsAsFactors = FALSE, warn = 1, error = function() {
  traceback(3)
  quit(status = 1)
})
suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})
cfg <- fromJSON(file.path("config", "identity.json"))
old <- file.path(ROOT, "results/wgcna")
snapshot <- file.path(ROOT, "data/reference/gene_identity")
result <- function(name) file.path(PACK, "results", name)
stamp <- function(...) message(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), " ", ...)
check <- function(name, ok) if (!isTRUE(ok)) stop("Check failed: ", name)
use <- function(p, role) normalizePath(p, mustWork = TRUE)
read <- function(p, role) fread(use(p, role), colClasses = "character", na.strings = c("", "NA"))
stamp("Entrez identity of network features")

ref <- function(n) file.path(snapshot, "references", n)
info <- fread(ref("gene_info.gz"), colClasses = "character", na.strings = NULL, quote = "")
history_table <- fread(ref("gene_history.tsv.gz"), colClasses = "character", na.strings = NULL, quote = "")
fwrite(info[, .(records = .N), by = "#tax_id"], file.path(PACK, "provenance/gene_info_species_counts.csv"))
info <- info[info[["#tax_id"]] == "9606"]
history_table <- history_table[history_table[[1]] == "9606"]
check("human_gene_info_only", all(info[["#tax_id"]] == "9606"))
check("unique_GeneID", !anyDuplicated(info$GeneID))
current <- list2env(as.list(setNames(info$Symbol, info$GeneID)), parent = emptyenv(), hash = TRUE)
history <- list2env(as.list(setNames(history_table$GeneID, history_table$Discontinued_GeneID)), parent = emptyenv(), hash = TRUE)
cache <- new.env(parent = emptyenv(), hash = TRUE)
resolve <- function(id) {
  if (exists(id, cache, inherits = FALSE)) {
    return(get(id, cache, inherits = FALSE))
  }
  raw <- id
  chain <- id
  seen <- character()
  while (!exists(id, current, inherits = FALSE) && exists(id, history, inherits = FALSE) && !id %in% seen) {
    seen <- c(seen, id)
    id <- get(id, history, inherits = FALSE)
    chain <- c(chain, id)
  }
  ans <- list(id = if (exists(id, current, inherits = FALSE)) id else NA_character_, history = paste(chain, collapse = ">"))
  assign(raw, ans, cache)
  ans
}
field <- function(s) {
  if (is.na(s) || trimws(s) %in% c("", "NA", "---", "-")) {
    return(list(ids = character(), bad = FALSE, history = ""))
  }
  tt <- unique(strsplit(trimws(s), "[[:space:]/;,|]+")[[1]])
  tt <- tt[nzchar(tt)]
  valid <- grepl("^[1-9][0-9]*$", tt)
  rr <- lapply(tt[valid], resolve)
  ii <- vapply(rr, `[[`, character(1), "id")
  list(ids = sort(unique(ii[!is.na(ii)])), bad = any(!valid) || anyNA(ii), history = paste(vapply(rr, `[[`, character(1), "history"), collapse = "|"))
}
annotate <- function(organ, probes) {
  platform <- if (organ == "lung") "GPL10558" else "GPL4372"
  original <- fread(ref(paste0(platform, "_original.tsv.gz")), colClasses = "character", quote = "", na.strings = c("", "NA", "---"))
  cached <- fread(ref(paste0(platform, ".annot.gz")), skip = "ID\t", colClasses = "character", quote = "", na.strings = c("", "NA", "---"))
  check(paste0(organ, "_unique_platform_probe_ids"), !anyDuplicated(original$ID) && !anyDuplicated(cached$ID))
  original_column <- if (organ == "lung") "Entrez_Gene_ID" else "EntrezGeneID"
  ann <- data.table(
    selected_probe_id = probes, original_entrez = original[[original_column]][match(probes, original$ID)],
    cached_entrez = cached[["Gene ID"]][match(probes, cached$ID)]
  )
  rr <- lapply(seq_len(nrow(ann)), function(i) {
    a <- field(ann$original_entrez[i])
    b <- field(ann$cached_entrez[i])
    id <- NA_character_
    status <- "missing_stable_ID"
    if (a$bad || b$bad) {
      status <- "unresolved_nonempty_ID"
    } else if (length(a$ids) > 1 || length(b$ids) > 1) {
      status <- "multiple_gene_targets"
    } else if (length(a$ids) == 1 && length(b$ids) == 1 && !identical(a$ids, b$ids)) {
      status <- "annotation_version_conflict"
    } else if (length(a$ids) == 1 && length(b$ids) == 1) {
      id <- a$ids
      status <- "original_cached_agree_after_history"
    } else if (length(a$ids) == 1) {
      id <- a$ids
      status <- "original_unique_cached_missing"
    } else if (length(b$ids) == 1) {
      id <- b$ids
      status <- "cached_unique_original_missing"
    }
    data.table(
      entrez_id = id, current_symbol = if (!is.na(id)) get(id, current, inherits = FALSE) else NA_character_,
      identity_status = status, original_history = a$history, cached_history = b$history, identity_confirmed = !is.na(id)
    )
  })
  cbind(ann, rbindlist(rr))
}
ledgers <- list()
memberships <- list()
summary_rows <- list()
trait_maps <- list()
for (organ in c("lung", "brain")) {
  p <- file.path(old, if (organ == "lung") "lung/COPD_WGCNA" else "brain/AD_WGCNA_unadjusted")
  m <- read(file.path(p, "module_assignment.csv"), paste0(organ, "_original_modules"))
  mp <- read(file.path(p, "training_probe_to_gene_map.csv"), paste0(organ, "_original_selected_probes"))
  use(file.path(p, "validation_probe_to_gene_map.csv"), paste0(organ, "_original_validation_selected_probes"))
  use(file.path(p, "gene_filter_table.csv"), paste0(organ, "_original_gene_filter"))
  use(file.path(p, "sample_table_used.csv"), paste0(organ, "_original_samples"))
  envfile <- use(file.path(p, if (organ == "lung") "copd_wgcna_env.RData" else "ad_wgcna_env.RData"), paste0(organ, "_original_saved_network"))
  check(paste0(organ, "_unique_original_features"), !anyDuplicated(m$gene_id))
  check(paste0(organ, "_unique_probe_map_keys"), !anyDuplicated(mp$gene_id))
  check(paste0(organ, "_feature_count"), nrow(m) == cfg[[paste0(organ, "_network_features")]])
  check(paste0(organ, "_module_count"), uniqueN(m[module != "grey", module]) == cfg[[paste0(organ, "_non_grey_modules")]])
  x <- data.table(organ = organ, original_gene_id = m$gene_id, original_symbol = m$gene_symbol, module = m$module)
  x[, selected_probe_id := mp$selected_probe_id[match(original_gene_id, mp$gene_id)]]
  check(paste0(organ, "_all_original_probes_found"), !anyNA(x$selected_probe_id))
  ann <- annotate(organ, x$selected_probe_id)
  x <- cbind(x, ann[, setdiff(names(ann), "selected_probe_id"), with = FALSE])
  x[, modules_per_ID := NA_integer_]
  x[identity_confirmed == TRUE, modules_per_ID := uniqueN(module), by = entrez_id]
  x[, overlap_eligible := identity_confirmed & !is.na(modules_per_ID) & modules_per_ID == 1L]
  x[, overlap_status := ifelse(!identity_confirmed, identity_status, ifelse(modules_per_ID > 1, "ID_in_multiple_modules", "confirmed_unique_module_membership"))]
  mem <- unique(x[overlap_eligible == TRUE, .(entrez_id, current_symbol, module)])
  check(paste0(organ, "_one_confirmed_module_per_ID"), !anyDuplicated(mem$entrez_id))
  fwrite(x, file.path(result("identity"), paste0(organ, "_all_original_features_identity.csv")))
  fwrite(x[overlap_eligible == FALSE], file.path(result("identity"), paste0(organ, "_excluded_from_overlap.csv")))
  fwrite(x[identity_confirmed == TRUE & (duplicated(entrez_id) | duplicated(entrez_id, fromLast = TRUE))], file.path(result("identity"), paste0(organ, "_duplicate_ID_features.csv")))
  fwrite(mem, file.path(result("identity"), paste0(organ, "_confirmed_module_gene_ids.csv")))
  e <- new.env(parent = baseenv())
  load(envfile, envir = e)
  check(paste0(organ, "_saved_network_module_order"), identical(colnames(e$dat_expr), x$original_gene_id))
  check(paste0(organ, "_saved_colors_match"), identical(as.character(e$module_colors), x$module))
  check(paste0(organ, "_original_sample_count"), nrow(e$dat_expr) == cfg[[paste0(organ, "_samples")]])
  check(paste0(organ, "_original_validation_feature_order"), identical(colnames(e$valid_expr), colnames(e$dat_expr)))
  rm(e)
  gc()
  traitfile <- file.path(p, if (organ == "lung") "copd_module_trait_tidy.csv" else "ad_module_trait_tidy.csv")
  tr <- read(traitfile, paste0(organ, "_original_module_disease_association"))
  tm <- tr[trait == if (organ == "lung") "copd_status" else "Disease_State", .(module, r = as.numeric(correlation))]
  check(paste0(organ, "_unique_disease_trait_rows"), !anyDuplicated(tm$module))
  trait_maps[[organ]] <- setNames(tm$r, tm$module)
  ledgers[[organ]] <- x
  memberships[[organ]] <- mem
  summary_rows[[organ]] <- data.table(
    organ = organ, original_features = nrow(x), identity_confirmed_features = sum(x$identity_confirmed),
    unresolved_annotation_features = sum(!x$identity_confirmed), multiple_module_features = sum(x$identity_confirmed & x$modules_per_ID > 1, na.rm = TRUE),
    confirmed_overlap_IDs = nrow(mem), same_module_duplicate_copies = sum(x$overlap_eligible) - nrow(mem)
  )
  stamp(organ, " identity complete: ", nrow(x), " original features, ", nrow(mem), " eligible unique IDs")
}
all_ledger <- rbindlist(ledgers)
fwrite(all_ledger, file.path(result("identity"), "all_17680_original_features_identity.csv"))
fwrite(rbindlist(summary_rows), file.path(result("identity"), "identity_summary.csv"))
shared <- merge(memberships$lung, memberships$brain, by = "entrez_id", suffixes = c("_lung", "_brain"))
check("same_ID_same_display_label", all(shared$current_symbol_lung == shared$current_symbol_brain))
shared[, current_symbol := current_symbol_lung]
shared[, c("current_symbol_lung", "current_symbol_brain") := NULL]
check("unique_shared_ID", !anyDuplicated(shared$entrez_id))
bg <- shared[module_lung != "grey" & module_brain != "grey"]
fwrite(shared, file.path(result("cross"), "shared_all_module_gene_ids.csv"))
fwrite(bg, file.path(result("cross"), "shared_non_grey_background.csv"))
lm <- sort(unique(ledgers$lung[module != "grey", module]))
bm <- sort(unique(ledgers$brain[module != "grey", module]))
N <- nrow(bg)
cells <- list()
for (l in lm) {
  for (b in bm) {
    L <- bg[module_lung == l, entrez_id]
    B <- bg[module_brain == b, entrez_id]
    a <- length(intersect(L, B))
    bb <- length(L) - a
    cc <- length(B) - a
    dd <- N - a - bb - cc
    tab <- matrix(c(a, bb, cc, dd), nrow = 2, byrow = TRUE)
    ft <- stats::fisher.test(tab, alternative = "greater")
    lr <- unname(trait_maps$lung[l])
    br <- unname(trait_maps$brain[b])
    cells[[length(cells) + 1L]] <- data.table(
      lung_module = l, brain_module = b, overlap = a, lung_genes = length(L), brain_genes = length(B), background = N,
      odds_ratio = unname(ft$estimate), p_value = ft$p.value, fisher_a = a, fisher_b = bb, fisher_c = cc, fisher_d = dd,
      lung_disease_r = lr, brain_disease_r = br, direction = if (lr > 0 && br > 0) "both_positive" else if (lr < 0 && br < 0) "both_negative" else "opposite_or_zero"
    )
  }
}
pairs <- rbindlist(cells)
pairs[, fdr := stats::p.adjust(p_value, method = "BH")]
pairs[, BH_family_size := .N]
setorder(pairs, fdr, -overlap, lung_module, brain_module)
check("original_128_pair_testing_family", nrow(pairs) == 128L)
check("valid_2x2_counts", all(pairs[, fisher_a + fisher_b + fisher_c + fisher_d] == N) && all(pairs[, fisher_d] >= 0))
check("finite_Fisher_and_BH", all(is.finite(pairs$p_value)) && all(is.finite(pairs$fdr)))
check("Fisher_matches_direct_hypergeometric_tail", isTRUE(all.equal(pairs$p_value, stats::phyper(pairs$overlap - 1, pairs$lung_genes, N - pairs$lung_genes, pairs$brain_genes, lower.tail = FALSE), tolerance = 1e-10)))
fwrite(pairs, file.path(result("cross"), "all128_module_pairs_fisher_BH.csv"))
axis_rows <- list()
traces <- list()
for (axis in c("AxisA", "AxisB")) {
  mods <- cfg$axes[[axis]]
  l <- mods[1]
  b <- mods[2]
  gs <- bg[module_lung == l & module_brain == b]
  pair <- pairs[lung_module == l & brain_module == b]
  check(paste0(axis, "_gene_set_matches_Fisher_overlap"), nrow(gs) == pair$overlap)
  fwrite(gs, file.path(result("cross"), paste0(axis, "_confirmed_genes.csv")))
  axis_rows[[axis]] <- cbind(data.table(axis = axis, confirmed_count = nrow(gs), fixed_axis = TRUE), pair)
}
axes <- rbindlist(axis_rows)
fwrite(axes, file.path(result("cross"), "axis_summary.csv"))
fwrite(data.table(
  all_shared_confirmed_IDs = nrow(shared), non_grey_shared_confirmed_IDs = N, pairs = 128,
  significant_pairs = sum(pairs$fdr < .05), fixed_axis = TRUE
), file.path(result("cross"), "background_summary.csv"))
writeLines(capture.output(sessionInfo()), file.path(PACK, "logs/sessionInfo.txt"))
print(axes[, .(axis, confirmed_count, fdr)])
stamp("128 module pairs tested; background = ", N)
