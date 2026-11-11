source("R/runtime.R")
.libPaths(c("env/library", .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(GO.db)
  library(dplyr)
  library(stringr)
  library(readr)
  library(tidyr)
})
stopifnot(as.character(packageVersion("org.Hs.eg.db")) == "3.22.0", as.character(packageVersion("GO.db")) == "3.22.0")
out <- "results/enrichment"
cross <- "results/module_overlap/results/cross"
ids <- function(n) as.character(fread(file.path(cross, n))$entrez_id)
axes <- list(AxisA = ids("AxisA_confirmed_genes.csv"), AxisB = ids("AxisB_confirmed_genes.csv"))
bgs <- list(primary = ids("shared_non_grey_background.csv"), sensitivity = ids("shared_all_module_gene_ids.csv"))
stopifnot(lengths(axes) == c(137, 434), lengths(bgs) == c(3623, 4832), all(vapply(axes, function(x) all(x %in% bgs$primary), logical(1))))
clean_symbols <- function(x) sort(unique(na.omit(as.character(x))))
source("R/overrepresentation.R")
valid <- intersect(bgs$sensitivity, AnnotationDbi::keys(org.Hs.eg.db, keytype = "ENTREZID"))
g <- as.data.table(AnnotationDbi::select(org.Hs.eg.db, keys = valid, keytype = "ENTREZID", columns = c("GOALL", "ONTOLOGYALL")))[ONTOLOGYALL == "BP" & !is.na(GOALL)]
gn <- as.data.table(AnnotationDbi::select(GO.db, keys = unique(g$GOALL), keytype = "GOID", columns = "TERM"))
g <- merge(g, gn, by.x = "GOALL", by.y = "GOID")
go <- unique(g[, .(term_id = GOALL, term_name = TERM, collection = "GO_Biological_Process", db_source = "org.Hs.eg.db_3.22.0", gene = as.character(ENTREZID))])
cache <- "data/reference/enrichment/full_convergence_enrichment/annotation_cache"
convert_cache <- function(file) {
  d <- fread(file.path(cache, file))
  keys <- intersect(unique(d$gene), AnnotationDbi::keys(org.Hs.eg.db, keytype = "SYMBOL"))
  mp <- as.data.table(AnnotationDbi::select(org.Hs.eg.db, keys = keys, keytype = "SYMBOL", columns = "ENTREZID"))
  mp <- unique(mp[!is.na(ENTREZID)])
  mp[, n_ids := uniqueN(ENTREZID), by = SYMBOL]
  fwrite(mp, file.path(out, "provenance", paste0(file, ".identity_ledger.csv")))
  d <- merge(d, mp[n_ids == 1, .(gene = SYMBOL, entrez = as.character(ENTREZID))], by = "gene")
  d[, gene := entrez]
  d[, entrez := NULL]
  unique(d)
}
resources <- list(GO_Biological_Process = list(data = go, min = 5L, max = 1500L), KEGGREST_hsa = list(data = convert_cache("kegg_hsa_term2gene_cached.csv"), min = 5L, max = 800L), MSigDB_Hallmark = list(data = fread(file.path(out, "provenance/MSigDB_Hallmark_2026.1.Hs_entrez.csv"), colClasses = c(gene = "character")), min = 5L, max = 1000L))
rows <- list()
for (bg in names(bgs)) {
  for (ax in names(axes)) {
    for (db in names(resources)) {
      s <- resources[[db]]
      d <- as.data.table(complete_family_ora(s$data, axes[[ax]], bgs[[bg]], db, s$min, s$max))
      d[, `:=`(background_role = bg, gene_set = ax, family_id = paste(bg, ax, db, sep = "::"))]
      rows[[length(rows) + 1]] <- d
    }
  }
}
a <- rbindlist(rows)
stopifnot(uniqueN(a$family_id) == 12, all(a[zero_query_hit == TRUE]$p_value == 1), all(a[zero_query_hit == TRUE]$fdr_bh_complete_family == 1))
fwrite(a, file.path(out, "results/complete_family_all_terms.csv"))
fwrite(a[, .(eligible_terms = .N, zero_hit_terms = sum(zero_query_hit), significant_terms = sum(fdr_bh_complete_family < .05), min_fdr = min(fdr_bh_complete_family)), by = .(gene_set, source, background_role, family_id)], file.path(out, "results/family_summary.csv"))
for (db in names(resources)) fwrite(resources[[db]]$data, file.path(out, "provenance", paste0(db, "_term2entrez.csv")))
saveRDS(resources, file.path(out, "results/annotation_resources.rds"))
writeLines(capture.output(sessionInfo()), file.path(out, "provenance/sessionInfo.txt"))
print(a[, .(terms = .N, significant = sum(fdr_bh_complete_family < .05)), by = .(gene_set, source, background_role)])
