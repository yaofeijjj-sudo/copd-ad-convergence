# Platform tables through GEOquery, KEGG pathway membership through KEGGREST and MSigDB Hallmark through msigdbr.
# KEGG is not redistributed; the study used the release current on 2026-06-25, so a later download can differ.
# Run from the repository root after 01_fetch.py.
suppressPackageStartupMessages({
  library(data.table)
  library(GEOquery)
  library(KEGGREST)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(msigdbr)
})
paths <- yaml::read_yaml("config/paths.yml")
root <- if (startsWith(paths$data_root, "/")) paths$data_root else file.path(getwd(), paths$data_root)
tmp <- file.path(root, "reference/geo_soft")
dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
out <- function(p) {
  f <- file.path(root, p)
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  f
}

# GPL96 and GPL570 use the GEO annotation tables; GPL14550 uses the submitter table, which carries the GENE column.
for (g in c("GPL96", "GPL570")) {
  fwrite(as.data.table(Table(getGEO(g, AnnotGPL = TRUE, destdir = tmp))), out(sprintf("raw/platform_annotations/%s_annot_table.csv", g)))
}
fwrite(as.data.table(Table(getGEO("GPL14550", destdir = tmp))), out("raw/platform_annotations/GPL14550_annot_table.csv"))
for (g in c("GPL10558", "GPL4372")) {
  t <- as.data.table(Table(getGEO(g, destdir = tmp)))
  for (k in names(t)) set(t, which(t[[k]] == ""), k, NA)
  fwrite(t, out(sprintf("reference/gene_identity/references/%s_original.tsv.gz", g)), sep = "\t", na = "")
}

cache <- "reference/enrichment/full_convergence_enrichment/annotation_cache"
links <- keggLink("pathway", "hsa")
names_kegg <- keggList("pathway", "hsa")
link <- unique(data.table(entrez = sub("^hsa:", "", names(links)), term_id = sub("^path:", "", as.character(links))))
term <- data.table(
  term_id = sub("^path:", "", names(names_kegg)),
  term_name = sub(" - Homo sapiens \\(human\\)$", "", as.character(names_kegg))
)
sym <- as.data.table(AnnotationDbi::select(org.Hs.eg.db, keys = unique(link$entrez), keytype = "ENTREZID", columns = "SYMBOL"))
sym <- unique(sym[!is.na(SYMBOL), .(entrez = ENTREZID, gene = toupper(SYMBOL))])
kegg <- merge(merge(link, sym, by = "entrez"), term, by = "term_id")
kegg <- unique(kegg[, .(term_id = toupper(term_id), term_name, collection = "KEGGREST_hsa", db_source = "KEGG", gene)])
fwrite(kegg, out(file.path(cache, "kegg_hsa_term2gene_cached.csv")))
cat("KEGG rows", nrow(kegg), "terms", uniqueN(kegg$term_id), "(2026-06-25 release: 39523 rows, 372 terms)\n")

h <- as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))
stopifnot(all(h$db_version == "2026.1.Hs"))
name <- tolower(gsub("_", " ", sub("^HALLMARK_", "", h$gs_name)))
h <- unique(data.table(
  term_id = h$gs_name, term_name = paste0(toupper(substr(name, 1, 1)), substring(name, 2)),
  collection = "MSigDB_Hallmark_2026.1.Hs", db_source = "HALLMARK", gene = h$gene_symbol
))
fwrite(h, out(file.path(cache, "msigdb_hallmark_human_term2gene_cached.csv")))
