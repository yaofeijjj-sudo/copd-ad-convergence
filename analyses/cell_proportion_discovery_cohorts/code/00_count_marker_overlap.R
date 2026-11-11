.libPaths(c("_workspace/r_libs_deconv_20261002","_workspace/r_libs",.libPaths()))
suppressMessages({library(BRETIGEA); library(org.Hs.eg.db)})
cat("BRETIGEA", as.character(packageVersion("BRETIGEA")), "; org.Hs.eg.db", as.character(packageVersion("org.Hs.eg.db")), "\n")
m <- markers_df_human_brain
cat("columns:", colnames(m), "; rows:", nrow(m), "\n"); print(table(m$cell))
ax <- read.csv("03_candidate_analysis/WGCNA_original_modules_Entrez_correction_20261001/results/cross/AxisA_confirmed_genes.csv", colClasses="character")
# symbol -> Entrez: unique official symbol only, no alias
sym2eg <- AnnotationDbi::select(org.Hs.eg.db, keys=unique(m$markers), keytype="SYMBOL", columns="ENTREZID")
n_per <- table(sym2eg$SYMBOL); amb <- names(n_per)[n_per>1]
sym2eg <- sym2eg[!sym2eg$SYMBOL %in% amb, ]
m$entrez <- sym2eg$ENTREZID[match(m$markers, sym2eg$SYMBOL)]
cat("markers unmapped (not official symbol or ambiguous):", sum(is.na(m$entrez)), "of", nrow(m), "\n")
for (cl in unique(m$cell)) {
  mm <- m[m$cell==cl,]  # ranked order as in package
  for (k in c(50, 100, 1000)) {
    top <- head(mm, k); ov <- top$markers[!is.na(top$entrez) & top$entrez %in% ax$entrez_id]
    cat(sprintf("%s top%-4d overlap=%d %s\n", cl, k, length(ov), if (k==50) paste(ov, collapse=",") else ""))
  }
}
cat("findCells default nMarker:", formals(findCells)$nMarker, "\n")
mic <- m[m$cell=="mic",]
cat("mic top50 unmapped:", paste(head(mic,50)$markers[is.na(head(mic,50)$entrez)], collapse=","), "\n")
cat("mic top50 raw-symbol matches to Axis A current symbols:", sum(head(mic,50)$markers %in% ax$current_symbol), "\n")
cat("mic top1000 unmapped count:", sum(is.na(mic$entrez)), "\n")
br <- read.csv("03_candidate_analysis/WGCNA_original_modules_Entrez_correction_20261001/results/identity/brain_confirmed_module_gene_ids.csv", colClasses="character")
blue <- br$entrez_id[br$module=="blue"]; cat("blue confirmed IDs:", length(blue), "\n")
for (k in c(50,100,1000)) { top <- head(mic,k); cat(sprintf("mic top%d overlap with blue module=%d\n", k, sum(!is.na(top$entrez) & top$entrez %in% blue))) }
nb <- mic[!is.na(mic$entrez) & !(mic$entrez %in% blue),]; cat("mic markers remaining after removing blue (of 1000):", nrow(nb), "; first-50 cutoff reaches original rank", match(nb$markers[50], mic$markers), "\n")
na <- mic[!(!is.na(mic$entrez) & mic$entrez %in% ax$entrez_id),]; cat("version B: 50th remaining marker at original rank", match(na$markers[50], mic$markers), "\n")
