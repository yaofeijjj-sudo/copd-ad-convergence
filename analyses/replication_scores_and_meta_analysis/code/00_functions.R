# 2026-10-04 rerun: scoring switched to the official singscore package (rankGenes + simpleScore); everything else unchanged.
# Shared functions, identical in logic to 02_run_specificity.R (which keeps its own copy because it was frozen and run first).
# fit_comparison() additionally returns the standard error and residual df, needed for meta-analysis.
suppressPackageStartupMessages({
  library(data.table); library(limma); library(edgeR); library(AnnotationDbi); library(org.Hs.eg.db)
})
ROOT <- normalizePath(".")
P <- function(x) file.path(ROOT, x)
axisA <- fread(P("03_candidate_analysis/WGCNA_original_modules_Entrez_correction_20261001/results/cross/AxisA_confirmed_genes.csv"), colClasses = "character")$entrez_id
stopifnot(length(axisA) == 137)
hist <- fread(cmd = paste("gzip -cd", shQuote(P("03_candidate_analysis/WGCNA_Entrez_training_only_20260930/references/gene_history.tsv.gz"))), colClasses = "character")
setnames(hist, 1, "tax_id"); hist <- hist[tax_id == "9606" & GeneID != "-"]
repl <- setNames(hist$GeneID, hist$Discontinued_GeneID)
current_id <- function(x) { for (i in 1:5) { hit <- !is.na(x) & x %in% names(repl); if (!any(hit)) break; x[hit] <- repl[x[hit]] }; x }
read_matrix <- function(file) {
  tab <- fread(cmd = paste("gzip -cd", shQuote(P(file)), "| sed -n '/^!series_matrix_table_begin/,$p' | sed '1d;/^!series_matrix_table_end/,$d'"), check.names = FALSE)
  m <- as.matrix(tab[, -1L, with = FALSE]); storage.mode(m) <- "numeric"; rownames(m) <- as.character(tab[[1L]]); m
}
keep_max_mad <- function(m, gene) {
  ok <- !is.na(gene); m <- m[ok, , drop = FALSE]; gene <- gene[ok]
  v <- apply(m, 1, mad, na.rm = TRUE); o <- order(gene, -v); o <- o[is.finite(v[o])]
  o <- o[!duplicated(gene[o])]; out <- m[o, , drop = FALSE]; rownames(out) <- gene[o]
  out[rowSums(is.finite(out)) == ncol(out), , drop = FALSE]
}
collapse_array <- function(m, annot, idcol = "Gene ID") {
  a <- if (grepl("\\.annot\\.gz$", annot)) {
    fread(cmd = paste("gzip -cd", shQuote(P(annot)), "| sed -n '/^ID\t/,$p' | grep -v '^!platform_table_end'"), colClasses = "character", select = c("ID", idcol))
  } else fread(P(annot), colClasses = "character", select = c("ID", idcol))
  setnames(a, c("probe", "gene")); a <- a[probe %in% rownames(m) & !is.na(gene) & nzchar(gene) & !grepl("/", gene)]
  keep_max_mad(m[a$probe, , drop = FALSE], current_id(a$gene))
}
unique_map <- function(keys, keytype) {
  mp <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = unique(keys), keytype = keytype, columns = "ENTREZID"))
  mp <- as.data.table(mp)[!is.na(ENTREZID)]; setnames(mp, 1, "key")
  mp <- mp[, if (.N == 1) .SD, by = key]
  current_id(mp$ENTREZID[match(keys, mp$key)])
}
series_chars <- function(file, ch = "ch1") {  # GEO sample characteristics -> data.table; key read from every cell (rows can mix keys)
  l <- system2("sh", c("-c", shQuote(sprintf("gzip -cd %s | grep -E '^!Sample_(geo_accession|title|characteristics_%s)'", shQuote(P(file)), ch))), stdout = TRUE)
  sp <- lapply(strsplit(l, "\t"), function(x) trimws(gsub('"', "", x[-1])))
  key <- sub("^!Sample_", "", vapply(strsplit(l, "\t"), `[`, "", 1))
  d <- data.table(gsm = sp[[which(key == "geo_accession")]], title = sp[[which(key == "title")]])
  long <- rbindlist(lapply(which(startsWith(key, "characteristics")), function(i) {
    v <- sp[[i]]; has <- grepl(":", v)
    data.table(gsm = d$gsm[has], k = make.names(trimws(sub(":.*$", "", v[has]))), val = trimws(sub("^[^:]*:", "", v[has])))
  }))
  stopifnot(!anyDuplicated(long[, .(gsm, k)]))
  merge(d, dcast(long, gsm ~ k, value.var = "val"), by = "gsm", all.x = TRUE, sort = FALSE)
}
fit_comparison <- function(id, family, role, dataset, tissue, contrast, expr, dat, covars, kind = c("array", "counts", "normcounts"), nrand = 1000) {
  kind <- match.arg(kind)
  dat <- as.data.table(dat); stopifnot(all(dat$sample %in% colnames(expr)))
  dat <- dat[complete.cases(dat[, c("disease", covars), with = FALSE])]
  for (v in covars) if (is.character(dat[[v]])) dat[[v]] <- factor(dat[[v]])
  covars <- covars[vapply(covars, function(v) length(unique(dat[[v]])) > 1, TRUE)]
  e <- expr[, dat$sample, drop = FALSE]
  rhs <- paste("disease", if (length(covars)) paste("+", paste(covars, collapse = " + ")) else "")
  X <- model.matrix(as.formula(paste("~", rhs)), dat); stopifnot(qr(X)$rank == ncol(X))
  if (kind == "counts") {
    y <- DGEList(round(e)); keep <- filterByExpr(y, group = dat$disease); y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
    loge <- cpm(y, log = TRUE, prior.count = 2); fit <- eBayes(lmFit(voom(y, X), X))
  } else if (kind == "normcounts") {
    keep <- filterByExpr(DGEList(round(e)), group = dat$disease); loge <- log2(e[keep, , drop = FALSE] + 1)
    fit <- eBayes(lmFit(loge, X), trend = TRUE)
  } else { loge <- e; fit <- eBayes(lmFit(loge, X)) }
  n <- nrow(loge); R <- singscore::rankGenes(loge); idx <- which(rownames(loge) %in% axisA)
  score_of <- function(ii) singscore::simpleScore(R, upSet = rownames(R)[ii])$TotalScore  # official singscore (Foroutan 2018)
  eff <- function(s) { d <- copy(dat); d$score_z <- as.numeric(scale(s)); f <- lm(as.formula(paste("score_z ~", rhs)), d)
    co <- summary(f)$coefficients["disease", ]; list(est = co[[1]], se = co[[2]], df = f$df.residual, ci = confint(f)["disease", ], p = co[[4]]) }
  obs <- eff(score_of(idx)); cam <- cameraPR(fit$t[, "disease"], index = list(AxisA = idx))
  mu <- rowMeans(loge); bins <- cut(rank(mu, ties.method = "first"), 20, labels = FALSE)
  need <- table(bins[idx]); rest <- setdiff(seq_len(n), idx); pool <- split(rest, bins[rest])
  null <- vapply(seq_len(nrand), function(b) eff(score_of(unlist(lapply(names(need), function(k) sample(pool[[k]], need[[k]])))))$est, 0)
  emp <- min(1, 2 * min((1 + sum(null >= obs$est)) / (nrand + 1), (1 + sum(null <= obs$est)) / (nrand + 1)))
  data.table(id, family, role, dataset, tissue, contrast, n_case = sum(dat$disease == 1), n_control = sum(dat$disease == 0),
             axisA_covered = length(idx), genes_ranked = n, covariates = paste(covars, collapse = ";"),
             effect_SD = obs$est, se = obs$se, df = obs$df, ci_low = obs$ci[[1]], ci_high = obs$ci[[2]], p = obs$p,
             camera_direction = cam$Direction, camera_p = cam$PValue, random_empirical_p = emp)
}
