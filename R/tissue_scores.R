source("R/runtime.R")
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(singscore)
})

axis_a <- fread("results/module_overlap/results/cross/AxisA_confirmed_genes.csv", colClasses = "character")$entrez_id
stopifnot(length(axis_a) == params$axis_a_size)

history <- fread(cmd = "gzip -cd data/reference/gene_identity/references/gene_history.tsv.gz", colClasses = "character")
setnames(history, 1, "tax_id")
history <- history[tax_id == "9606" & GeneID != "-"]
replacement <- setNames(history$GeneID, history$Discontinued_GeneID)

# Only explicit NCBI replacement records are followed; symbols are never used to rescue an ID.
current_id <- function(x) {
  for (i in 1:5) {
    hit <- !is.na(x) & x %in% names(replacement)
    if (!any(hit)) break
    x[hit] <- replacement[x[hit]]
  }
  x
}

read_series_matrix <- function(file) {
  cmd <- paste(
    "gzip -cd", shQuote(file),
    "| sed -n '/^!series_matrix_table_begin/,$p' | sed '1d;/^!series_matrix_table_end/,$d'"
  )
  tab <- fread(cmd = cmd, check.names = FALSE)
  m <- as.matrix(tab[, -1L, with = FALSE])
  storage.mode(m) <- "numeric"
  rownames(m) <- as.character(tab[[1L]])
  m
}

# One feature per Entrez ID: the one with the largest MAD across samples.
keep_max_mad <- function(m, gene) {
  ok <- !is.na(gene)
  m <- m[ok, , drop = FALSE]
  gene <- gene[ok]
  v <- apply(m, 1, mad, na.rm = TRUE)
  o <- order(gene, -v)
  o <- o[is.finite(v[o])]
  o <- o[!duplicated(gene[o])]
  out <- m[o, , drop = FALSE]
  rownames(out) <- gene[o]
  out[rowSums(is.finite(out)) == ncol(out), , drop = FALSE]
}

collapse_array <- function(m, annot, idcol = "Gene ID") {
  a <- if (grepl("\\.annot\\.gz$", annot)) {
    fread(
      cmd = paste("gzip -cd", shQuote(annot), "| sed -n '/^ID\t/,$p' | grep -v '^!platform_table_end'"),
      colClasses = "character", select = c("ID", idcol)
    )
  } else {
    fread(annot, colClasses = "character", select = c("ID", idcol))
  }
  setnames(a, c("probe", "gene"))
  # Probes annotated to more than one gene are dropped.
  a <- a[probe %in% rownames(m) & !is.na(gene) & nzchar(gene) & !grepl("/", gene)]
  keep_max_mad(m[a$probe, , drop = FALSE], current_id(a$gene))
}

# Keys that map to more than one Entrez ID are left unmapped.
unique_map <- function(keys, keytype) {
  mp <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = unique(keys), keytype = keytype, columns = "ENTREZID"))
  mp <- as.data.table(mp)[!is.na(ENTREZID)]
  setnames(mp, 1, "key")
  mp <- mp[, if (.N == 1) .SD, by = key]
  current_id(mp$ENTREZID[match(keys, mp$key)])
}

# GEO sample characteristics as one column per key; a key is read from every cell because rows can mix keys.
series_chars <- function(file, ch = "ch1") {
  pattern <- sprintf("'^!Sample_(geo_accession|title|characteristics_%s)'", ch)
  l <- system2("sh", c("-c", shQuote(paste("gzip -cd", shQuote(file), "| grep -E", pattern))), stdout = TRUE)
  fields <- strsplit(l, "\t")
  key <- sub("^!Sample_", "", vapply(fields, `[`, "", 1))
  values <- lapply(fields, function(x) trimws(gsub('"', "", x[-1])))
  d <- data.table(gsm = values[[which(key == "geo_accession")]], title = values[[which(key == "title")]])
  long <- rbindlist(lapply(which(startsWith(key, "characteristics")), function(i) {
    v <- values[[i]]
    has <- grepl(":", v)
    data.table(gsm = d$gsm[has], k = make.names(trimws(sub(":.*$", "", v[has]))), val = trimws(sub("^[^:]*:", "", v[has])))
  }))
  stopifnot(!anyDuplicated(long[, .(gsm, k)]))
  merge(d, dcast(long, gsm ~ k, value.var = "val"), by = "gsm", all.x = TRUE, sort = FALSE)
}

# Axis A singscore per sample, standardized within the dataset, regressed on disease and covariates.
# kind: 'array' (log-scale matrix), 'counts' (raw counts -> TMM log-CPM), 'normcounts' (submitter-normalized counts -> log2(x + 1)).
fit_score <- function(id, expr, dat, covars, kind = c("array", "counts", "normcounts")) {
  kind <- match.arg(kind)
  dat <- as.data.table(dat)
  stopifnot(all(dat$sample %in% colnames(expr)))
  dat <- dat[complete.cases(dat[, c("disease", covars), with = FALSE])]
  for (v in covars) if (is.character(dat[[v]])) dat[[v]] <- factor(dat[[v]])
  covars <- covars[vapply(covars, function(v) length(unique(dat[[v]])) > 1, TRUE)]
  e <- expr[, dat$sample, drop = FALSE]
  rhs <- paste("disease", if (length(covars)) paste("+", paste(covars, collapse = " + ")) else "")
  X <- model.matrix(as.formula(paste("~", rhs)), dat)
  stopifnot(qr(X)$rank == ncol(X))
  loge <- switch(kind,
    counts = {
      y <- DGEList(round(e))
      keep <- filterByExpr(y, group = dat$disease)
      y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
      cpm(y, log = TRUE, prior.count = 2)
    },
    normcounts = {
      keep <- filterByExpr(DGEList(round(e)), group = dat$disease)
      log2(e[keep, , drop = FALSE] + 1)
    },
    array = e
  )
  ranks <- rankGenes(loge)
  measured <- rownames(loge)[rownames(loge) %in% axis_a]
  dat$score_z <- as.numeric(scale(simpleScore(ranks, upSet = measured)$TotalScore))
  f <- lm(as.formula(paste("score_z ~", rhs)), dat)
  co <- summary(f)$coefficients["disease", ]
  ci <- confint(f)["disease", ]
  data.table(
    id = id, n_case = sum(dat$disease == 1), n_control = sum(dat$disease == 0),
    axis_a_measured = length(measured), genes_ranked = nrow(loge), covariates = paste(covars, collapse = ";"),
    effect_SD = co[[1]], se = co[[2]], df = f$df.residual, ci_low = ci[[1]], ci_high = ci[[2]], p = co[[4]]
  )
}
