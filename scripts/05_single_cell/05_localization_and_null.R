source("R/runtime.R")
source("scripts/05_single_cell/common.R")
OUT <- file.path(ROOT, "results/single_cell/results/models")
for (d in c("tables", "qa", "objects")) dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)
axis <- AXIS
matched_null <- function(gain, profiles, props, axis, target, cohort, prefix) {
  m <- length(axis)
  den <- m * RMAX - m * (m + 1) / 2
  profiles[, mean_ucell := colSums(gain[axis, , drop = FALSE]) / den]
  model <- fit_model(profiles, sort(unique(profiles$annotation)))
  res <- target_contrasts(model, target, cohort, prefix)
  keys <- paste(model$data$donor, model$data$annotation, sep = "|||")
  idx <- match(keys, paste(profiles$donor, profiles$annotation, sep = "|||"))
  stopifnot(!anyNA(idx))
  gc <- drop(gain[, idx, drop = FALSE] %*% res$main$profile_weights)
  names(gc) <- rownames(gain)
  obs <- sum(gc[axis]) / den
  stopifnot(abs(obs - res$main$row$effect) < 1e-12)
  props <- copy(props)
  props[, axis_member := gene %chin% axis]
  props[, mean_bin := pmax(1L, pmin(8L, ceiling(frank(mean_expression, ties.method = "average") / .N * 8)))]
  props[, variance_bin := pmax(1L, pmin(8L, ceiling(frank(expression_variance, ties.method = "average") / .N * 8)))]
  props[, `:=`(z_mean = as.numeric(scale(log1p(mean_expression))), z_variance = as.numeric(scale(log1p(expression_variance))))]
  targetprop <- props[axis_member == TRUE]
  bg <- props[axis_member == FALSE]
  stopifnot(nrow(targetprop) == m)
  pools <- split(bg$gene, paste(bg$mean_bin, bg$variance_bin, sep = "|"))
  nearest <- lapply(seq_len(m), function(i) bg$gene[order((bg$z_mean - targetprop$z_mean[i])^2 + (bg$z_variance - targetprop$z_variance[i])^2)])
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(params$matched_null_seed)
  sets <- vector("list", params$matched_sets)
  matched <- vector("list", params$matched_sets)
  null <- numeric(params$matched_sets)
  for (b in seq_len(params$matched_sets)) {
    chosen <- character()
    rules <- character(m)
    for (i in seq_len(m)) {
      pool <- setdiff(pools[[paste(targetprop$mean_bin[i], targetprop$variance_bin[i], sep = "|")]], chosen)
      if (length(pool)) {
        pick <- sample(pool, 1)
        rule <- "same_mean_and_variance_bin"
      } else {
        pick <- setdiff(nearest[[i]], chosen)[1]
        rule <- "nearest_unused_in_standardized_mean_variance"
      }
      chosen <- c(chosen, pick)
      rules[i] <- rule
    }
    sets[[b]] <- chosen
    null[b] <- sum(gc[chosen]) / den
    matched[[b]] <- data.table(random_signature = sprintf("RandomMatched%04d", b), target_gene = targetprop$gene, random_gene = chosen, match_rule = rules)
  }
  matched <- rbindlist(matched)
  stopifnot(!any(matched$random_gene %chin% axis), all(vapply(sets, uniqueN, integer(1)) == m))
  draws <- data.table(cohort = cohort, random_signature = sprintf("RandomMatched%04d", seq_len(params$matched_sets)), effect = null, observed_effect = obs, ge_observed = null >= obs)
  prefix0 <- sub("_broad$", "", prefix)
  summary <- data.table(cohort = cohort, target = target, observed_effect = obs, measured_axis_genes = m, universe_size = nrow(props), B = params$matched_sets, seed = params$matched_null_seed, exceedances = sum(null >= obs), empirical_p = (1 + sum(null >= obs)) / (params$matched_sets + 1), null_median = median(null), null_q95 = quantile(null, .95), null_max = max(null))
  write_tsv(draws, paste0(prefix0, "_matched_null_draws.tsv"))
  write_tsv(summary, paste0(prefix0, "_matched_null_summary.tsv"))
  fwrite(matched, file.path(OUT, "tables", paste0(prefix0, "_matched_null_gene_sets.tsv.gz")), sep = "\t")
  write_tsv(props, paste0(prefix0, "_matching_properties.tsv"))
  balance <- merge(matched, props[, .(gene, mean_expression, expression_variance, detection_rate)], by.x = "random_gene", by.y = "gene")[, .(mean_expression = mean(mean_expression), expression_variance = mean(expression_variance), detection_rate = mean(detection_rate)), by = random_signature]
  balance[, `:=`(axis_mean_expression = mean(targetprop$mean_expression), axis_expression_variance = mean(targetprop$expression_variance), axis_detection_rate = mean(targetprop$detection_rate))]
  write_tsv(balance, paste0(prefix0, "_matched_null_balance.tsv"))
  write_tsv(matched[, .N, by = match_rule], paste0(prefix0, "_matched_null_match_rules.tsv"))
  qa <- rbindlist(lapply(c(1, 2500, params$matched_sets), function(b) {
    z <- copy(profiles)
    z[, mean_ucell := colSums(gain[sets[[b]], , drop = FALSE]) / den]
    mod <- fit_model(z, model$types)
    data.table(random_draw = b, error = abs(sum(res$main$cvec * coef(mod$fit)) - null[b]))
  }))
  stopifnot(max(qa$error) < 1e-12)
  write_tsv(qa, paste0(prefix0, "_null_numeric_QA.tsv"))
  saveRDS(list(model = model, primary = res$main), file.path(OUT, "objects", paste0(prefix0, "_models.rds")))
  print(summary)
  res$main$row
}
lung <- readRDS(file.path(ROOT, "results/single_cell/results/lung/objects/lung_rank_aggregate.rds"))
if (file.exists(file.path(OUT, "tables/lung_matched_null_summary.tsv"))) lres <- fread(file.path(OUT, "tables/lung_broad_primary_contrast.tsv")) else lres <- matched_null(lung$mean_gain, copy(lung$profiles), lung$gene_properties, axis, "Myeloid", "GSE302339", "lung_broad")
if ("lung-only" %in% commandArgs(trailingOnly = TRUE)) quit(save = "no", status = 0)
bwork <- file.path(ROOT, "results/single_cell/results/brain_work")
a <- readRDS(file.path(bwork, "broad_ucell/profile_gene_rank_aggregate.rds"))
p <- copy(a$profiles)
p[, n_cells := a$profile_n_cells]
mp <- fread(file.path(ROOT, "resources/brain_subclass_to_group.tsv"))
p[, group := mp$group[match(annotation, mp$subclass)]]
stopifnot(!anyNA(p$group))
g <- p[, .(n_cells = sum(n_cells)), by = .(donor, annotation = group)]
setorder(g, donor, annotation)
p[, group_index := match(paste(donor, group), paste(g$donor, g$annotation))]
A <- sparseMatrix(i = seq_len(nrow(p)), j = p$group_index, x = p$n_cells / g$n_cells[p$group_index], dims = c(nrow(p), nrow(g)))
gain <- RMAX - as.matrix(a$mean_truncated_rank %*% A)
colnames(gain) <- paste(g$donor, g$annotation, sep = "|||")
baxis <- a$axis_genes
bprop <- fread(file.path(bwork, "broad_ucell/gene_properties.tsv"))
bres <- matched_null(gain, g, bprop, baxis, "Immune", "SEA-AD", "brain_broad")
h <- rbind(lres, bres)
h[, Holm_two_organs := p.adjust(p, "holm")]
write_tsv(h, "harmonized_primary_localization.tsv")
sub <- readRDS(file.path(bwork, "subtype_ucell/profile_gene_rank_aggregate.rds"))
s <- copy(sub$profiles)
s[, n_cells := sub$profile_n_cells]
m <- length(sub$axis_genes)
den <- m * RMAX - m * (m + 1) / 2
s[, mean_ucell := colSums(RMAX - sub$mean_truncated_rank[sub$axis_genes, , drop = FALSE]) / den]
cov <- s[n_cells >= 20, .(n_donors = uniqueN(donor), n_cells = sum(n_cells)), by = annotation]
types <- sort(cov[n_donors >= 5]$annotation)
mod <- fit_model(s, types)
sr <- rbindlist(lapply(types, function(target) {
  w <- setNames(rep(-1 / (length(types) - 1), length(types)), types)
  w[target] <- 1
  contrast(mod, w, "SEA-AD", target, "subtype")$row
}))
sr[, BH_FDR := p.adjust(p, "BH")]
sr <- merge(sr, cov, by.x = "contrast", by.y = "annotation", suffixes = c("_model", "_target"))
sr[, coverage_tier := fifelse(n_donors_target >= 40, "confirmatory", fifelse(n_donors_target >= 16, "exploratory", "descriptive"))]
p[, mean_ucell := colSums(RMAX - a$mean_truncated_rank[a$axis_genes, , drop = FALSE]) / (length(a$axis_genes) * RMAX - length(a$axis_genes) * (length(a$axis_genes) + 1) / 2)]
write_tsv(p[n_cells >= 20], "supplementary24_brain_broad_donor_profiles.tsv")
write_tsv(sr, "brain_subtype_CR2.tsv")
write_tsv(s, "brain_subtype_donor_profiles.tsv")
writeLines(capture.output(sessionInfo()), file.path(OUT, "qa/model_sessionInfo.txt"))
