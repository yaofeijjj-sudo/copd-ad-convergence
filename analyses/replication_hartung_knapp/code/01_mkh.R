# Sensitivity: modified Knapp-Hartung (mKH; Roever, Knapp & Friede 2015) for the six-cohort AD pool.
# Input: frozen per-cohort effects after the 10-19 sample exclusion. metafor test = "adhoc" implements the ad hoc
# correction (variance scaling factor truncated at 1), so the CI is never narrower than the unscaled t-based interval.
suppressPackageStartupMessages({library(data.table);library(metafor)})
ROOT <- '.'
IN <- file.path(ROOT,'03_candidate_analysis/AxisA_replication_after_exclusion_20261019/results')
OUT <- file.path(ROOT,'03_candidate_analysis/AxisA_meta_mKH_20261028/results')
b <- fread(file.path(IN,'by_cohort.tsv'))
ad <- b[cohort %in% c('GSE84422','GSE118553','GSE313307','GSE140841','GSE48350','GSE36980')]
stopifnot(nrow(ad)==6)
fit <- function(label, method, test) {
  f <- rma(yi=ad$effect_SD, sei=ad$se, method=method, test=test)
  data.table(analysis=label, effect_SD=as.numeric(f$b), ci_low=f$ci.lb, ci_high=f$ci.ub, p=f$pval, I2=f$I2, tau2=f$tau2)
}
res <- rbind(fit('DL (primary)','DL','z'), fit('REML + Hartung-Knapp','REML','knha'),
             fit('REML + modified Hartung-Knapp (mKH)','REML','adhoc'))
# check against the frozen pooled_sensitivity.tsv
ref <- fread(file.path(IN,'pooled_sensitivity.tsv'))
for (a in c('DL (primary)','REML + Hartung-Knapp'))
  stopifnot(abs(res[analysis==a]$ci_low-ref[analysis==a]$ci_low)<1e-9, abs(res[analysis==a]$p-ref[analysis==a]$p)<1e-9)
# Knapp-Hartung scaling factor q (mKH uses max(q, 1))
f0 <- rma(yi=ad$effect_SD, sei=ad$se, method='REML'); w <- 1/(ad$se^2+f0$tau2)
q <- sum(w*(ad$effect_SD-as.numeric(f0$b))^2)/(nrow(ad)-1)
res[, KH_scaling_q := c(NA, q, q)]
fwrite(res, file.path(OUT,'pooled_mKH.tsv'), sep='\t'); print(res); cat('metafor', as.character(packageVersion('metafor')), '\n')
