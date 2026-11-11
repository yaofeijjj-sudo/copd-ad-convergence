# Advice 3.3: GSE150910 diagnosis x processing batch (date) and x plate cross-tables (chronic HP vs control).
# Run from project root: R_LIBS=_workspace/r_libs _workspace/envs/r/bin/Rscript 03_candidate_analysis/AxisA_auxiliary_checks_20261019/code/04_chp_batch_crosstab.R
suppressPackageStartupMessages(library(data.table))
m <- fread("03_candidate_analysis/AxisA_disease_specificity_redesign_20261003/metadata/GSE150910_sample_metadata.tsv")
dx <- grep("diagnosis|disease", names(m), value = TRUE, ignore.case = TRUE)[1]; bt <- grep("^batch", names(m), value = TRUE, ignore.case = TRUE)[1]; pl <- grep("plate", names(m), value = TRUE, ignore.case = TRUE)[1]
cat("columns:", dx, bt, pl, "\n"); print(table(m[[dx]]))
s <- m[tolower(m[[dx]]) %in% c("chp", "control")]
tb <- as.data.table(table(batch = s[[bt]], diagnosis = s[[dx]])); tw <- dcast(tb, batch ~ diagnosis, value.var = "N")
tw[, both_present := get(names(tw)[2]) > 0 & get(names(tw)[3]) > 0]; print(tw)
cat("batches containing both CHP and control:", sum(tw$both_present), "of", nrow(tw), "\n")
tp <- dcast(as.data.table(table(plate = s[[pl]], diagnosis = s[[dx]])), plate ~ diagnosis, value.var = "N"); print(tp)
fwrite(tw, "03_candidate_analysis/AxisA_auxiliary_checks_20261019/results/GSE150910_diagnosis_by_batch.tsv", sep = "\t")
fwrite(tp, "03_candidate_analysis/AxisA_auxiliary_checks_20261019/results/GSE150910_diagnosis_by_plate.tsv", sep = "\t")
