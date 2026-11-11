source("R/runtime.R")
suppressPackageStartupMessages(library(data.table))
p <- "results/single_cell/preprocessing/lung"
q <- fread(file.path(p, "objects/cell_qc_all.tsv.gz"))[passes_QC == TRUE]
a <- fread(file.path(p, "objects/celltypist_all.tsv.gz"))
l <- fread(file.path(p, "tables/sample_ledger.tsv"))[matrix_available == TRUE]
d <- merge(q, a, by = c("gsm", "barcode"), sort = FALSE)
fields <- c("gsm", "subject_id", "disease_group", "emphysema_group", "lobe", "age", "sex", "smoke", "gold", "fev1", "fev1_fvc", "lobe_laa_percent", "packyear", "steroidstype", "steroids")
d <- merge(d, l[, ..fields], by = "gsm", sort = FALSE, all.x = TRUE)
d[, cell_id := paste(gsm, barcode, sep = "|")]
stopifnot(nrow(d) == 211464L, !anyDuplicated(d$cell_id))
fwrite(d, file.path(p, "objects/axisA_cell_scores.tsv.gz"), sep = "\t")
