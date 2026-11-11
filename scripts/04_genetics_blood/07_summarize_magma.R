source("R/runtime.R")
.libPaths(c("env/library", .libPaths()))
library(data.table)
out <- "results/genetics_blood"
res <- rbindlist(lapply(c("AD", "COPD"), function(id) {
  d <- fread(file.path(out, "results/magma", paste0(id, "_axes.gsa.out")), skip = "VARIABLE", header = TRUE)
  d[, disease := id]
  d
}))
stopifnot(nrow(res) == 4, all(res$P > 0 & res$P <= 1))
res[, `:=`(holm_p_across_2x2_four_tests = p.adjust(P, "holm"), lo = BETA - 1.96 * SE, hi = BETA + 1.96 * SE, axis = ifelse(VARIABLE == "AxisA", "Axis A", "Axis B"), source_axis_size = ifelse(VARIABLE == "AxisA", 137L, 434L))]
fwrite(res, file.path(out, "results/magma_four_tests.csv"))
print(res[, .(disease, axis, NGENES, BETA, SE, P, holm_p_across_2x2_four_tests)])
