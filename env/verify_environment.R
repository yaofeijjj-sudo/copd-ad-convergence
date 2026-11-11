expected <- read.delim("env/r-versions.tsv", stringsAsFactors = FALSE)
for (i in seq_len(nrow(expected))) stopifnot(as.character(packageVersion(expected$package[i])) == expected$version[i])
cat(nrow(expected), "R package versions matched\n")
