if (dir.exists("env/library")) .libPaths(c("env/library", .libPaths()))
params <- yaml::read_yaml("config/params.yml")
options(stringsAsFactors = FALSE)
