# CollecTRI network built exactly as decoupleR 2.16.0 get_collectri(split_complexes = FALSE), from the OmniPath interactions
# endpoint downloaded on 2026-10-19 (OmnipathR 3.18.4 failed in unnest_evidences). TF-miRNA interactions are not added.
# Run from project root: R_LIBS=_workspace/r_libs_decoupler_20261019:_workspace/r_libs _workspace/envs/r/bin/Rscript 03_candidate_analysis/SingleCell_pseudobulk_supplements_20261019/code/00_collectri.R
suppressPackageStartupMessages({ library(data.table); library(dplyr); library(stringr) })
D <- "03_candidate_analysis/SingleCell_pseudobulk_supplements_20261019/source_data/"
x <- fread(paste0(D, "omnipath_collectri_raw_20261019.tsv")); x[, `:=`(is_stimulation = as.integer(as.logical(is_stimulation)), is_inhibition = as.integer(as.logical(is_inhibition)))]
cols <- c("source_genesymbol", "target_genesymbol", "is_stimulation", "is_inhibition")
ci <- as.data.frame(x[!str_detect(source, "COMPLEX"), ..cols]); cc <- as.data.frame(x[str_detect(source, "COMPLEX"), ..cols])
cc <- cc %>% mutate(source_genesymbol = case_when(str_detect(source_genesymbol, "JUN") | str_detect(source_genesymbol, "FOS") ~ "AP1",
                                                  str_detect(source_genesymbol, "REL") | str_detect(source_genesymbol, "NFKB") ~ "NFKB"))
net <- rbind(ci, cc) %>% distinct(source_genesymbol, target_genesymbol, .keep_all = TRUE) %>%
  mutate(weight = case_when(is_stimulation == 1 ~ 1, is_stimulation == 0 ~ -1)) %>% select(source_genesymbol, target_genesymbol, weight) %>%
  rename(source = source_genesymbol, target = target_genesymbol, mor = weight)
fwrite(net, paste0(D, "collectri_human_20261019.tsv"), sep = "\t")
cat("interactions", nrow(net), "| TFs", length(unique(net$source)), "| complexes kept as AP1/NFKB:", sum(net$source %in% c("AP1", "NFKB")), "\n"); print(table(net$mor))
print(net %>% filter(source %in% c("STAT3", "CEBPB", "FOS", "LYL1", "JUNB", "ELF1", "RREB1", "KLF6", "KLF7", "MEF2A")) %>% count(source))
