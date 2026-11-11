# Advice 2.4 / 2.5 / section 4 on saved donor pseudobulk counts.
# Run from project root:
#   LC_ALL=en_US.UTF-8 R_LIBS=_workspace/r_libs_decoupler_20261019:_workspace/r_libs _workspace/envs/r/bin/Rscript 03_candidate_analysis/SingleCell_pseudobulk_supplements_20261019/code/01_pseudobulk_tf.R
suppressPackageStartupMessages({ library(data.table); library(edgeR); library(limma); library(Matrix); library(decoupleR); library(AnnotationDbi); library(org.Hs.eg.db) })
O <- "03_candidate_analysis/SingleCell_pseudobulk_supplements_20261019"; PB <- "03_candidate_analysis/AxisA_pseudobulk_disease_20261014"
SYN <- "cell_localization/SYN01_Entrez_downstream_20261001/results"; T6 <- "01_current_analysis/01_storyline_current/07_virtual_ko_celloracle/results/T06E_SEAAD_independent_state_evaluation_20260827_candidate/inputs"
z <- function(x) (x - mean(x)) / sd(x)
net <- fread(file.path(O, "source_data/collectri_human_20261019.tsv"))
TOP10 <- c("STAT3", "FOS", "LYL1", "JUNB", "CEBPB", "ELF1", "RREB1", "KLF6", "KLF7", "MEF2A")
stopifnot(setequal(TOP10, fread("03_candidate_analysis/CellOracle_ranking_scale_20261019/results/all357_TF_ranking.tsv")[rank_shared_down <= 10]$TF))

# ---------------- inputs ----------------
lc <- fread(file.path(PB, "source_data/lung_myeloid_pseudobulk_counts.tsv.gz")); lm <- fread(file.path(PB, "source_data/lung_myeloid_donors.tsv"))
led <- unique(fread("01_current_analysis/01_storyline_current/06_single_cell_projection/results/GSE302339_axisA_official_ucell_20260903_candidate/tables/sample_ledger.tsv")[, .(subject_id, emphysema_group, fev1)])
stopifnot(!anyDuplicated(led$subject_id)); lm <- merge(lm, led, by = "subject_id", sort = FALSE)
L <- as.matrix(lc[, -1]); rownames(L) <- lc$gene_id; L <- L[, lm$subject_id]
lm[, `:=`(disease = factor(disease_group, levels = c("Control", "COPD")), age_z = z(age), sex = factor(sex), smoke = factor(smoke, levels = c("Never", "Former", "Current")),
          emph = factor(emphysema_group, levels = c("Control", "E-COPD", "NE-COPD")), fev1_z = z(as.numeric(fev1)))]
stopifnot(!anyNA(lm[, .(disease, age_z, sex, smoke, emph, fev1_z)]), all(lm[disease == "Control"]$emph == "Control"))
lung_axis <- readLines(file.path(SYN, "lung_axis_measured.txt"))

map <- fread(file.path(PB, "source_data/brain_symbol_to_ensembl.tsv")); brain_axis <- readLines(file.path(SYN, "brain_axis_measured.txt"))
adnc <- fread(file.path(O, "source_data/seaad_donor_ADNC.tsv"))
bm <- fread(file.path(SYN, "disease/brain_CPS_donor_data.tsv"))[Supertype == "Micro-PVM_2"]; setorder(bm, `Donor ID`)
bc <- fread("virtual_knockout/AxisA_disease_scale_check_20260925/source_data/brain_donor_pseudobulk_counts.tsv.gz"); stopifnot(identical(bc$gene, map$gene))
B2 <- as.matrix(bc[, -1]); rownames(B2) <- map$gene_id; B2 <- B2[, bm$`Donor ID`]
mg <- t(as.matrix(readMM(file.path(T6, "SEAAD_DFC_microglia_pseudobulk_counts.mtx.gz")))); gg <- readLines(file.path(T6, "SEAAD_DFC_microglia_pseudobulk_genes.txt")); dd <- readLines(file.path(T6, "SEAAD_DFC_microglia_pseudobulk_donors.txt"))
if (nrow(mg) != length(gg)) mg <- t(mg); stopifnot(nrow(mg) == length(gg), ncol(mg) == length(dd), identical(gg, map$gene))
rownames(mg) <- map$gene_id; colnames(mg) <- dd
mm <- fread(file.path(T6, "SEAAD_DFC_microglia_donor_metadata.tsv")); setnames(mm, c("donor_id", "cps_local", "pmi", "rin"), c("Donor ID", "cps_local", "pmi", "rin")); setorder(mm, `Donor ID`); BM <- mg[, mm$`Donor ID`]
prep_brain <- function(d) { d <- merge(copy(d), adnc, by = "Donor ID", all.x = TRUE, sort = TRUE); stopifnot(!anyNA(d$ADNC))
  d[, `:=`(CPS_Local_z = z(cps_local), age_z = z(age), PMI_z = z(pmi), RIN_z = z(rin), sex = factor(sex), study = factor(study),
           ADNC_bin = factor(fifelse(ADNC %in% c("Intermediate", "High"), "ADNC", "noADNC"), levels = c("noADNC", "ADNC")))]; d }
bm <- prep_brain(bm); mm <- prep_brain(mm)

# symbols for decoupleR
ens2sym <- function(ids) { s <- mapIds(org.Hs.eg.db, ids, "SYMBOL", "ENSEMBL", multiVals = "first"); s }
lung_sym <- ens2sym(rownames(L)); brain_sym <- setNames(map$gene, map$gene_id)

# ---------------- core runner ----------------
run <- function(organ, cls, model, counts, design, coef, axis, symbols, keep_cols = colnames(counts)) {
  counts <- counts[, keep_cols, drop = FALSE]; design <- design[keep_cols, , drop = FALSE]
  stopifnot(identical(colnames(counts), rownames(design)), qr(design)$rank == ncol(design))
  y <- DGEList(counts); keep <- filterByExpr(y, design = design); y <- normLibSizes(y[keep, , keep.lib.sizes = FALSE])
  members <- intersect(axis, rownames(y)); j <- match(coef, colnames(design)); stopifnot(!is.na(j))
  y <- estimateDisp(y, design); fr <- fry(y, index = list(AxisA = members), design = design, contrast = j)
  v <- voom(y, design); cm <- camera(v, index = list(AxisA = which(rownames(v) %in% members)), design = design, contrast = j, inter.gene.cor = NA)
  row <- data.table(organ, class = cls, model, coefficient = coef, n_donors = ncol(counts), axis_after_filter = length(members),
                    fry_direction = fr$Direction, fry_p = fr$PValue, camera_direction = cm$Direction, camera_p = cm$PValue, camera_corr = cm$Correlation)
  # decoupleR ulm on logCPM with gene symbols
  lc <- cpm(y, log = TRUE, prior.count = 2); sym <- symbols[rownames(lc)]; ok <- !is.na(sym) & !duplicated(sym); lc <- lc[ok, ]; rownames(lc) <- sym[ok]
  axis_sym <- unname(symbols[members]); axis_sym <- axis_sym[!is.na(axis_sym)]
  score <- singscore::simpleScore(singscore::rankGenes(lc), upSet = intersect(axis_sym, rownames(lc)))$TotalScore
  tf_one <- function(nt, label) {
    act <- run_ulm(mat = lc, network = as.data.frame(nt), .source = "source", .target = "target", .mor = "mor", minsize = 5)
    A <- as.matrix(dcast(as.data.table(act), source ~ condition, value.var = "score"), rownames = "source")[, colnames(lc), drop = FALSE]
    f <- eBayes(lmFit(A, design)); tt <- as.data.table(topTable(f, coef = j, number = Inf, sort.by = "none"), keep.rownames = "TF")
    tt[, r_with_axisA_score := apply(A[TF, , drop = FALSE], 1, function(a) cor(a, score))]
    ntar <- as.data.table(nt)[target %in% rownames(lc), .N, by = source]; tt[, n_targets_measured := ntar$N[match(TF, ntar$source)]]
    tt[, `:=`(organ = organ, class = cls, model = model, network = label)]
    tt[TF %in% TOP10, FDR_top10 := p.adjust(P.Value, "BH")][, FDR_all := p.adjust(P.Value, "BH")][] }
  tf <- rbind(tf_one(net, "CollecTRI"), tf_one(net[!target %in% axis_sym], "CollecTRI minus Axis A targets"))
  list(row = row, tf = tf)
}

# ---------------- models ----------------
dl1 <- model.matrix(~ disease + age_z + sex + smoke, lm); rownames(dl1) <- lm$subject_id
dlf <- model.matrix(~ fev1_z + age_z + sex + smoke, lm); rownames(dlf) <- lm$subject_id
bdes <- function(d, f) { x <- model.matrix(f, d); rownames(x) <- d$`Donor ID`; x }
R <- list()
R$lung_primary <- run("lung", "Myeloid", "primary: ~ disease + age + sex + smoking", L, dl1, "diseaseCOPD", lung_axis, lung_sym)
R$lung_fev1 <- run("lung", "Myeloid", "FEV1 % pred (per SD) + age + sex + smoking", L, dlf, "fev1_z", lung_axis, lung_sym)
# emphysema: reparametrise so that one coefficient is E-COPD minus NE-COPD (direct test), others vs control
lm[, emph2 := relevel(emph, ref = "NE-COPD")]; de1 <- model.matrix(~ emph2 + age_z + sex + smoke, lm); rownames(de1) <- lm$subject_id
R$lung_E_vs_NE <- run("lung", "Myeloid", "E-COPD vs NE-COPD (3-group model) + age + sex + smoking", L, de1, "emph2E-COPD", lung_axis, lung_sym)
de2 <- model.matrix(~ emph + age_z + sex + smoke, lm); rownames(de2) <- lm$subject_id
R$lung_E_vs_C <- run("lung", "Myeloid", "E-COPD vs control (3-group model)", L, de2, "emphE-COPD", lung_axis, lung_sym)
R$lung_NE_vs_C <- run("lung", "Myeloid", "NE-COPD vs control (3-group model)", L, de2, "emphNE-COPD", lung_axis, lung_sym)
fcps <- ~ CPS_Local_z + age_z + sex + PMI_z + RIN_z + study; fad <- ~ ADNC_bin + age_z + sex + PMI_z + RIN_z + study
for (cls in c("Microglia 2", "Microglia overall")) {
  d <- if (cls == "Microglia 2") bm else mm; X <- if (cls == "Microglia 2") B2 else BM; ncol_n <- d$n_nuclei
  k50 <- d$`Donor ID`[ncol_n >= 50]
  R[[paste(cls, "CPS")]] <- run("brain", cls, "primary: CPS_Local (per SD) + age + sex + PMI + RIN + study", X, bdes(d, fcps), "CPS_Local_z", brain_axis, brain_sym)
  R[[paste(cls, "ADNC")]] <- run("brain", cls, "ADNC intermediate/high vs not/low + covariates", X, bdes(d, fad), "ADNC_binADNC", brain_axis, brain_sym)
  d50 <- d[`Donor ID` %in% k50]; d50[, `:=`(CPS_Local_z = z(cps_local), age_z = z(age), PMI_z = z(pmi), RIN_z = z(rin))]
  R[[paste(cls, "CPS50")]] <- run("brain", cls, sprintf("sensitivity: >=50 nuclei (%d donors)", nrow(d50)), X, bdes(d50, fcps), "CPS_Local_z", brain_axis, brain_sym, keep_cols = d50$`Donor ID`)
}
rows <- rbindlist(lapply(R, `[[`, "row")); tfs <- rbindlist(lapply(R, `[[`, "tf"))
prim <- rows$model %like% "^primary"; rows[prim, `:=`(fry_BH_within_organ = p.adjust(fry_p, "BH"), camera_BH_within_organ = p.adjust(camera_p, "BH")), by = organ]
ref <- fread(file.path(PB, "results/fry_summary.tsv"))
chk <- c(lung = abs(rows[organ == "lung" & model %like% "^primary"]$fry_p - ref[organ == "lung" & model %like% "^primary"]$fry_p),
         brain = abs(rows[class == "Microglia 2" & model %like% "^primary"]$fry_p - ref[organ == "brain"]$fry_p)); print(chk); stopifnot(all(chk < 1e-10))
fwrite(rows, file.path(O, "results/pseudobulk_fry_camera.tsv"), sep = "\t"); fwrite(tfs, file.path(O, "results/tf_activity_all.tsv.gz"), sep = "\t")
print(rows[, .(organ, class, model = substr(model, 1, 45), n_donors, axis_after_filter, fry_direction, fry_p = signif(fry_p, 3), camera_direction, camera_p = signif(camera_p, 3), camera_corr = round(camera_corr, 3))])
top <- tfs[TF %in% TOP10 & (model %like% "^primary")]
print(top[, .(organ, class, network, TF, n_targets_measured, t = round(t, 2), P = signif(P.Value, 2), FDR_top10 = signif(FDR_top10, 2), r_axisA = round(r_with_axisA_score, 2))][order(organ, class, network, P)])
print(tfs[model %like% "^primary" & network == "CollecTRI", .(n_TF = .N, n_FDR_all_lt_0.05 = sum(FDR_all < 0.05)), by = .(organ, class)])
writeLines(capture.output(sessionInfo()), file.path(O, "qa/sessionInfo_01.txt"))
