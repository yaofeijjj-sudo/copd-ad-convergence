source("R/runtime.R")
theme_label <- function(term_name, term_id = term_name) {
  nm <- toupper(paste(term_id, term_name))
  case_when(
    str_detect(
      nm,
      "IMMUNE|INFLAM|CYTOKINE|INTERLEUKIN|INTERFERON|NF.?KB|TOLL|TLR|COMPLEMENT|LEUKOCYTE|MYELOID|NEUTROPHIL|MACROPHAGE|MICROGLIA|PHAGOCYT|B_CELL|T_CELL|TNF|ANTIGEN|HLA|CD40|FCGR|INTEGRIN"
    ) ~ "immune_inflammatory_myeloid",
    str_detect(
      nm,
      "MITOCHON|OXIDATIVE[ _]PHOSPHORYLATION|RESPIRATORY|ELECTRON|ATP|NADH|TCA|CITRIC|METABOL|GLYCOL|FATTY|LIPID|CHOLESTEROL|OXIDATIVE[ _]STRESS|HYPOXIA|ROS"
    ) ~ "mitochondrial_metabolic_oxidative",
    str_detect(
      nm,
      "PROTEIN|UBIQUITIN|PROTEASOME|CHAPERONE|FOLD|TRANSLATION|RIBOSOM|RNA|MRNA|SPLICE|EIF"
    ) ~ "protein_rna_translation",
    str_detect(
      nm,
      "LYSOSOM|AUTOPHAG|ENDOSOM|VESICLE|PHAGOSOM|VACUOL|ENDOCYT"
    ) ~ "lysosome_autophagy_vesicle",
    str_detect(
      nm,
      "SYNAP|NEURON|AXON|DENDR|NEUROTRANSMITTER|VESICLE CYCLE"
    ) ~ "neuronal_synaptic",
    str_detect(
      nm,
      "ACTIN|CYTOSKELETON|MIGRATION|ADHESION|EXTRACELLULAR[ _]MATRIX|COLLAGEN|VASCULAR|PLATELET|COAGULATION|ENDOTHELIAL"
    ) ~ "cytoskeleton_ecm_vascular",
    str_detect(
      nm,
      "STRESS|DNA DAMAGE|APOPT|SENESC|P53|RESPONSE TO|NECROSIS"
    ) ~ "cell_stress_death",
    TRUE ~ "other"
  )
}

complete_family_ora <- function(
  term2gene, query_genes, universe_genes, source_name,
  min_size, max_size
) {
  universe_genes <- clean_symbols(universe_genes)
  query_genes <- intersect(clean_symbols(query_genes), universe_genes)

  term_stats <- term2gene |>
    transmute(
      term_id = as.character(term_id),
      term_name = as.character(term_name),
      collection = as.character(collection),
      db_source = as.character(db_source),
      gene = toupper(trimws(as.character(gene)))
    ) |>
    filter(gene %in% universe_genes) |>
    distinct(term_id, term_name, collection, db_source, gene) |>
    group_by(term_id, term_name, collection, db_source) |>
    summarise(
      bg_count = n_distinct(gene),
      overlap_genes = paste(
        sort(intersect(unique(gene), query_genes)),
        collapse = "/"
      ),
      .groups = "drop"
    ) |>
    mutate(
      gene_count = if_else(
        nzchar(overlap_genes),
        as.integer(str_count(overlap_genes, fixed("/")) + 1L),
        0L
      )
    ) |>
    filter(bg_count >= min_size, bg_count <= max_size)

  n_universe <- length(universe_genes)
  n_query <- length(query_genes)

  term_stats |>
    mutate(
      source = source_name,
      query_size = n_query,
      universe_size = n_universe,
      fisher_a_overlap = gene_count,
      fisher_b_query_only = n_query - gene_count,
      fisher_c_term_only = bg_count - gene_count,
      fisher_d_neither = (
        n_universe - gene_count - fisher_b_query_only - fisher_c_term_only
      ),
      odds_ratio = if_else(
        fisher_b_query_only * fisher_c_term_only == 0,
        Inf,
        fisher_a_overlap * fisher_d_neither /
          (fisher_b_query_only * fisher_c_term_only)
      ),
      gene_ratio = paste0(gene_count, "/", n_query),
      bg_ratio = paste0(bg_count, "/", n_universe),
      p_value = phyper(
        gene_count - 1L, bg_count, n_universe - bg_count, n_query,
        lower.tail = FALSE
      ),
      zero_query_hit = gene_count == 0L,
      theme = theme_label(term_name, term_id)
    ) |>
    mutate(
      fdr_bh_complete_family = p.adjust(p_value, method = "BH"),
      neg_log10_fdr_complete_family = -log10(pmax(
        fdr_bh_complete_family, .Machine$double.xmin
      ))
    ) |>
    arrange(
      fdr_bh_complete_family, p_value, desc(gene_count), term_id, term_name
    ) |>
    mutate(rank_within_axis_database_background = row_number()) |>
    select(
      source, collection, db_source, term_id, term_name, theme,
      gene_count, zero_query_hit, overlap_genes, odds_ratio,
      fisher_a_overlap, fisher_b_query_only,
      fisher_c_term_only, fisher_d_neither,
      p_value, fdr_bh_complete_family,
      neg_log10_fdr_complete_family,
      gene_ratio, bg_ratio, query_size, universe_size, bg_count,
      rank_within_axis_database_background
    )
}
