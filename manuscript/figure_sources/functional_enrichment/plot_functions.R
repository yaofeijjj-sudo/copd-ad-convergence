#!/usr/bin/env Rscript

# Standalone, side-by-side Sankey-bubble figure for the current Axis A/B ORA.
#
# Geometry is adapted from the public code accompanying:
# Wang et al. Longitudinal metabolome profiling from pregnancy through
# childhood and risk of neurodevelopmental disorders at age 10.
# Nature Communications, DOI: 10.1038/s41467-025-68115-3, Figure 5.
# Public source (Apache-2.0):
# https://github.com/tiwa1125/Longitudinal_metabolome/blob/
# 15128b41ed9beebd83576e0a0cfdb3d422fb5679/
# Enrichment/Pathway_Enrichment_NDD.R
#
# Adapted elements: neutral membership flows, coloured nodes, cumulative
# term-centre coordinates, and an adjacent aligned enrichment-bubble layer.
# The original ggsankey/cowplot geometry is translated to ggalluvial/patchwork.
# No ORA statistic, term, gene subset, or background is changed here.

# 2026-09-17 style adaptation: Arial 8 pt; original scientific encodings preserved.
arial_match <- systemfonts::match_fonts("Arial")
arial_path <- arial_match$path[[1]]
if (
  is.na(arial_path) ||
    !nzchar(arial_path) ||
    !file.exists(arial_path) ||
    !grepl("Arial", basename(arial_path), ignore.case = TRUE)
) {
  stop("Arial could not be resolved to an installed system font.")
}
font_family <- "Arial"

# Modern, high-contrast palette. Neutral links follow the published example;
# colour is reserved for axis identity, annotation database, and significance.
palette <- c(
  axis_a = "#F04438",
  axis_b = "#168BEA",
  go = "#2DD4BF",
  kegg = "#FBBF24",
  hallmark = "#A78BFA",
  ribbon = "#CBD5E1",
  panel = "#F8FAFC",
  panel_border = "#CBD5E1",
  text = "#111827",
  text_muted = "#64748B",
  nonsig = "#94A3B8"
)
q_colors <- c("#EEF2FF", "#C7D2FE", "#818CF8", "#6366F1", "#4338CA")

wrap_one <- function(x, width = 22L) {
  paste(strwrap(x, width = width), collapse = "\n")
}

database_key <- function(x) {
  ifelse(
    x == "GO BP", "GO",
    ifelse(x == "KEGG", "KEGG", "Hallmark")
  )
}

database_colour <- function(x) {
  unname(c(
    "GO BP" = palette[["go"]],
    "KEGG" = palette[["kegg"]],
    "Hallmark" = palette[["hallmark"]]
  )[x])
}

prepare_axis_data <- function(axis_name) {
  label_wrap_width <- if (axis_name == "Axis A") 18L else 22L
  axis_terms <- terms |>
    filter(axis == axis_name) |>
    arrange(display_order) |>
    mutate(
      database_short = database_key(database),
      label_modern = vapply(
        paste0(term_name, " (", database_short, ")"),
        wrap_one,
        character(1),
        width = label_wrap_width
      ),
      node_fill = database_colour(database)
    )

  axis_links <- links |>
    filter(axis == axis_name) |>
    arrange(gene_order, display_order) |>
    mutate(link_id = paste0(axis_name, "_", row_number()))

  gene_order <- axis_links |>
    distinct(gene, gene_order) |>
    arrange(gene_order)
  term_order <- axis_terms |>
    select(term_id, display_order) |>
    arrange(display_order)

  gene_counts <- axis_links |>
    count(gene, name = "displayed_memberships") |>
    right_join(gene_order, by = "gene") |>
    arrange(gene_order)
  term_counts <- axis_links |>
    count(term_id, name = "displayed_memberships") |>
    right_join(term_order, by = "term_id") |>
    arrange(display_order)

  total_links <- nrow(axis_links)
  gene_centres <- gene_counts |>
    mutate(
      ymax = total_links - lag(cumsum(displayed_memberships), default = 0),
      ymin = ymax - displayed_memberships,
      y = (ymin + ymax) / 2
    )
  term_centres <- term_counts |>
    mutate(
      ymax = total_links - lag(cumsum(displayed_memberships), default = 0),
      ymin = ymax - displayed_memberships,
      y = (ymin + ymax) / 2
    ) |>
    left_join(axis_terms, by = c("term_id", "display_order"))

  gene_ids <- paste0("gene:", gene_order$gene)
  term_ids <- paste0("term:", term_order$term_id)
  stratum_levels <- c(rev(gene_ids), rev(term_ids))

  gene_lodes <- axis_links |>
    transmute(
      link_id,
      x = 1,
      stratum = paste0("gene:", gene),
      node_type = "gene",
      database,
      weight = 1
    )
  term_lodes <- axis_links |>
    transmute(
      link_id,
      x = 2,
      stratum = paste0("term:", term_id),
      node_type = database,
      database,
      weight = 1
    )
  lodes <- bind_rows(gene_lodes, term_lodes) |>
    mutate(stratum = factor(stratum, levels = stratum_levels))

  list(
    terms = axis_terms,
    links = axis_links,
    lodes = lodes,
    gene_centres = gene_centres,
    term_centres = term_centres,
    total_links = total_links
  )
}

make_axis_plot <- function(
  axis_name,
  label_side = c("left", "right"),
  show_plot_legend = axis_name == "Axis A",
  term_node_x = 2,
  legend_layout = c("bottom", "right")
) {
  label_side <- match.arg(label_side)
  legend_layout <- match.arg(legend_layout)
  axis_data <- prepare_axis_data(axis_name)
  axis_data$lodes <- axis_data$lodes |>
    mutate(x = ifelse(node_type == "gene", 1, term_node_x))
  is_axis_a <- axis_name == "Axis A"
  axis_colour <- if (is_axis_a) palette[["axis_a"]] else palette[["axis_b"]]
  title_text <- if (is_axis_a) {
    "Axis A  |  disease-increased"
  } else {
    "Axis B  |  disease-decreased"
  }
  subtitle_text <- if (is_axis_a) {
    "COPD red \u2229 AD blue  \u00b7  135 genes"
  } else {
    "COPD turquoise \u2229 AD turquoise  \u00b7  418 genes"
  }
  ratio_limit <- if (is_axis_a) 0.60 else 0.06
  ratio_breaks <- if (is_axis_a) {
    c(0, 0.20, 0.40, 0.60)
  } else {
    c(0, 0.02, 0.04, 0.06)
  }

  labels_on_right <- identical(label_side, "right")
  term_node_left <- term_node_x - 0.0375
  term_node_right <- term_node_x + 0.0375
  bubble_left <- if (labels_on_right) {
    term_node_x + 1.02
  } else {
    term_node_x + 0.30
  }
  bubble_right <- if (labels_on_right) {
    term_node_x + 2.36
  } else {
    term_node_x + 1.38
  }
  bubble_axis_left <- if (labels_on_right) {
    term_node_x + 1.12
  } else {
    term_node_x + 0.37
  }
  bubble_axis_right <- if (labels_on_right) {
    term_node_x + 2.26
  } else {
    term_node_x + 1.31
  }
  term_label_x <- if (labels_on_right) {
    term_node_right + 0.0175
  } else {
    term_node_left - 0.0125
  }
  term_label_hjust <- if (labels_on_right) 0 else 1
  legend_on_right <- identical(legend_layout, "right")
  plot_xlim <- if (labels_on_right) {
    c(0.55, bubble_right + 0.04)
  } else {
    c(0.61, bubble_right + 0.05)
  }
  ratio_to_x <- function(x) {
    bubble_axis_left +
      (x / ratio_limit) * (bubble_axis_right - bubble_axis_left)
  }

  term_plot_data <- axis_data$term_centres |>
    mutate(
      bubble_x = ratio_to_x(gene_ratio),
      q_fill = ifelse(
        significant_primary,
        pmin(neg_log10_fdr, 8),
        NA_real_
      ),
      significance = factor(
        ifelse(
          significant_primary,
          "BH q < 0.05",
          "BH q \u2265 0.05"
        ),
        levels = c("BH q < 0.05", "BH q \u2265 0.05")
      )
    )
  tick_data <- data.frame(
    ratio = ratio_breaks,
    x = ratio_to_x(ratio_breaks),
    label = percent(ratio_breaks, accuracy = 1)
  )

  p <- ggplot(
    axis_data$lodes,
    aes(
      x = x,
      stratum = stratum,
      alluvium = link_id,
      y = weight
    )
  ) +
    annotate(
      "rect",
      xmin = bubble_left, xmax = bubble_right,
      ymin = 0, ymax = axis_data$total_links,
      fill = palette[["panel"]],
      colour = palette[["panel_border"]],
      linewidth = 0.42
    ) +
    geom_alluvium(
      fill = palette[["ribbon"]],
      alpha = 0.72,
      colour = NA,
      width = 0.075,
      knot.pos = 0.45,
      curve_type = "xspline",
      reverse = FALSE,
      show.legend = FALSE
    ) +
    geom_rect(
      data = axis_data$gene_centres,
      aes(
        xmin = 0.9625, xmax = 1.0375,
        ymin = ymin, ymax = ymax
      ),
      inherit.aes = FALSE,
      fill = axis_colour,
      colour = "white",
      linewidth = 0.34,
      show.legend = FALSE
    ) +
    geom_rect(
      data = filter(term_plot_data, database == "GO BP"),
      aes(
        xmin = term_node_left, xmax = term_node_right,
        ymin = ymin, ymax = ymax
      ),
      inherit.aes = FALSE,
      fill = palette[["go"]],
      colour = "white",
      linewidth = 0.34,
      show.legend = FALSE
    ) +
    geom_rect(
      data = filter(term_plot_data, database == "KEGG"),
      aes(
        xmin = term_node_left, xmax = term_node_right,
        ymin = ymin, ymax = ymax
      ),
      inherit.aes = FALSE,
      fill = palette[["kegg"]],
      colour = "white",
      linewidth = 0.34,
      show.legend = FALSE
    ) +
    geom_rect(
      data = filter(term_plot_data, database == "Hallmark"),
      aes(
        xmin = term_node_left, xmax = term_node_right,
        ymin = ymin, ymax = ymax
      ),
      inherit.aes = FALSE,
      fill = palette[["hallmark"]],
      colour = "white",
      linewidth = 0.34,
      show.legend = FALSE
    ) +
    geom_segment(
      data = term_plot_data,
      aes(
        x = bubble_left + 0.01,
        xend = bubble_x,
        y = y,
        yend = y
      ),
      inherit.aes = FALSE,
      colour = "#E2E8F0",
      linewidth = 0.48,
      lineend = "round"
    ) +
    geom_text(
      data = axis_data$gene_centres,
      aes(x = 0.94, y = y, label = gene),
      inherit.aes = FALSE,
      family = font_family,
      colour = palette[["text"]],
      size = 2.85,
      hjust = 1,
      lineheight = 0.92
    ) +
    geom_label(
      data = term_plot_data,
      aes(x = term_label_x, y = y, label = label_modern),
      inherit.aes = FALSE,
      family = font_family,
      colour = palette[["text"]],
      fill = alpha("white", 0.94),
      linewidth = 0,
      label.padding = grid::unit(0.65, "pt"),
      label.r = grid::unit(0, "pt"),
      size = 2.55,
      hjust = term_label_hjust,
      lineheight = 0.90
    ) +
    geom_point(
      data = term_plot_data,
      aes(
        x = bubble_x,
        y = y,
        size = gene_count,
        fill = q_fill,
        colour = significance
      ),
      inherit.aes = FALSE,
      shape = 21,
      stroke = 0.62,
      show.legend = show_plot_legend
    ) +
    geom_segment(
      data = tick_data,
      aes(
        x = x, xend = x,
        y = 0, yend = -0.65
      ),
      inherit.aes = FALSE,
      colour = palette[["text"]],
      linewidth = 0.38
    ) +
    geom_text(
      data = tick_data,
      aes(x = x, y = -1.45, label = label),
      inherit.aes = FALSE,
      family = font_family,
      colour = palette[["text"]],
      size = 2.65,
      vjust = 1
    ) +
    annotate(
      "text",
      x = (bubble_axis_left + bubble_axis_right) / 2,
      y = -5.0,
      label = "Gene ratio",
      family = font_family,
      colour = palette[["text"]],
      size = 2.90,
      vjust = 1
    ) +
    scale_size_continuous(
      name = "Overlap n",
      limits = c(9, 73),
      breaks = c(10, 30, 60),
      range = c(2.4, 7.8)
    ) +
    scale_fill_gradientn(
      name = "\u2212log10(BH q)",
      colours = q_colors,
      limits = c(1, 8),
      breaks = c(1, 4, 8),
      labels = c("1", "4", "\u22658"),
      oob = squish,
      na.value = "white"
    ) +
    scale_colour_manual(
      name = NULL,
      values = c(
        "BH q < 0.05" = palette[["text"]],
        "BH q \u2265 0.05" = palette[["nonsig"]]
      ),
      drop = FALSE
    ) +
    guides(
      size = guide_legend(
        order = 1,
        title.position = "top",
        ncol = if (legend_on_right) 1 else 3,
        override.aes = list(
          fill = "#C7D2FE",
          colour = palette[["text"]],
          shape = 21
        )
      ),
      fill = guide_colourbar(
        order = 2,
        display = "rectangles",
        title.position = "top",
        title.hjust = if (legend_on_right) 0 else 0.5,
        direction = if (legend_on_right) "vertical" else "horizontal",
        label.position = if (legend_on_right) "right" else "bottom",
        barwidth = grid::unit(
          if (legend_on_right) 2.8 else 28,
          "mm"
        ),
        barheight = grid::unit(
          if (legend_on_right) 28 else 2.6,
          "mm"
        ),
        ticks.colour = "white",
        frame.colour = NA
      ),
      colour = guide_legend(
        order = 3,
        title.position = "top",
        ncol = if (legend_on_right) 1 else 2,
        override.aes = list(
          fill = c("#818CF8", "white"),
          shape = 21,
          size = 3.8
        )
      )
    ) +
    coord_cartesian(
      xlim = plot_xlim,
      ylim = c(-7.2, axis_data$total_links),
      clip = "off",
      expand = FALSE
    ) +
    labs(title = title_text, subtitle = subtitle_text) +
    theme_void(base_family = font_family, base_size = 8) +
    theme(
      plot.title = element_text(
        family = font_family,
        face = "bold",
        colour = axis_colour,
        size = 10.5,
        hjust = 0,
        margin = margin(b = 1.5)
      ),
      plot.subtitle = element_text(
        family = font_family,
        colour = palette[["text_muted"]],
        size = 8.2,
        hjust = 0,
        margin = margin(b = 5)
      ),
      plot.margin = margin(5, 6, 2, 8),
      legend.position = if (legend_on_right) "right" else "bottom",
      legend.direction = if (legend_on_right) "vertical" else "horizontal",
      legend.box = if (legend_on_right) "vertical" else "horizontal",
      legend.box.just = if (legend_on_right) "left" else "center",
      legend.justification = "center",
      legend.title = element_text(
        family = font_family,
        colour = palette[["text"]],
        size = 8.5,
        face = "bold"
      ),
      legend.text = element_text(
        family = font_family,
        colour = palette[["text_muted"]],
        size = 8.0
      ),
      legend.key = element_blank(),
      legend.spacing.x = grid::unit(4, "pt"),
      legend.spacing.y = grid::unit(6, "pt"),
      legend.margin = if (legend_on_right) {
        margin(t = 2, r = 2, b = 2, l = 8)
      } else {
        margin(t = 1, r = 2, b = 0, l = 2)
      }
    )

  list(plot = p, data = axis_data, bubble_data = term_plot_data)
}

